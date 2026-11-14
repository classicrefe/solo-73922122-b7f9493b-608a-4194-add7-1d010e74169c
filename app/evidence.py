"""SAM 读段证据核查（仅依赖标准库）。

流程：
1. 复用 bio.parse_fasta / bio.parse_snv_sites 解析参考与 VCF，
   只核查被调用的单碱基替换位点；被调用的插缺在 VCF 解析阶段即报错。
2. 解析同染色体 SAM（CIGAR 仅接受 M/I/D/S）：
   - 忽略未比对（0x4）、次级（0x100）、补充（0x800）、PCR 重复（0x400）读段；
   - 整段映射质量 MAPQ < 20 的读段不采纳；
   - 只采纳碱基质量 Phred >= 20 且为 A/C/G/T 的碱基。
3. 按 QNAME 把两端合为一个片段：
   - 同位点两端碱基冲突 -> 记一次冲突，不计覆盖；
   - 片段在该位点落在 D（缺失）上 -> 不计覆盖；
   - 其余情况每个片段在该位点贡献一个碱基。
4. 按位点统计 A/C/G/T 片段数与冲突数，结合 GT 判定是否相符：
   - 有效覆盖 < 10 -> 证据不足；
   - 否则 GT 含有的每种碱基占比 >= 20%，其余碱基占比严格 < 20% 才算相符。

任何输入格式错误都抛出 CheckError，代表整次核查失败。
"""

import re
from dataclasses import dataclass, field

from .bio import CheckError, Reference, SiteCall, parse_fasta, parse_snv_sites

# 上传体积上限（字节）
MAX_SAM_BYTES = 256 * 1024

# 质量阈值
MIN_MAPQ = 20
MIN_BASEQ = 20

# 判定阈值
MIN_COVERAGE = 10
MIN_FRACTION = 0.20

_BASES = ("A", "C", "G", "T")

# SAM FLAG 位
FLAG_REVERSE = 0x10
FLAG_UNMAPPED = 0x4
FLAG_PCR_DUP = 0x400
FLAG_SECONDARY = 0x100
FLAG_SUPPLEMENTARY = 0x800
FLAG_READ1 = 0x40
FLAG_READ2 = 0x80
_IGNORE_FLAGS = FLAG_UNMAPPED | FLAG_SECONDARY | FLAG_SUPPLEMENTARY | FLAG_PCR_DUP

_CIGAR_RE = re.compile(r"(\d+)([MIDS])")
_CIGAR_FULL_RE = re.compile(r"(?:\d+[MIDS])+")


@dataclass
class AlignedRead:
    """一条通过过滤的 SAM 读段，坐标均为参考 0-based。"""

    qname: str
    end: int | None  # 1=READ1, 2=READ2, None=未标记
    calls: dict[int, str] = field(default_factory=dict)   # 参考位置 -> 观测碱基
    deleted: set[int] = field(default_factory=set)        # D 操作覆盖的参考位置


@dataclass
class SiteEvidence:
    """单个 VCF SNV 位点的片段级证据与判定。"""

    pos: int
    ref_base: str
    alt_bases: list[str]
    alt_dosage: dict[str, int]      # ALT 碱基 -> GT 剂量
    gt_bases: list[str]
    counts: dict[str, int]          # A/C/G/T 各自的片段数
    conflicts: int
    deletion_fragments: int         # 落在缺失上、未计覆盖的片段数
    line_numbers: list[int]
    raw_descriptions: list[str]

    @property
    def coverage(self) -> int:
        """有效覆盖：四种碱基片段数之和（冲突与缺失片段均不计）。"""
        return sum(self.counts.values())

    @property
    def fractions(self) -> dict[str, float]:
        cov = self.coverage
        return {b: (self.counts[b] / cov if cov else 0.0) for b in _BASES}

    @property
    def status(self) -> str:
        """insufficient（证据不足）| match（相符）| mismatch（不符）。"""
        if self.coverage < MIN_COVERAGE:
            return "insufficient"
        return "match" if not self.mismatch_reasons else "mismatch"

    @property
    def mismatch_reasons(self) -> list[str]:
        """覆盖达标时给出所有不满足 20% 规则的具体原因。"""
        cov = self.coverage
        if cov < MIN_COVERAGE:
            return []
        reasons: list[str] = []
        for b in _BASES:
            pct = self.counts[b] / cov
            if b in self.gt_bases:
                if pct < MIN_FRACTION:
                    reasons.append(f"GT 碱基 {b} 占比 {pct:.1%} 低于 20% 下限")
            elif pct >= MIN_FRACTION:
                reasons.append(f"非 GT 碱基 {b} 占比 {pct:.1%}，未严格低于 20%")
        return reasons


@dataclass
class EvidenceReport:
    ref_name: str
    ref_length: int
    sample_name: str
    fragments_total: int
    reads_ignored: int
    reads_used: int
    sites: list[SiteEvidence]


# --------------------------------------------------------------------------- #
# SAM 解析
# --------------------------------------------------------------------------- #
def _decode(token: str, what: str, line_no: int) -> str:
    # parse_sam 已将整份 SAM 解码为 str，这里保留行号上下文的包装函数以便扩展
    if not token:
        raise CheckError(f"SAM 第 {line_no} 行：{what} 为空。")
    return token


def _walk_cigar(
    cigar: str,
    seq: str,
    qual: str,
    start0: int,
    line_no: int,
) -> tuple[dict[int, str], set[int]]:
    """按 CIGAR 把查询碱基铺到参考坐标，返回 (位置->碱基, 缺失位置集合)。"""
    if not cigar or cigar == "*" or not _CIGAR_FULL_RE.fullmatch(cigar):
        raise CheckError(
            f"SAM 第 {line_no} 行：CIGAR {cigar!r} 非法，本工具仅接受 M/I/D/S 操作。"
        )

    calls: dict[int, str] = {}
    deleted: set[int] = set()
    q = 0
    r = start0
    for length_s, op in _CIGAR_RE.findall(cigar):
        length = int(length_s)
        if length == 0:
            raise CheckError(f"SAM 第 {line_no} 行：CIGAR 存在长度为 0 的操作。")
        if op == "M":
            if q + length > len(seq):
                raise CheckError(
                    f"SAM 第 {line_no} 行：CIGAR 超出 SEQ 长度（{len(seq)}）。"
                )
            for i in range(length):
                base = seq[q + i].upper()
                if base in _BASES and ord(qual[q + i]) - 33 >= MIN_BASEQ:
                    calls[r + i] = base
            q += length
            r += length
        elif op == "I":
            if q + length > len(seq):
                raise CheckError(
                    f"SAM 第 {line_no} 行：CIGAR 的 I 操作超出 SEQ 长度。"
                )
            q += length  # 插入不占参考坐标
        elif op == "D":
            deleted.update(range(r, r + length))
            r += length
        else:  # S：软剪切只消耗查询碱基
            if q + length > len(seq):
                raise CheckError(
                    f"SAM 第 {line_no} 行：CIGAR 的 S 操作超出 SEQ 长度。"
                )
            q += length

    if q != len(seq):
        raise CheckError(
            f"SAM 第 {line_no} 行：CIGAR 消耗的查询长度 {q} 与 SEQ 长度 {len(seq)} 不一致。"
        )
    return calls, deleted


def parse_sam(data: bytes, refobj: Reference) -> tuple[dict[str, list[AlignedRead]], int, int]:
    """解析 SAM，返回 (QNAME -> 读段列表, 被忽略读段数, 采纳读段数)。"""
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise CheckError(f"SAM 不是 UTF-8 文本（字节 {exc.start} 处解码失败）。")

    grouped: dict[str, list[AlignedRead]] = {}
    ignored = 0
    used = 0

    for line_no, raw in enumerate(text.splitlines(), start=1):
        line = raw.rstrip("\r")
        if not line or line.startswith("@"):
            continue

        fields = line.split("\t")
        if len(fields) < 11:
            raise CheckError(
                f"SAM 第 {line_no} 行：字段数 {len(fields)} 不足 11 列，不是合法 SAM 记录。"
            )

        qname = _decode(fields[0], "QNAME", line_no)
        try:
            flag = int(fields[1])
            pos = int(fields[3])
            mapq = int(fields[4])
        except ValueError:
            raise CheckError(f"SAM 第 {line_no} 行：FLAG/POS/MAPQ 存在非整数字段。")

        if flag & _IGNORE_FLAGS:
            # 未比对、次级、补充及 PCR 重复读段整体忽略
            ignored += 1
            continue

        rname = _decode(fields[2], "RNAME", line_no)
        if rname != refobj.name:
            raise CheckError(
                f"SAM 第 {line_no} 行：RNAME={rname!r} 与参考 FASTA 名称 "
                f"{refobj.name!r} 不一致，本工具仅接受同染色体 SAM。"
            )
        if pos < 1:
            raise CheckError(f"SAM 第 {line_no} 行：已比对读段 POS={pos} 非法（须 >= 1）。")

        if mapq < MIN_MAPQ:
            ignored += 1  # 映射质量不足 20，整段不采纳
            continue

        cigar = _decode(fields[5], "CIGAR", line_no)
        seq = _decode(fields[9], "SEQ", line_no)
        qual = _decode(fields[10], "QUAL", line_no)
        if seq == "*" or qual == "*":
            raise CheckError(
                f"SAM 第 {line_no} 行：已比对读段缺少 SEQ/QUAL，无法核查碱基证据。"
            )
        if len(seq) != len(qual):
            raise CheckError(
                f"SAM 第 {line_no} 行：SEQ 长度 {len(seq)} 与 QUAL 长度 {len(qual)} 不一致。"
            )

        # 反向比对读段的 SEQ/QUAL 在 SAM 中已按参考方向反向互补，可直接铺坐标
        calls, deleted = _walk_cigar(cigar, seq, qual, pos - 1, line_no)
        end = 1 if flag & FLAG_READ1 else 2 if flag & FLAG_READ2 else None
        grouped.setdefault(qname, []).append(
            AlignedRead(qname=qname, end=end, calls=calls, deleted=deleted)
        )
        used += 1

    # 一个 QNAME 至多对应一对读段，且不能同为 R1 或同为 R2
    for qname, reads in grouped.items():
        if len(reads) > 2:
            raise CheckError(
                f"SAM：QNAME {qname!r} 对应 {len(reads)} 段已采纳读段，"
                "双端片段至多包含两段。"
            )
        if len(reads) == 2 and reads[0].end is not None and reads[0].end == reads[1].end:
            raise CheckError(
                f"SAM：QNAME {qname!r} 的两段读段带有相同的 R{reads[0].end} 标记，"
                "无法合为一个双端片段。"
            )

    return grouped, ignored, used


# --------------------------------------------------------------------------- #
# 片段合并与位点判定
# --------------------------------------------------------------------------- #
def evaluate(
    refobj: Reference,
    sites: dict[int, SiteCall],
    fragments: dict[str, list[AlignedRead]],
) -> list[SiteEvidence]:
    """把片段证据汇总到每个 VCF SNV 位点。"""
    results: list[SiteEvidence] = []

    for pos in sorted(sites):
        site = sites[pos]
        pos0 = pos - 1
        counts = {b: 0 for b in _BASES}
        conflicts = 0
        deletion_fragments = 0

        for reads in fragments.values():
            if any(pos0 in read.deleted for read in reads):
                deletion_fragments += 1
                continue
            observed = {
                base for read in reads for p, base in read.calls.items() if p == pos0
            }
            if len(observed) > 1:
                conflicts += 1  # 两端冲突：记冲突，不计覆盖
            elif observed:
                counts[next(iter(observed))] += 1

        results.append(
            SiteEvidence(
                pos=pos,
                ref_base=site.ref_base,
                alt_bases=list(site.alt_dosage),
                alt_dosage=dict(site.alt_dosage),
                gt_bases=list(site.gt_bases),
                counts=counts,
                conflicts=conflicts,
                deletion_fragments=deletion_fragments,
                line_numbers=list(site.line_numbers),
                raw_descriptions=list(site.raw_descriptions),
            )
        )

    return results


def check_evidence(fasta_bytes: bytes, vcf_bytes: bytes, sam_bytes: bytes) -> EvidenceReport:
    """证据核查总入口：任何输入错误抛 CheckError（整次失败）。"""
    refobj = parse_fasta(fasta_bytes)
    sites, sample_name = parse_snv_sites(vcf_bytes, "待核查 VCF", refobj)
    fragments, ignored, used = parse_sam(sam_bytes, refobj)
    site_results = evaluate(refobj, sites, fragments)
    return EvidenceReport(
        ref_name=refobj.name,
        ref_length=len(refobj.sequence),
        sample_name=sample_name,
        fragments_total=len(fragments),
        reads_ignored=ignored,
        reads_used=used,
        sites=site_results,
    )
