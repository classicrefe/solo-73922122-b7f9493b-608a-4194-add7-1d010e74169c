"""从 samples/demo.fa 生成两份演示 VCF（坐标由参考序列自动计算，保证 REF 合法）。

运行：python3 scripts/make_demo.py
场景覆盖：
  SNV 一致 / 同一插入的不同写法 / 同一缺失的不同写法 / 多 ALT /
  剂量不符 / 仅单侧出现 / 剂量 0 被忽略

同时生成 samples/evidence_demo.sam 与 samples/evidence_demo.vcf：
用合成双端读段演示 GT 证据核查的相符、证据不足、两端冲突、缺失片段与
MAPQ/碱基质量/FLAG 过滤。
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app.bio import parse_fasta  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
HEADER = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
)


def find(seq: str, motif: str, start: int = 1) -> int:
    """返回 motif 在参考中的 1-based 起点。"""
    idx = seq.find(motif, start - 1)
    if idx < 0:
        raise SystemExit(f"motif {motif!r} not found in reference")
    return idx + 1


def row(pos, ref, alt, gt, extra_fmt=None):
    fmt = "GT" if extra_fmt is None else "GT:DP"
    sample = gt if extra_fmt is None else f"{gt}:{extra_fmt}"
    return f"chrDemo\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\t{fmt}\t{sample}"


def main() -> None:
    ref = parse_fasta((ROOT / "samples" / "demo.fa").read_bytes())
    seq = ref.sequence

    # 找一段纯 A 重复区，演示插入/缺失在重复区内的等价写法
    a_run = find(seq, "AAAA")
    # 找一个 AT 二碱基重复区，演示不同锚点写法的缺失
    at_run = find(seq, "ATATATAT")

    rows_a: list[str] = []
    rows_b: list[str] = []

    # 1) SNV 两侧一致（0/1）
    p = find(seq, "ACGTAC", 1)
    snv_pos = p + 1  # 参考碱基 C -> G
    rows_a.append(row(snv_pos, "C", "G", "0/1"))
    rows_b.append(row(snv_pos, "C", "G", "0/1"))

    # 2) 同一插入的两种重复区写法：A 在 A 串前插入 A（左对齐后应相同）
    anchor_pos = a_run - 1          # A 串前一个碱基的 1-based 位置
    anchor = seq[anchor_pos - 1]    # 该锚定碱基本身
    rows_a.append(row(anchor_pos, anchor, anchor + "A", "0/1"))
    rows_b.append(row(a_run, "A", "AA", "0/1"))

    # 3) 同一缺失的两种写法（AT 重复区，删一个 AT）
    rows_a.append(row(at_run - 1, seq[at_run - 2] + "AT", seq[at_run - 2], "0/1"))
    rows_b.append(row(at_run + 1, "TAT", "T", "0/1"))

    # 4) 多 ALT 记录：A 侧 1/2 两个等位各剂量 1；B 侧只报第一个 ALT
    p = find(seq, "GGATCC")
    multi_pos = p  # 参考 G
    rows_a.append(row(multi_pos, "G", "A,T", "1/2"))
    rows_b.append(row(multi_pos, "G", "A", "0/1"))

    # 5) 剂量不符：A 纯合 1/1，B 杂合 0/1
    p = find(seq, "TTTT")
    mm_pos = p + 1  # 参考 T
    rows_a.append(row(mm_pos, "T", "C", "1/1", "42"))
    rows_b.append(row(mm_pos, "T", "C", "0/1", "40"))

    # 6) 仅 A 侧出现
    p = find(seq, "GCGCGCGC")
    a_only_pos = p  # 参考 G
    rows_a.append(row(a_only_pos, "G", "T", "0/1"))

    # 7) 仅 B 侧出现
    b_only_pos = p + 2  # 同区域内另一个 G
    rows_b.append(row(b_only_pos, "G", "A", "0/1"))

    # 8) B 侧一条 0/0 记录：剂量 0，应被完全忽略
    rows_b.append(row(snv_pos, "C", "G", "0/0"))

    (ROOT / "samples" / "caller_a.vcf").write_text(
        HEADER + "\n".join(rows_a) + "\n", encoding="utf-8"
    )
    (ROOT / "samples" / "caller_b.vcf").write_text(
        HEADER + "\n".join(rows_b) + "\n", encoding="utf-8"
    )
    print(f"wrote {len(rows_a)} records -> samples/caller_a.vcf")
    print(f"wrote {len(rows_b)} records -> samples/caller_b.vcf")

    write_evidence_demo(ROOT, seq)


def write_evidence_demo(root: Path, seq: str) -> None:
    """生成证据核查演示 VCF + 同染色体 SAM。"""
    # 选三个互不重叠的 SNV 位点
    p = find(seq, "ACGTAC", 1)
    pos_match = p + 1      # 参考 C -> G，杂合 0/1，构造为相符
    pos_multi = p + 2      # 参考 G -> A,T，1/2，两个 ALT 共用片段分母
    pos_weak = find(seq, "GCGCGCGC")  # 参考 G -> T，覆盖不足

    vcf_rows = [
        f"chrDemo\t{pos_match}\t.\t{seq[pos_match-1]}\tG\t.\tPASS\t.\tGT\t0/1",
        f"chrDemo\t{pos_multi}\t.\t{seq[pos_multi-1]}\tA,T\t.\tPASS\t.\tGT\t1/2",
        f"chrDemo\t{pos_weak}\t.\t{seq[pos_weak-1]}\tT\t.\tPASS\t.\tGT\t0/1",
    ]
    (root / "samples" / "evidence_demo.vcf").write_text(
        HEADER + "\n".join(vcf_rows) + "\n", encoding="utf-8"
    )

    q20 = chr(33 + 20)  # '5'
    q19 = chr(33 + 19)  # '4'

    def sam_row(qname, flag, pos, cigar, read_seq, qual=None, mapq=60):
        if qual is None:
            qual = q20 * len(read_seq)
        return (
            f"{qname}\t{flag}\tchrDemo\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t"
            f"{read_seq}\t{qual}"
        )

    def window(target: int, span: int = 6):
        """返回覆盖 target（1-based）的读段起点与参考窗口序列。"""
        start = max(1, target - span // 2 + 1)
        return start, seq[start - 1:start - 1 + span]

    sam_lines = ["@HD\tVN:1.6", f"@SQ\tSN:chrDemo\tLN:{len(seq)}"]

    # 位点 1（相符）：12 个片段，参考碱基 C 与 ALT G 各占一半（均 >=20%，A/T 为 0）
    start, win = window(pos_match)
    idx = pos_match - start
    for i in range(6):
        s = win[:idx] + "C" + win[idx + 1:]
        sam_lines.append(sam_row(f"mC{i}", 0, start, f"{len(win)}M", s))
    for i in range(6):
        s = win[:idx] + "G" + win[idx + 1:]
        sam_lines.append(sam_row(f"mG{i}", 0, start, f"{len(win)}M", s))

    # 双端冲突片段：两端在位点 1 观测到不同碱基 -> 记冲突、不计覆盖
    sam_lines.append(sam_row(
        "conf", 0x40, start, f"{len(win)}M",
        win[:idx] + "G" + win[idx + 1:],
    ))
    sam_lines.append(sam_row(
        "conf", 0x80, start, f"{len(win)}M",
        win[:idx] + "T" + win[idx + 1:],
    ))

    # 落在缺失上的片段：位点 1 处于 CIGAR 的 D 区间 -> 不计覆盖
    d_anchor = pos_match - 1
    remaining = len(win) - 1
    sam_lines.append(sam_row(
        "delfrag", 0, d_anchor,
        f"1M1D{remaining}M",
        seq[d_anchor - 1] + seq[pos_match:pos_match + remaining],
    ))

    # 位点 2（多 ALT，相符）：10 个片段，A 与 T 各 5；GT 1/2 只含 A、T
    start2, win2 = window(pos_multi)
    idx2 = pos_multi - start2
    for i in range(5):
        s = win2[:idx2] + "A" + win2[idx2 + 1:]
        sam_lines.append(sam_row(f"mA{i}", 0, start2, f"{len(win2)}M", s))
    for i in range(5):
        s = win2[:idx2] + "T" + win2[idx2 + 1:]
        sam_lines.append(sam_row(f"mT{i}", 0, start2, f"{len(win2)}M", s))

    # 位点 3（证据不足）：仅 6 个片段，覆盖 < 10
    start3, win3 = window(pos_weak)
    idx3 = pos_weak - start3
    for i in range(6):
        s = win3[:idx3] + "T" + win3[idx3 + 1:]
        sam_lines.append(sam_row(f"w{i}", 0, start3, f"{len(win3)}M", s))

    # 过滤演示：PCR 重复、未比对、次级、补充、低 MAPQ、低碱基质量
    sam_lines.append(sam_row("dup", 0x400, start, f"{len(win)}M", win))
    sam_lines.append(sam_row("unmapped", 0x4, start, f"{len(win)}M", win))
    sam_lines.append(sam_row("secondary", 0x100, start, f"{len(win)}M", win))
    sam_lines.append(sam_row("supp", 0x800, start, f"{len(win)}M", win))
    sam_lines.append(sam_row("lowmapq", 0, start, f"{len(win)}M", win, mapq=19))
    sam_lines.append(sam_row(
        "lowbaseq", 0, start, f"{len(win)}M", win,
        qual=q19 * len(win),
    ))

    (root / "samples" / "evidence_demo.sam").write_text(
        "\n".join(sam_lines) + "\n", encoding="utf-8"
    )
    print("wrote samples/evidence_demo.vcf + samples/evidence_demo.sam")


if __name__ == "__main__":
    main()
