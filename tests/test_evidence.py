"""SAM 证据核查测试，仅依赖标准库：
PYTHONPATH=/tmp/pylibs python3 -I -m unittest -v tests.test_evidence
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app.bio import CheckError, parse_fasta  # noqa: E402
from app.evidence import (  # noqa: E402
    check_evidence,
    evaluate,
    parse_sam,
)

# 1-based: 1A 2C 3G 4T | 5A 6A 7A 8A | 9C 10G 11T 12T
SEQ = "ACGTAAAACGTT"
FASTA = b">chrT\nACGTAAAACGTT\n"
HEADER = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
)


def vcf(rows):
    return (HEADER + "".join(r + "\n" for r in rows)).encode()


def sam(rows, headers=()):
    out = list(headers)
    out.extend(rows)
    return ("\n".join(out) + "\n").encode()


# Q20 的 Phred 质量字符为 '5'（33+20）；Q19 为 '4'
Q20 = "5"
Q19 = "4"


def sam_row(qname, flag, rname, pos, mapq, cigar, seq, qual=None, rnext="*"):
    if qual is None:
        qual = Q20 * len(seq)
    return (
        f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t0\t0\t"
        f"{seq}\t{qual}"
    )


def run(fasta=FASTA, vcf_bytes=None, sam_bytes=None):
    return check_evidence(fasta, vcf_bytes, sam_bytes)


class SamFilterTests(unittest.TestCase):
    def setUp(self):
        self.ref = parse_fasta(FASTA)

    def test_unmapped_secondary_supplementary_duplicate_ignored(self):
        rows = [
            sam_row("u", 0x4, "chrT", 2, 60, "3M", "AAA"),       # 未比对
            sam_row("s", 0x100, "chrT", 2, 60, "3M", "AAA"),     # 次级
            sam_row("x", 0x800, "chrT", 2, 60, "3M", "AAA"),     # 补充
            sam_row("d", 0x400, "chrT", 2, 60, "3M", "AAA"),     # PCR 重复
            sam_row("ok", 0, "chrT", 2, 60, "3M", "GGG"),        # 采纳
        ]
        grouped, ignored, used = parse_sam(sam(rows), self.ref)
        self.assertEqual(ignored, 4)
        self.assertEqual(used, 1)
        self.assertEqual(set(grouped), {"ok"})

    def test_mapq_below_20_ignored_20_kept(self):
        rows = [
            sam_row("lo", 0, "chrT", 2, 19, "3M", "GGG"),
            sam_row("hi", 0, "chrT", 2, 20, "3M", "GGG"),
        ]
        grouped, ignored, used = parse_sam(sam(rows), self.ref)
        self.assertEqual(ignored, 1)
        self.assertEqual(set(grouped), {"hi"})

    def test_base_quality_below_20_base_dropped(self):
        # pos2 参考 C；读段在 pos2 的碱基质量 Q19，应被剔除
        rows = [sam_row("r", 0, "chrT", 1, 60, "3M", "GGG", Q20 + Q19 + Q20)]
        grouped, _i, _u = parse_sam(sam(rows), self.ref)
        calls = grouped["r"][0].calls
        self.assertNotIn(1, calls)          # 0-based pos1 = 1-based pos2 被剔除
        self.assertEqual(calls[0], "G")
        self.assertEqual(calls[2], "G")

    def test_non_acgt_base_ignored_but_segment_still_counts(self):
        # N 碱基不采纳；同片段其他位点仍可用
        rows = [sam_row("r", 0, "chrT", 1, 60, "3M", "GNG")]
        grouped, _i, _u = parse_sam(sam(rows), self.ref)
        calls = grouped["r"][0].calls
        self.assertEqual(calls, {0: "G", 2: "G"})

    def test_rname_mismatch_fails(self):
        rows = [sam_row("r", 0, "chrOther", 2, 60, "3M", "GGG")]
        with self.assertRaises(CheckError):
            parse_sam(sam(rows), self.ref)

    def test_bad_cigar_operator_fails(self):
        for cigar in ("3N", "2X", "abc", "3M1", ""):
            with self.subTest(cigar=cigar):
                rows = [sam_row("r", 0, "chrT", 2, 60, cigar, "GGG")]
                with self.assertRaises(CheckError):
                    parse_sam(sam(rows), self.ref)

    def test_cigar_mids_accepted_and_soft_clip_skipped(self):
        # 2S3M：前两个查询碱基被软剪切，pos2-4 对应查询第 3-5 个碱基
        rows = [sam_row("r", 0, "chrT", 2, 60, "2S3M", "TTGGG")]
        grouped, _i, used = parse_sam(sam(rows), self.ref)
        self.assertEqual(used, 1)
        calls = grouped["r"][0].calls
        self.assertEqual(calls, {1: "G", 2: "G", 3: "G"})  # 0-based 1,2,3

    def test_insertion_consumes_query_only(self):
        # 1M2I2M：参考 pos1 一个碱基，插入 2 bp，再覆盖 pos2-3
        rows = [sam_row("r", 0, "chrT", 1, 60, "1M2I2M", "GTTAC")]
        grouped, _i, _u = parse_sam(sam(rows), self.ref)
        calls = grouped["r"][0].calls
        self.assertEqual(calls, {0: "G", 1: "A", 2: "C"})

    def test_deletion_marks_reference_positions(self):
        # 1M2D1M：pos1 覆盖，pos2-3 为缺失，pos4 覆盖
        rows = [sam_row("r", 0, "chrT", 1, 60, "1M2D1M", "GC", Q20 * 2)]
        grouped, _i, _u = parse_sam(sam(rows), self.ref)
        read = grouped["r"][0]
        self.assertEqual(read.calls, {0: "G", 3: "C"})
        self.assertEqual(read.deleted, {1, 2})

    def test_cigar_length_mismatch_fails(self):
        rows = [sam_row("r", 0, "chrT", 1, 60, "4M", "GC")]
        with self.assertRaises(CheckError):
            parse_sam(sam(rows), self.ref)

    def test_short_line_and_missing_seq_fail(self):
        with self.assertRaises(CheckError):
            parse_sam(sam(["r\t0\tchrT\t1\t60\t1M"]), self.ref)
        rows = [sam_row("r", 0, "chrT", 1, 60, "1M", "*", "*")]
        with self.assertRaises(CheckError):
            parse_sam(sam(rows), self.ref)

    def test_three_primary_reads_same_qname_fails(self):
        rows = [
            sam_row("p", 0x40, "chrT", 1, 60, "1M", "G"),
            sam_row("p", 0x80, "chrT", 2, 60, "1M", "G"),
            sam_row("p", 0, "chrT", 3, 60, "1M", "G"),
        ]
        with self.assertRaises(CheckError):
            parse_sam(sam(rows), self.ref)

    def test_same_end_flag_twice_fails(self):
        rows = [
            sam_row("p", 0x40, "chrT", 1, 60, "1M", "G"),
            sam_row("p", 0x40, "chrT", 2, 60, "1M", "G"),
        ]
        with self.assertRaises(CheckError):
            parse_sam(sam(rows), self.ref)

    def test_headers_skipped(self):
        rows = [sam_row("r", 0, "chrT", 2, 60, "1M", "G")]
        grouped, ignored, used = parse_sam(
            sam(rows, headers=["@HD\tVN:1.6", "@SQ\tSN:chrT\tLN:12"]), self.ref
        )
        self.assertEqual(used, 1)
        self.assertEqual(ignored, 0)


class FragmentMergeTests(unittest.TestCase):
    def setUp(self):
        self.ref = parse_fasta(FASTA)
        # 位点 2（1-based，参考 C）的 SNV C>G
        self.sites_vcf = vcf(["chrT\t2\t.\tC\tG\t.\t.\t.\tGT\t0/1"])

    def _evidence(self, rows):
        grouped, _i, _u = parse_sam(sam(rows), self.ref)
        from app.bio import parse_snv_sites

        sites, _sample = parse_snv_sites(self.sites_vcf, "x.vcf", self.ref)
        result = evaluate(self.ref, sites, grouped)
        self.assertEqual(len(result), 1)
        return result[0]

    def test_paired_ends_merged_one_fragment_same_base(self):
        rows = [
            sam_row("p", 0x40, "chrT", 1, 60, "2M", "AG"),   # pos2=G
            sam_row("p", 0x80, "chrT", 2, 60, "2M", "GA"),   # pos2=G
        ]
        ev = self._evidence(rows)
        self.assertEqual(ev.counts["G"], 1)   # 两端只算一个片段
        self.assertEqual(ev.coverage, 1)
        self.assertEqual(ev.conflicts, 0)

    def test_paired_end_conflict_counts_not_coverage(self):
        rows = [
            sam_row("p", 0x40, "chrT", 1, 60, "2M", "AG"),   # pos2=G
            sam_row("p", 0x80, "chrT", 2, 60, "2M", "TA"),   # pos2=T
        ]
        ev = self._evidence(rows)
        self.assertEqual(ev.conflicts, 1)
        self.assertEqual(ev.coverage, 0)
        self.assertEqual(ev.counts, {"A": 0, "C": 0, "G": 0, "T": 0})

    def test_deletion_fragment_excluded_from_coverage(self):
        rows = [
            # 片段 p：R1 覆盖 pos2 为 G；R2 在 pos2 落在 D 上 -> 整个片段不计覆盖
            sam_row("p", 0x40, "chrT", 1, 60, "2M", "AG"),
            sam_row("p", 0x80, "chrT", 1, 60, "1M1D1M", "AC", Q20 * 2),
        ]
        ev = self._evidence(rows)
        self.assertEqual(ev.deletion_fragments, 1)
        self.assertEqual(ev.coverage, 0)

    def test_one_end_overlapping_site_still_one_fragment(self):
        # 只有 R1 覆盖位点，R2 不覆盖：仍按一个片段计
        rows = [
            sam_row("p", 0x40, "chrT", 1, 60, "2M", "AG"),
            sam_row("p", 0x80, "chrT", 9, 60, "2M", "CG"),
        ]
        ev = self._evidence(rows)
        self.assertEqual(ev.counts["G"], 1)
        self.assertEqual(ev.coverage, 1)


class DecisionTests(unittest.TestCase):
    REF = parse_fasta(FASTA)

    def _report(self, vcf_rows, sam_rows):
        return run(
            vcf_bytes=vcf(vcf_rows),
            sam_bytes=sam(sam_rows),
        )

    def test_insufficient_when_coverage_below_10(self):
        v = ["chrT\t2\t.\tC\tG\t.\t.\t.\tGT\t0/1"]
        rows = [
            sam_row(f"g{i}", 0, "chrT", 1, 60, "2M", "AG") for i in range(9)
        ]
        report = self._report(v, rows)
        ev = report.sites[0]
        self.assertEqual(ev.coverage, 9)
        self.assertEqual(ev.status, "insufficient")

    def test_heterozygous_match_at_20_percent_boundary(self):
        # 10 个片段：G=2（20%）、C=8（80%）；GT 0/1 含 C、G，均 >=20%；A、T 为 0 <20%
        v = ["chrT\t2\t.\tC\tG\t.\t.\t.\tGT\t0/1"]
        rows = [sam_row(f"g{i}", 0, "chrT", 1, 60, "2M", "AG") for i in range(2)]
        rows += [sam_row(f"c{i}", 0, "chrT", 1, 60, "2M", "AC") for i in range(8)]
        ev = self._report(v, rows).sites[0]
        self.assertEqual(ev.coverage, 10)
        self.assertEqual(ev.counts, {"A": 0, "C": 8, "G": 2, "T": 0})
        self.assertEqual(ev.status, "match")

    def test_non_gt_base_at_20_percent_is_mismatch(self):
        # T 占 2/10 = 20%，T 不在 GT {C,G} 中 -> 不符
        v = ["chrT\t2\t.\tC\tG\t.\t.\t.\tGT\t0/1"]
        rows = [sam_row(f"g{i}", 0, "chrT", 1, 60, "2M", "AG") for i in range(2)]
        rows += [sam_row(f"c{i}", 0, "chrT", 1, 60, "2M", "AC") for i in range(6)]
        rows += [sam_row(f"t{i}", 0, "chrT", 1, 60, "2M", "AT") for i in range(2)]
        ev = self._report(v, rows).sites[0]
        self.assertEqual(ev.status, "mismatch")
        self.assertTrue(any("非 GT 碱基 T" in r for r in ev.mismatch_reasons))

    def test_homozygous_alt_other_base_below_20_ok(self):
        # GT 1/1 仅含 G：G 90% >=20%，C 10% <20% -> 相符
        v = ["chrT\t2\t.\tC\tG\t.\t.\t.\tGT\t1/1"]
        rows = [sam_row(f"g{i}", 0, "chrT", 1, 60, "2M", "AG") for i in range(9)]
        rows += [sam_row("c0", 0, "chrT", 1, 60, "2M", "AC")]
        ev = self._report(v, rows).sites[0]
        self.assertEqual(ev.gt_bases, ["G"])
        self.assertEqual(ev.alt_dosage, {"G": 2})
        self.assertEqual(ev.status, "match")

    def test_homozygous_gt_with_ref_base_at_20_percent_mismatch(self):
        # GT 1/1 只含 G，但 C 占 20% -> 不符
        v = ["chrT\t2\t.\tC\tG\t.\t.\t.\tGT\t1/1"]
        rows = [sam_row(f"g{i}", 0, "chrT", 1, 60, "2M", "AG") for i in range(8)]
        rows += [sam_row(f"c{i}", 0, "chrT", 1, 60, "2M", "AC") for i in range(2)]
        ev = self._report(v, rows).sites[0]
        self.assertEqual(ev.status, "mismatch")

    def test_multi_alt_shared_denominator_and_dosage(self):
        # 多 ALT 1/2：C(ref)/G/T，两个 ALT 共用同一片段分母，剂量各 1、总和 2
        v = ["chrT\t2\t.\tC\tG,T\t.\t.\t.\tGT\t1/2"]
        rows = [sam_row(f"g{i}", 0, "chrT", 1, 60, "2M", "AG") for i in range(5)]
        rows += [sam_row(f"t{i}", 0, "chrT", 1, 60, "2M", "AT") for i in range(5)]
        ev = self._report(v, rows).sites[0]
        self.assertEqual(ev.alt_dosage, {"G": 1, "T": 1})
        self.assertEqual(sum(ev.alt_dosage.values()), 2)
        self.assertEqual(ev.gt_bases, ["G", "T"])
        self.assertEqual(ev.coverage, 10)
        self.assertEqual(ev.status, "match")  # G/T 各 50%；C、A 为 0

    def test_conflicts_excluded_before_fraction_decision(self):
        # 10 个 G 片段 + 3 个冲突片段：覆盖 10 且全 G；GT 1/1 含 G -> 相符
        v = ["chrT\t2\t.\tC\tG\t.\t.\t.\tGT\t1/1"]
        rows = [sam_row(f"g{i}", 0, "chrT", 1, 60, "2M", "AG") for i in range(10)]
        for i in range(3):
            rows += [
                sam_row(f"x{i}", 0x40, "chrT", 1, 60, "2M", "AG"),
                sam_row(f"x{i}", 0x80, "chrT", 2, 60, "2M", "TA"),
            ]
        ev = self._report(v, rows).sites[0]
        self.assertEqual(ev.conflicts, 3)
        self.assertEqual(ev.coverage, 10)
        self.assertEqual(ev.status, "match")

    def test_line_numbers_reported(self):
        v = ["chrT\t2\t.\tC\tG\t.\t.\t.\tGT\t0/1"]
        rows = [sam_row(f"g{i}", 0, "chrT", 1, 60, "2M", "AG") for i in range(10)]
        ev = self._report(v, rows).sites[0]
        self.assertEqual(ev.line_numbers, [3])
        self.assertEqual(ev.raw_descriptions, ["2:C>G"])

    def test_two_sites_sorted_and_counted_independently(self):
        v = [
            "chrT\t2\t.\tC\tG\t.\t.\t.\tGT\t0/1",
            "chrT\t3\t.\tG\tT\t.\t.\t.\tGT\t0/1",
        ]
        # 每条读段 3M 同时覆盖 pos2 与 pos3：10 个片段
        rows = [sam_row(f"r{i}", 0, "chrT", 1, 60, "3M", "AGT") for i in range(10)]
        report = self._report(v, rows)
        self.assertEqual([s.pos for s in report.sites], [2, 3])
        self.assertTrue(all(s.coverage == 10 for s in report.sites))


class VcfRestrictionTests(unittest.TestCase):
    def test_called_insertion_fails(self):
        v = ["chrT\t4\t.\tT\tTA\t.\t.\t.\tGT\t0/1"]
        rows = [sam_row("r", 0, "chrT", 1, 60, "4M", "ACGT")]
        with self.assertRaises(CheckError) as ctx:
            run(vcf_bytes=vcf(v), sam_bytes=sam(rows))
        self.assertIn("单碱基替换", str(ctx.exception))

    def test_called_deletion_fails(self):
        v = ["chrT\t4\t.\tTA\tT\t.\t.\t.\tGT\t0/1"]
        rows = [sam_row("r", 0, "chrT", 1, 60, "2M", "AC")]
        with self.assertRaises(CheckError):
            run(vcf_bytes=vcf(v), sam_bytes=sam(rows))

    def test_zero_zero_site_not_analyzed(self):
        v = ["chrT\t2\t.\tC\tG\t.\t.\t.\tGT\t0/0"]
        rows = [sam_row("r", 0, "chrT", 1, 60, "2M", "AG")]
        report = run(vcf_bytes=vcf(v), sam_bytes=sam(rows))
        self.assertEqual(report.sites, [])

    def test_ref_mismatch_still_fails(self):
        v = ["chrT\t2\t.\tG\tT\t.\t.\t.\tGT\t0/1"]  # pos2 参考实为 C
        rows = [sam_row("r", 0, "chrT", 1, 60, "2M", "AG")]
        with self.assertRaises(CheckError):
            run(vcf_bytes=vcf(v), sam_bytes=sam(rows))

    def test_dosage_over_two_fails(self):
        # 两条记录在同位点累计剂量 3
        v = [
            "chrT\t2\t.\tC\tG\t.\t.\t.\tGT\t1/1",
            "chrT\t2\t.\tC\tG,T\t.\t.\t.\tGT\t0/2",
        ]
        rows = [sam_row("r", 0, "chrT", 1, 60, "2M", "AG")]
        with self.assertRaises(CheckError):
            run(vcf_bytes=vcf(v), sam_bytes=sam(rows))

    def test_merged_heterozygous_records_dosage_two(self):
        # 同位点两条 0/1（不同 ALT 碱基）合并：剂量 1+1=2，合法
        v = [
            "chrT\t2\t.\tC\tG\t.\t.\t.\tGT\t0/1",
            "chrT\t2\t.\tC\tT\t.\t.\t.\tGT\t0/1",
        ]
        rows = [sam_row(f"r{i}", 0, "chrT", 1, 60, "2M", "AG") for i in range(10)]
        report = run(vcf_bytes=vcf(v), sam_bytes=sam(rows))
        ev = report.sites[0]
        self.assertEqual(ev.alt_dosage, {"G": 1, "T": 1})
        self.assertEqual(ev.line_numbers, [3, 4])


class IntegrationTests(unittest.TestCase):
    def test_report_metadata(self):
        v = ["chrT\t2\t.\tC\tG\t.\t.\t.\tGT\t0/1"]
        rows = [
            sam_row("p", 0x40, "chrT", 1, 60, "2M", "AG"),
            sam_row("p", 0x80, "chrT", 2, 60, "2M", "GA"),
            sam_row("dup", 0x400, "chrT", 1, 60, "2M", "AC"),
        ]
        report = run(vcf_bytes=vcf(v), sam_bytes=sam(rows))
        self.assertEqual(report.ref_name, "chrT")
        self.assertEqual(report.ref_length, 12)
        self.assertEqual(report.sample_name, "S1")
        self.assertEqual(report.fragments_total, 1)
        self.assertEqual(report.reads_used, 2)
        self.assertEqual(report.reads_ignored, 1)

    def test_empty_sam_is_valid_zero_coverage(self):
        report = run(
            vcf_bytes=vcf(["chrT\t2\t.\tC\tG\t.\t.\t.\tGT\t0/1"]),
            sam_bytes=sam([]),
        )
        self.assertEqual(report.sites[0].coverage, 0)
        self.assertEqual(report.sites[0].status, "insufficient")


if __name__ == "__main__":
    unittest.main()
