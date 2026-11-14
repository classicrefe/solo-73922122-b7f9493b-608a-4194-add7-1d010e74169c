"""Web 层端到端测试，需要 fastapi/httpx 依赖：
PYTHONPATH=/tmp/pylibs python -I -m unittest -v tests.test_web
"""

import io
import sys
import unittest
from pathlib import Path

# 本地无 pip 环境时，依赖安装在隔离目录；容器内该目录不存在则跳过
for _cand in (Path("/tmp/pylibs"), Path("/workspace/.pylibs")):
    if _cand.is_dir() and str(_cand) not in sys.path:
        sys.path.insert(0, str(_cand))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from fastapi.testclient import TestClient  # noqa: E402

from app.main import app  # noqa: E402
from app.bio import MAX_FASTA_BYTES, MAX_VCF_BYTES  # noqa: E402
from app.evidence import MAX_SAM_BYTES  # noqa: E402

client = TestClient(app)

FASTA = b">chrT\nACGTAAAACGTT\n"
HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"


def vcf(rows):
    return (HEADER + "".join(r + "\n" for r in rows)).encode()


def upload(fasta=FASTA, a=None, b=None):
    return client.post(
        "/compare",
        files={
            "fasta": ("ref.fa", fasta, "text/plain"),
            "vcf_a": ("a.vcf", a if a is not None else b"", "text/plain"),
            "vcf_b": ("b.vcf", b if b is not None else b"", "text/plain"),
        },
    )


class WebTests(unittest.TestCase):
    def test_index_renders(self):
        resp = client.get("/")
        self.assertEqual(resp.status_code, 200)
        self.assertIn("双 VCF 变异核对".encode(), resp.content)

    def test_happy_path_groups_and_line_numbers(self):
        a = vcf([
            "chrT\t2\t.\tC\tG\t.\t.\t.\tGT\t0/1",       # SNV 一致
            "chrT\t4\t.\tT\tTA\t.\t.\t.\tGT\t0/1",      # 插入（写法一）
            "chrT\t10\t.\tG\tT\t.\t.\t.\tGT\t1/1",      # 仅 A 有
        ])
        b = vcf([
            "chrT\t2\t.\tC\tG\t.\t.\t.\tGT\t0/1",
            "chrT\t5\t.\tA\tAA\t.\t.\t.\tGT\t0/1",      # 同一插入的另一种写法
            "chrT\t3\t.\tG\tC\t.\t.\t.\tGT\t0/1",       # 仅 B 有
        ])
        resp = upload(a=a, b=b)
        self.assertEqual(resp.status_code, 200)
        body = resp.text
        self.assertIn("两侧一致", body)
        self.assertIn("4:TinsTA", body)          # 归一变异
        self.assertIn("原始行号：4", body)        # 插入记录的物理行号
        self.assertIn("5:A&gt;AA", body)         # B 侧保留原始写法
        self.assertIn("仅单侧出现", body)
        self.assertIn("剂量不符", body)

    def test_ref_mismatch_shows_failure_page_without_partial_results(self):
        a = vcf(["chrT\t2\t.\tG\tT\t.\t.\t.\tGT\t0/1"])  # pos2 参考是 C
        b = vcf(["chrT\t2\t.\tC\tG\t.\t.\t.\tGT\t0/1"])
        resp = upload(a=a, b=b)
        self.assertEqual(resp.status_code, 422)
        self.assertIn("核对失败", resp.text)
        self.assertNotIn("核对概览", resp.text)

    def test_size_limit_enforced(self):
        big = b">chrT\n" + b"A" * (MAX_FASTA_BYTES + 10)
        resp = upload(fasta=big, a=vcf([]), b=vcf([]))
        self.assertEqual(resp.status_code, 422)
        self.assertIn("超过大小上限", resp.text)

        big_vcf = HEADER.encode() + b"#\n" * (MAX_VCF_BYTES // 2)
        resp = upload(a=big_vcf, b=vcf([]))
        self.assertEqual(resp.status_code, 422)

    def test_empty_file_rejected(self):
        resp = upload(a=b"", b=vcf([]))
        self.assertEqual(resp.status_code, 422)
        self.assertIn("为空", resp.text)

    def test_api_json_happy_path(self):
        a = vcf([
            "chrT\t2\t.\tC\tG,T\t.\t.\t.\tGT\t1/2",
            "chrT\t4\t.\tT\tTA\t.\t.\t.\tGT\t1/1",
        ])
        b = vcf([
            "chrT\t2\t.\tC\tG\t.\t.\t.\tGT\t0/1",  # G 剂量不符（1 vs 2）
            "chrT\t4\t.\tT\tTA\t.\t.\t.\tGT\t1/1",
        ])
        resp = client.post(
            "/api/compare",
            files={
                "fasta": ("ref.fa", FASTA, "text/plain"),
                "vcf_a": ("a.vcf", a, "text/plain"),
                "vcf_b": ("b.vcf", b, "text/plain"),
            },
        )
        self.assertEqual(resp.status_code, 200)
        data = resp.json()
        self.assertTrue(data["ok"])
        labels = {e["variant"]["label"] for e in data["result"]["consistent"]}
        self.assertIn("4:TinsTA", labels)
        # 多 ALT 记录 1/2：G、T 各剂量 1；B 侧只有 G 剂量 1 → T 仅 A 有，G 一致
        mismatch = data["result"]["dosage_mismatch"]
        self.assertEqual(mismatch, [])
        left_only = {e["variant"]["label"]: e for e in data["result"]["left_only"]}
        self.assertIn("2:C>T", left_only)
        self.assertEqual(left_only["2:C>T"]["left"]["dosage"], 1)

    def test_api_json_dosage_mismatch(self):
        a = vcf(["chrT\t2\t.\tC\tG\t.\t.\t.\tGT\t1/1"])
        b = vcf(["chrT\t2\t.\tC\tG\t.\t.\t.\tGT\t0/1"])
        resp = client.post(
            "/api/compare",
            files={
                "fasta": ("ref.fa", FASTA, "text/plain"),
                "vcf_a": ("a.vcf", a, "text/plain"),
                "vcf_b": ("b.vcf", b, "text/plain"),
            },
        )
        self.assertEqual(resp.status_code, 200)
        mismatch = resp.json()["result"]["dosage_mismatch"]
        self.assertEqual(len(mismatch), 1)
        self.assertEqual(mismatch[0]["left"]["dosage"], 2)
        self.assertEqual(mismatch[0]["right"]["dosage"], 1)

    def test_api_json_failure(self):
        resp = client.post(
            "/api/compare",
            files={
                "fasta": ("ref.fa", FASTA, "text/plain"),
                "vcf_a": ("a.vcf", vcf(["chrT\t2\t.\tG\tT\t.\t.\t.\tGT\t0/1"]), "text/plain"),
                "vcf_b": ("b.vcf", vcf([]), "text/plain"),
            },
        )
        self.assertEqual(resp.status_code, 422)
        self.assertFalse(resp.json()["ok"])


# --------------------------------------------------------------------------- #
# 读段证据核查
# --------------------------------------------------------------------------- #
EVIDENCE_VCF = vcf([
    "chrT\t2\t.\tC\tG\t.\t.\t.\tGT\t0/1",
    "chrT\t3\t.\tG\tT\t.\t.\t.\tGT\t1/1",
])

Q20 = "5"


def sam_line(qname, flag, pos, cigar, seq, qual=None):
    if qual is None:
        qual = Q20 * len(seq)
    return (
        f"{qname}\t{flag}\tchrT\t{pos}\t60\t{cigar}\t*\t0\t0\t{seq}\t{qual}"
    )


def evidence_sam(lines, headers=("@HD\tVN:1.6", "@SQ\tSN:chrT\tLN:12")):
    return ("\n".join(list(headers) + lines) + "\n").encode()


def upload_evidence(vcf_bytes=EVIDENCE_VCF, sam_bytes=None, fasta=FASTA):
    return client.post(
        "/evidence",
        files={
            "fasta": ("ref.fa", fasta, "text/plain"),
            "vcf": ("c.vcf", vcf_bytes, "text/plain"),
            "sam": ("r.sam", sam_bytes if sam_bytes is not None else b"", "text/plain"),
        },
    )


class EvidenceWebTests(unittest.TestCase):
    def test_evidence_entry_renders(self):
        resp = client.get("/evidence")
        self.assertEqual(resp.status_code, 200)
        self.assertIn("读段证据核查".encode(), resp.content)
        self.assertIn("256 KiB".encode(), resp.content)

    def test_index_links_to_evidence(self):
        resp = client.get("/")
        self.assertIn(b"/evidence", resp.content)

    def test_evidence_happy_path_page(self):
        # 位点 2：8 个 C + 2 个 G（C 80%、G 20%，GT 0/1 含 C/G -> 相符）
        # 位点 3：8 个读段给出参考 G、2 个给出 T（G 80%、T 20%，GT 0/1 -> 相符）
        lines = []
        for i in range(8):
            lines.append(sam_line(f"c{i}", 0, 1, "3M", "ACG"))  # pos2=C, pos3=G
        for i in range(2):
            lines.append(sam_line(f"g{i}", 0, 1, "3M", "AGT"))  # pos2=G, pos3=T
        resp = upload_evidence(
            vcf_bytes=vcf([
                "chrT\t2\t.\tC\tG\t.\t.\t.\tGT\t0/1",
                "chrT\t3\t.\tG\tT\t.\t.\t.\tGT\t0/1",
            ]),
            sam_bytes=evidence_sam(lines),
        )
        self.assertEqual(resp.status_code, 200)
        body = resp.text
        self.assertIn("位点证据明细", body)
        self.assertIn("有效覆盖", body)
        self.assertIn("原始行号", body)
        # 两个位点都应判相符：表格状态列出现两次“相符”徽标
        self.assertEqual(body.count('color:var(--ok)">相符'), 2)
        self.assertIn("3:G", body)

    def test_evidence_conflict_and_insufficient_shown(self):
        # 位点 2：一个双端冲突片段 + 8 个 G 片段 -> 覆盖 8 < 10，证据不足
        lines = [
            sam_line("x", 0x40, 1, "2M", "AG"),
            sam_line("x", 0x80, 2, "2M", "TA"),
        ]
        for i in range(8):
            lines.append(sam_line(f"g{i}", 0, 1, "2M", "AG"))
        resp = upload_evidence(
            vcf_bytes=vcf(["chrT\t2\t.\tC\tG\t.\t.\t.\tGT\t0/1"]),
            sam_bytes=evidence_sam(lines),
        )
        self.assertEqual(resp.status_code, 200)
        self.assertIn("证据不足", resp.text)
        self.assertIn(">1<", resp.text)  # 冲突数 1

    def test_evidence_indel_failure_is_whole_failure(self):
        bad_vcf = vcf(["chrT\t4\t.\tT\tTA\t.\t.\t.\tGT\t0/1"])
        resp = upload_evidence(
            vcf_bytes=bad_vcf,
            sam_bytes=evidence_sam([sam_line("r", 0, 1, "4M", "ACGT")]),
        )
        self.assertEqual(resp.status_code, 422)
        self.assertIn("证据核查失败", resp.text)
        self.assertNotIn("位点证据明细", resp.text)

    def test_evidence_sam_size_limit(self):
        big_sam = b"@HD\tVN:1.6\n" + b"#\n" * (MAX_SAM_BYTES // 2)
        resp = upload_evidence(sam_bytes=big_sam)
        self.assertEqual(resp.status_code, 422)
        self.assertIn("超过大小上限", resp.text)

    def test_evidence_empty_sam_file_rejected(self):
        resp = upload_evidence(sam_bytes=b"")
        self.assertEqual(resp.status_code, 422)
        self.assertIn("为空", resp.text)

    def test_evidence_bad_cigar_failure(self):
        resp = upload_evidence(
            sam_bytes=evidence_sam([sam_line("r", 0, 1, "3N", "ACG")])
        )
        self.assertEqual(resp.status_code, 422)
        self.assertIn("CIGAR", resp.text)

    def test_api_evidence_happy_path(self):
        # GT 0/1 含 C、G：8 G + 2 C，两种 GT 碱基均 >=20%
        lines = [sam_line(f"g{i}", 0, 1, "2M", "AG") for i in range(8)]
        lines += [sam_line(f"c{i}", 0, 1, "2M", "AC") for i in range(2)]
        resp = client.post(
            "/api/evidence",
            files={
                "fasta": ("ref.fa", FASTA, "text/plain"),
                "vcf": ("c.vcf", vcf(["chrT\t2\t.\tC\tG\t.\t.\t.\tGT\t0/1"]), "text/plain"),
                "sam": ("r.sam", evidence_sam(lines), "text/plain"),
            },
        )
        self.assertEqual(resp.status_code, 200)
        data = resp.json()
        self.assertTrue(data["ok"])
        site = data["sites"][0]
        self.assertEqual(site["pos"], 2)
        self.assertEqual(site["fragments"], {"A": 0, "C": 2, "G": 8, "T": 0})
        self.assertEqual(site["coverage"], 10)
        self.assertEqual(site["line_numbers"], [3])
        self.assertEqual(site["status"], "match")

    def test_api_evidence_failure_422(self):
        resp = client.post(
            "/api/evidence",
            files={
                "fasta": ("ref.fa", FASTA, "text/plain"),
                "vcf": ("c.vcf", vcf(["chrT\t4\t.\tT\tTA\t.\t.\t.\tGT\t0/1"]), "text/plain"),
                "sam": ("r.sam", evidence_sam([]), "text/plain"),
            },
        )
        self.assertEqual(resp.status_code, 422)
        self.assertFalse(resp.json()["ok"])


if __name__ == "__main__":
    unittest.main()
