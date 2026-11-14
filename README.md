# 双 VCF 变异核对（FastAPI）

同一参考、同一样本下，核对两套变异识别程序输出的单样本 VCF，按**归一化后的等位变异**
列出两侧一致、仅单侧出现与剂量不符的结果，并保留每条结论的 VCF 原始行号，供人工复查。

另提供**读段证据核查**入口：放行单样本 VCF 前，上传参考 FASTA、一份 VCF 与同染色体 SAM，
按 QNAME 合并双端读段，逐位点核查被调用单碱基替换的 GT 是否得到原始读段支持。

## 功能与判定规则

- **输入限制**
  - 单条染色体参考 FASTA，纯 A/C/G/T 序列，不超过 **1 MiB**；
  - 两份单样本 VCF（恰好 10 列、1 个样本），各不超过 **256 KiB**。
- **接受的变异**
  - A/C/G/T 构成的**单碱基替换**（SNV）；
  - 事件长度不超过 **20 bp** 的**纯插入 / 纯缺失**；
  - GT 必须是**非缺失、未分相的二倍体**（`0/1`、`1/1`，多 ALT 可为 `1/2`）。
  - 其他形态（MNV、符号化 ALT、`*`、分相 `|`、缺失 `.`、非二倍体等）一律判失败。
- **归一化**：对照参考序列对 indel 做左对齐（与 `bcftools norm` / `vt normalize` 等价），
  重复区域中位置写法不同的同一插入 / 缺失归为同一等位变异；SNV 位置不变。
- **剂量计算**：多 ALT 记录按 GT 分别计算每个 ALT 的 0/1/2 剂量，剂量为 0 的变异忽略；
  同一文件内同一归一变异的多条记录剂量累加，超过二倍体上限 2 即判定为**互相矛盾**。
- **整次失败**（返回 422，只显示错误、不展示任何局部结论）：
  - 任一记录 REF 与参考序列不符（含坐标越界、CHROM 不一致）；
  - 同一文件中同一归一变异存在互相矛盾的剂量；
  - 以及任何不满足上述输入规则的情形。
- **结果分组**（均按归一变异排序，并标注两侧原始行号与原始 `POS:REF>ALT` 写法）：
  两侧一致 / 仅 VCF A 或 B 出现 / 剂量不符。

## 读段证据核查（GT 放行前核查）

入口：首页“读段证据核查”卡片 → `/evidence`；JSON 接口 `POST /api/evidence`。

- **输入限制**：参考 FASTA 不超过 1 MiB；待核查单样本 VCF 与同染色体 SAM 各不超过
  **256 KiB**。沿用双 VCF 核对的文件大小、REF 与 GT 校验。
- **变异范围**：只分析**被调用（GT 剂量非 0）的单碱基替换**；任何被调用的插入 / 缺失
  （及其他非 SNV 形态）直接判整次失败；`0/0` 位点不产生核查任务。
- **SAM 过滤**：
  - CIGAR 仅接受 `M/I/D/S`，其他操作或 CIGAR 与 SEQ 长度不一致即整次失败；
  - 忽略未比对（FLAG 0x4）、次级（0x100）、补充（0x800）、PCR 重复（0x400）读段；
  - 整段 MAPQ 低于 20 的读段不采纳；仅采纳碱基质量 Phred ≥ 20 且为 A/C/G/T 的碱基。
- **片段合并**：按 **QNAME** 把两端合为一个片段。
  - 同一位点两端碱基冲突 → 记一次**冲突**，不计覆盖；
  - 片段在该位点落在 CIGAR 的 D（缺失）上 → 不计覆盖（单列“缺失片段”数）；
  - 其余情况每个片段在该位点贡献一个碱基。
- **位点统计**：按位点列出 A/C/G/T 四种碱基的**片段数**与占比、冲突数、缺失片段数及
  **原始 VCF 行号**；多个 ALT 共用同一片段分母，位点 ALT 剂量总和不得超过 2。
- **判定**：
  - 有效覆盖（四碱基片段数之和）**少于 10** → 证据不足；
  - 否则 GT 含有的每种碱基占比须 **≥ 20%**，其他碱基占比均**严格低于 20%**，
    满足判“相符”，否则判“不符”并列出具体原因。
- 任何输入错误均**整次失败**（422，只显示错误、不展示局部结论）；双 VCF 比对功能保持不变。

## 目录结构

```
app/            FastAPI 应用与核心逻辑（bio.py / evidence.py 仅用标准库）
  main.py       路由与上传校验（双 VCF 核对 + SAM 证据核查）
  bio.py        FASTA/VCF 解析、左对齐归一化、剂量与分组、SNV 位点提取
  evidence.py   SAM 解析、QNAME 片段合并、位点覆盖与 GT 相符判定
  templates/    上传页、双 VCF 结果页、证据核查页（Jinja2）
tests/          单元 / 端到端测试（93 个）
samples/        演示数据（demo.fa、两份 VCF、evidence_demo.vcf/.sam）
scripts/        演示数据生成器（坐标由参考自动计算）
Dockerfile, docker-compose.yml, requirements.txt
```

## 用 Docker Compose 启动（推荐）

```bash
docker compose up --build
```

浏览器入口：

- 上传页面：<http://localhost:8000/>
- 交互式 API 文档（Swagger UI）：<http://localhost:8000/docs>

停止：`docker compose down`。容器以只读根文件系统、`no-new-privileges` 运行，
内置 `/` 健康检查。

## 本地运行（无 Docker 时）

```bash
python3 -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt
uvicorn app.main:app --host 0.0.0.0 --port 8000
```

## 试用

打开 <http://localhost:8000/>，分别选择：

- 参考：`samples/demo.fa`
- VCF A：`samples/caller_a.vcf`
- VCF B：`samples/caller_b.vcf`

演示数据覆盖：SNV 一致、重复区内插入 / 缺失的两种等价写法归一、多 ALT 按 GT 拆分、
`0/0` 剂量 0 被忽略、剂量不符（1/1 vs 0/1）、仅单侧出现。

证据核查演示（打开 `/evidence`）：

- 参考：`samples/demo.fa`
- VCF：`samples/evidence_demo.vcf`
- SAM：`samples/evidence_demo.sam`

该 SAM 覆盖三种判定（相符 / 不符 / 证据不足），并含一对冲突双端读段、一个落在缺失上的
片段，以及 PCR 重复、未比对、次级、补充、低 MAPQ、低碱基质量等应被忽略的读段。

也可以走 JSON 接口：

```bash
curl -F 'fasta=@samples/demo.fa' \
     -F 'vcf_a=@samples/caller_a.vcf' \
     -F 'vcf_b=@samples/caller_b.vcf' \
     http://localhost:8000/api/compare

curl -F 'fasta=@samples/demo.fa' \
     -F 'vcf=@samples/evidence_demo.vcf' \
     -F 'sam=@samples/evidence_demo.sam' \
     http://localhost:8000/api/evidence
```

重新生成全部演示数据：`python3 -I scripts/make_demo.py`。

## 测试

```bash
# 核心逻辑（仅标准库）
python3 -I -m unittest discover -s tests -t . -p 'test_bio.py' -v

# 全部测试（Web 层需要先安装 requirements.txt）
python3 -I -m unittest discover -s tests -t . -v
```
