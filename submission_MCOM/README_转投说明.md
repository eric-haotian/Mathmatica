# Mathematics of Computation 投稿包（2026-10-01 修订版）

基线：2026-09-29 的转投包（`author_checks/previous_versions/` 保存了上一版三个源文件）。期刊、标题、作者信息均不变。

## 本轮针对 desk reject 风险做了什么

编辑初筛主要看前两页、是否看得懂问题与贡献、是否符合期刊范围（"significant computational interest and substantial analysis or computational methodology"）。本轮改动全部针对这几点：

1. **摘要重写（274 词，低于 AMS 300 词上限）。** 第一句给出物理观测律 $Z(\omega)=R+i\omega L+\int(1+i\omega\tau)^{-1}d\mu(\tau)$，第二句说明计算问题，然后用平实语言列出四项结果。每一句都对应正文定理（命题 4.6、定理 4.2/4.3/5.2/6.2、推论 5.3、命题 6.3），没有新增主张。
2. **引言重构。** 新增 1.1 "The problem"：阻抗谱与弛豫时间分布（DRT）背景、未知校准参数 $R,L$、归一化模型、信息宽度与数值超额的定义，让编辑在半页内明白"算什么、为什么难、本文保证什么"。1.2 保留原主要结果概述；1.3 新增经典 Markov–Krein 矩理论（Karlin–Studden、Krein–Nudel'man）与正则化 DRT 反演文献作为对照；1.4 新增记号段，补定义了原稿未定义就使用的 $\mathcal P_n$、系数范数、各范数、"raw feasible" 与常数约定。
3. **各节导语。** 第 2–6 节开头各加一段说明本节做什么、为什么；定理 4.2 用到的 "baseline prediction distance" 在定理前明确定义；4.3 小节标题与两处行文中的 "worst-class" 改为 "worst-case"。
4. **数值部分更像计算数学论文。** 新增 Table 1（12 个等努力设计效率的认证区间）、Table 2（认证的期望努力下界）、Table 4（二次守卫网格节点数对比）和 Figure 1（精确响应差 $|\Delta Z(\omega)|$ 及认证极大点与一致上界）。所有数值均取自补充材料中已认证的证书；图由闭式公式绘制，脚本先与归档区间交叉核对再出图。**没有任何新实验。** 原停止路径表逐字节不变，现为 Table 3。
5. **结论增加局限性段落**（已知噪声水平与区间、类调常数、合成数据），并说明命题 4.4 的推广范围。
6. **AI 使用声明**措辞调整但范围不变，并明确作者核验全部结果并负全责（符合 AMS 披露要求）。
7. **参考文献新增 5 条**（均已核对出版信息），Wasilkowski–Woźniakowski 条目格式与其他期刊条目统一；原有 25 条全部保留。

## 没有改动的内容（已由脚本逐字节核验）

- 4 个定理、3 个引理、5 个命题、2 个推论、2 个注记、14 段证明；
- 原有 60 个显示公式块（顺序与内容）；算法 1；停止路径表；
- 补充材料 `code/` 下 210 个代码、输入、证书与结果文件（仅 `SECTION_MAP.txt`、`TRANSFER_NOTE.txt` 更新表号说明，并重新生成清单）。

记录见 `author_checks/source_integrity.json`、`author_checks/MCOM_2026-09-29_to_2026-10-01.diff`、`author_checks/supplement_transfer.json`。

## 交付文件

`upload/`：

1. `01_Main_Manuscript.pdf`：25 页，AMS `amsart` 10pt 默认版式（上一版 21 页；新增引言 1.1/1.4、三张表、一张图与结论段落）。未找到该刊页数上限。
2. `02_Cover_Letter.pdf`：1 页，日期 2026-10-01，面向 Mathematics of Computation，建议 Markus Bachmayr 处理（其副主编任期至 2030 年 1 月，已核对 AMS 编辑委员会页面）。
3. `03_Reproducibility_Materials.zip`：235 个文件（清单覆盖其余 234 个）。
4. `04_Reproducibility_Index.pdf`：1 页，补充材料导航，已按新表号更新。

`submission_fields.txt`：标题、作者、MSC、关键词、建议编辑与纯文本摘要，供投稿表单粘贴。

## 重建与复验

需要 TeX Live（`amsart`、`algorithm`、`algpseudocode`、`microtype`、`cleveref`、`hyperref`、`cm-super` 字体）与 Python 3.10+；完整核验另需 `pymupdf`；重绘图需 `matplotlib`、`numpy`。

```bash
cd submission_MCOM
python3 tools/assemble_package.py      # 编译三个 PDF、隔离重建比对、预检、完整性核验、重建补充包并复验、生成两级清单、verify_delivery.py
python3 verify_delivery.py             # 仅核对交付文件哈希
cd source && bash build.sh             # 仅编译
cd source/figures && python3 make_figure_response_difference.py   # 重绘 Figure 1（先与归档区间交叉核对）
```

补充材料解压后在 `Supplementary_Materials/` 内运行 `python3 run_all.py`（标准库即可）。本轮已在未解压目录与新 ZIP 的全新解压目录各复验一次：80 项证书任务、154 项检查（含 53 项拒绝检查）全部通过，无新实验、无新优化运行。

## 投稿前仍需作者本人确认

- 本稿未同时在其他期刊审理；FoCM 在审稿件的独立性如需披露由作者据实处理。
- 该刊任意 12 个月内每位作者最多 3 篇的限制。
- AI 使用声明与实际情况一致（本轮未改变声明范围）。
- 投稿系统的补充材料上传方式（PDF 初投为首选；ZIP 依系统实际栏目上传）。
- 本轮环境无法直接访问 ams.org，期刊政策条目来自官方页面的检索摘要，投稿前请在 ams.org 复核。

本轮没有操作投稿网站、邮箱或云盘；旧版文件保留在 `author_checks/previous_versions/`。
