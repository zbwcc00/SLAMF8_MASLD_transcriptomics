# Supplementary Tables S1–S14

本目录由 `scripts/053_build_supplementary_tables.py` 从冻结的结果文件汇总而成，不新增生物学分析。

- 单细胞 QC 是 library/cell 单位；样本级比较与相关分析以 patient/sample 为单位。
- GSE192741 单独列为 Table S14，与 Figure S12 对应；spot 是混合测量单位，不作为患者级独立重复。
- LIANA、CellChat 和 NicheNet 是转录组推断或候选排序，不能单独证明空间邻近、蛋白分泌或因果。
- Table S12 是反事实网络预测，不是实验性 SLAMF8 敲除。
- Table S11 的 SAMac-like signature specificity audit（leave-one-gene-out 和表达量十分位匹配随机基因集）是敏感性分析，不是独立队列。
- 每张表均保留 `source_file`、`analysis_unit` 与 `interpretation_note` 以保证可追溯性。

S1 队列；S2 QC/双细胞；S3 组成；S4 注释；S5 拟时序；S6 SLAMF8-high/low；S7 巨噬细胞–HSC 关联；S8 通讯；S9 SPP1 审计；S10 NicheNet；S11 bulk；S12 虚拟扰动；S13 冻结与单位审计；S14 GSE192741 Visium 空间共表达审计。
