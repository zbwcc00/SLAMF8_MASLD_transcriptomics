# 投稿前第 1–3 步终审与修复报告

日期：2026-09-29  
目标期刊：Scientific Reports  
终审范围：主文、补充材料、主图与补充图、Supplementary Tables S1–S14、Source Data、作者信息、投稿声明和 Cover Letter。

## 总体结论

第 1–3 步已完成，当前投稿包达到可上传状态。终审未发现页面裁切、空白异常页、图中文字重叠、缺失面板、关键数字冲突、作者顺序错误或声明缺项。修复过程未改变任何科学结果、统计量或图形内容。

## 第 1 步：逐页和图件终审

- 主文 PDF：15 页；所有页面文本块均位于页面范围内，无裁切或异常空白页。
- Supplementary Information PDF：14 页；Figure S1–S12 均已嵌入并完整显示。
- 独立图件：Figure 1–5 和 Figure S1–S12，共 17 幅。
- 17 幅 PDF 图件全部通过文字与碰撞审计：无低于 5 pt 的字体、无 text-text collision、无 text-stroke collision、无 clipping。
- 最低字体为 Figure S1 的 5.0 pt；其余图件均高于该阈值。
- 所有 TIFF 均转换为 RGB、LZW 无损压缩、600 × 600 dpi；逐图 RGB 像素 SHA-256 在压缩前后完全一致。
- TIFF 总体积由 1418.8 MB 降至 30.0 MB，最大单图 4.7 MB，显著降低投稿系统上传失败风险。

## 第 2 步：正文、图注、补充表和 Source Data 一致性

- 关键数字一致性审计：24/24 项通过，0 项失败。
- Figure legends 顺序完整：Figure 1–5、Figure S1–S12。
- Supplementary Tables S1–S14 均存在，工作簿可正常读取；Table S6 的说明已与正文同步。
- FARG95 共 95 个唯一基因，包含 SLAMF8，不包含 SPP1、GPNMB、TREM2 和 CD9。
- 已明确说明：预设的五基因 anchor-exclusion 规则中，实际只有 SLAMF8 属于 FARG95，因此 anchor-excluded 与 SLAMF8-excluded 分数数值完全相同。
- 已修正文献溯源：ferro-aging 概念来源为 Liu et al., Cell Metabolism 2026（PMID 41819088；DOI 10.1016/j.cmet.2026.02.010）；FARG95 清单转录自 Lin et al., Biology Direct 2026 的 Supplementary Table 1（DOI 10.1186/s13062-026-00956-4）。
- 已删除误用的勘误 DOI 10.1016/j.cmet.2026.08.011，并将参考文献扩展为连续的 1–41；全部参考文献均在正文中被引用，无缺号或未引用条目。

## 第 3 步：作者信息和声明终审

- 作者顺序：Bowen Zheng、Guanghua Xie、Wangde Jin、Hao Li（通讯作者）。
- 全体作者单位统一为：Division of Hepatobiliary Pancreatic Surgery, The Affiliated Hospital of Yanbian University, Yanji 133000, China。
- 四位作者邮箱与用户提供信息一致。
- 四个 ORCID 均通过 ISO 7064 MOD 11-2 校验。
- CRediT 贡献符合既定排序：Bowen Zheng 承担主要研究与写作工作；Hao Li 负责总体构思、监督与项目管理；Guanghua Xie 和 Wangde Jin 按作者顺序承担方法、数据、验证和修订工作。
- Competing interests：None；Funding：None；Acknowledgements：None。
- Cover Letter 已删除未经作者确认的 “MD, PhD” 学位后缀，保留 Hao Li 作为 corresponding author 的准确身份。

## 最终格式指标

- 标题：18 words。
- 摘要：175 words，非结构式。
- Introduction + Results + Discussion：2942 words。
- 五个主图图注：114、146、131、133、110 words，均低于 350-word 上限。
- 主图与补充图均提供 TIFF 和 PDF；Supplementary Information、Supplementary Tables、Source Data、代码与作者信息均已纳入投稿包。

## 最终判定

PASS。当前版本可作为 Scientific Reports 正式投稿包使用。上传时应以最新 `SUBMIT_v4` 压缩包为准，不再使用此前的 v1–v3 包。
