# Codex 项目记忆

根目录 `AGENTS.md` 要求启动时读取本文件、`PROJECT_STATE.md`、`NEXT_ACTIONS.md` 和 `WORKLOG.md`。因此前三份只保留当前入口所需信息；逐句审稿细节已从活动记忆中清理，当前复核要点见 `MATCH_NARRATIVE_LOGIC_REVIEW_2026-10-18.md`。验证原始计数时查 `validation/`。

2026-11-10 全文术语校准入口：`MATCH_TERMINOLOGY_CALIBRATION_2026-11-10.md`，追查见 `validation/MATCH_TERMINOLOGY_FOLLOWUP_2026-11-10.md`。已按作者决定修复 PLD/LCD、将十节点模板表移至附录，并保留 topology recovery 原文；其余建议不自动视为作者已接受。

剩余术语的第二轮复核见 `validation/MATCH_TERMINOLOGY_ROUND2_2026-11-10.md`；R2-02、R2-03 已于 2026-11-11 按作者批准落实，其余建议尚未应用。摘要九句与后文的术语/词频审查见 `validation/MATCH_ABSTRACT_TERM_AUDIT_2026-11-11.md`。

枚举算法的统一定义、状态、伪代码和正确性入口现为 `docs/subgroup_dfs_trie_algorithm.md`。2026-11-16 C++ 完整代表去重已切为压缩数组 Trie（ABI 3），实测见 `validation/COMPACT_TRIE_REVIEW_2026-11-16.md`；历史分阶段性能记录继续保留，不混用不同版本数字。

2026-11-23 精确ρ配额已接入优化Python/C++核（ABI 4），实现和验证见 `native/README.md`、`validation/RHO_NATIVE_REVIEW_2026-11-23.json`。论文已落实作者明确批准的1/2/3/5/9/11/13/14项（1/14为后续批准）；最新讨论状态位于原拟稿顶部。取消“普通调度”后的八处方法联动也已于2026-11-23按批准落稿，后续主要讨论CASM性能正文和表格。最新31项审读审批及15项对比见 `validation/MATCH_REVIEW_DISCUSSION_2026-11-23.md`。作者随后批准当前ABI-4/新标签池重跑，主表/附表/vba旧时间已全部替换；默认后端改为auto。沟通须通俗、短句、少字。

2026-11-24：上述31项审读已处理完毕，新批准十四项已落实；6/13/14/15/22/24以及此前7/9保留。最新作用范围入口为 `validation/MATCH_APPROVED_REVIEW_SCOPE_2026-11-24.json`。主稿33页，旧讨论稿的待定标签不再代表当前状态。

## 文档分工

- `PROJECT_STATE.md`：只记当前有效的方法、论文口径、结果边界和运行约束。
- `NEXT_ACTIONS.md`：只保留尚未完成且仍有价值的任务；完成后删除，不留“complete”清单。
- `WORKLOG.md`：阶段级历史摘要。具体逐轮差异以 Git 历史和论文源码为准。
- `validation/`：可复核的证据、审计脚本和原始材料。

更新记忆时优先修改现状和待办；工作日志只记录值得后续会话继承的阶段结果，避免把同一变化复制到多个文件。不要把生成物、图片或临时脚本当作启动材料。
