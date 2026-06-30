# Codex Project Memory

这个目录用于保存 MOF 结构枚举与对称性分析项目的可持续工作记忆，方便之后任何 Codex 会话快速接上当前进度。

建议每次开始新会话时先让 Codex 读取：

```text
codex/README.md
codex/PROJECT_STATE.md
codex/WORKLOG.md
```

## 文件说明

- `PROJECT_STATE.md`: 当前项目事实、代码状态、DFT 设计、服务器运行方式、论文口径。它是接续工作的主入口。
- `WORKLOG.md`: 按日期追加工作记录。每次完成关键修改、得到 DFT 结果、改变论文表述或修复代码后，都追加一条。
- `NEXT_ACTIONS.md`: 下一步任务清单和注意事项。

## 更新规则

1. 做完任何会影响项目方向、实验设计、服务器运行、论文叙述或代码行为的改动后，更新 `PROJECT_STATE.md` 中对应段落。
2. 每次阶段性工作结束时，在 `WORKLOG.md` 追加日期、改动、验证结果、下一步。
3. 如果下一步优先级变化，更新 `NEXT_ACTIONS.md`。
4. 保持论文与代码一致，尤其是 DFT 方法：PBE+D3、500/60 Ry、single-point、composition-controlled stratified sampling。

