"""
内嵌修复 Agent —— 通过 OpenAI API 诊断并修复 MOF 组装代码问题。

调用方式：
    agent = RepairAgent()
    success = agent.run(problem, retest_fn=my_retest_func)

每个问题独立一轮 agent 循环（Step 1→5），支持连续修复多个问题。
"""

import os
import sys
import re
import openai

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


def _load_env_file():
    """Load .env file from project root into os.environ if present."""
    env_path = os.path.join(PROJECT_ROOT, ".env")
    if not os.path.exists(env_path):
        return
    with open(env_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, _, value = line.partition("=")
            key = key.strip()
            value = value.strip().strip('"').strip("'")
            if key and key not in os.environ:
                os.environ[key] = value


_load_env_file()

SYSTEM_PROMPT = """\
You are a MOF crystal structure assembly debugging expert.
Your task is to diagnose and fix bugs in a Python codebase that assembles MOF (Metal-Organic Framework) structures.

Core pipeline:
  topology graph → symmetry analysis → enumeration → verification → assembly → CIF output

Key modules:
  - assembler.py: StructureAssembler — places nodes/edges in 3D, builds structure
  - main.py: run_enumeration, enrich_graph_data_with_topology — pipeline orchestration
  - enumerator.py: MOFEnumerator — enumerates valid colorings
  - symmetry.py: SymmetryEngine — automorphism group computation
  - structure_validator.py: StructureValidator — post-assembly validation

Rules:
- Only analyze code that has been explicitly provided to you
- Do not guess or invent implementation details
- When proposing a fix, include exact file path, old snippet, and new snippet
- Respond in Chinese (the user speaks Chinese)
- Keep responses concise and actionable"""


class RepairAgent:
    def __init__(self):
        api_key = os.environ.get("OPENAI_API_KEY")
        if not api_key:
            raise RuntimeError("OPENAI_API_KEY 环境变量未设置，无法初始化 RepairAgent")
        self.client = openai.OpenAI(api_key=api_key)
        self.model = "gpt-4o"
        self.messages = []          # 多轮对话历史
        self.extra_reads = 0        # "继续找" 计数
        self.max_extra_reads = 3

    # ── API 调用 ──

    def _call_api(self, user_message):
        """调用 API，追加到对话历史。返回响应文本或 None。"""
        self.messages.append({"role": "user", "content": user_message})
        try:
            response = self.client.chat.completions.create(
                model=self.model,
                max_tokens=2048,
                messages=[{"role": "system", "content": SYSTEM_PROMPT}] + self.messages,
            )
            text = response.choices[0].message.content
            self.messages.append({"role": "assistant", "content": text})
            return text
        except Exception as e:
            print(f"  [API 错误] {e}")
            return None

    # ── 代码读取 ──

    def _read_function(self, module_path, func_name, max_lines=60):
        """
        读取指定函数/方法的实现代码（最多 max_lines 行）。

        从 'def func_name' 或 'class ClassName' 开始，
        读到下一个同级定义或文件末尾。
        """
        if not os.path.exists(module_path):
            return None, f"文件不存在: {module_path}"

        with open(module_path, "r") as f:
            all_lines = f.readlines()

        # 查找定义行
        start_line = None
        pattern = re.compile(rf"^\s*(def|class)\s+{re.escape(func_name)}\b")
        for i, line in enumerate(all_lines):
            if pattern.match(line):
                start_line = i
                break

        if start_line is None:
            # 尝试部分匹配（如 "assemble_and_save" 匹配 "def assemble_and_save"）
            for i, line in enumerate(all_lines):
                stripped = line.rstrip()
                if (stripped.startswith("def ") or stripped.startswith("class ")) and func_name in stripped:
                    # 对方法名做更精确匹配
                    if re.search(rf"\b{re.escape(func_name)}\b", stripped):
                        start_line = i
                        break

        if start_line is None:
            return None, f"未找到定义: {func_name}"

        # 确定基准缩进
        base_indent = len(all_lines[start_line]) - len(all_lines[start_line].lstrip())

        # 读取到下一个同级定义或空行后的同级定义
        end_line = len(all_lines)
        for i in range(start_line + 1, len(all_lines)):
            line = all_lines[i]
            stripped = line.rstrip()
            if stripped == "" or stripped.startswith("#"):
                continue
            indent = len(line) - len(line.lstrip())
            # 遇到同级或更外层的 def/class
            if indent <= base_indent and (stripped.startswith("def ") or stripped.startswith("class ")):
                end_line = i
                break

        # 限制行数
        if end_line - start_line > max_lines:
            end_line = start_line + max_lines

        code = "".join(all_lines[start_line:end_line])
        return code, None

    def _get_function_signatures(self, module_path):
        """读取模块文件，只提取函数/类/方法签名行。"""
        sigs = []
        if not os.path.exists(module_path):
            return sigs
        with open(module_path, "r") as f:
            for line in f:
                stripped = line.rstrip()
                if stripped.startswith("class ") or stripped.startswith("    def ") or stripped.startswith("def "):
                    sigs.append(stripped)
        return sigs

    # ── 补丁应用 ──

    def _apply_patch(self, file_path, old_string, new_string):
        """在文件中执行字符串替换。返回 True 表示成功。"""
        abs_path = os.path.join(PROJECT_ROOT, file_path)
        if not os.path.exists(abs_path):
            print(f"  [Patch 错误] 文件不存在: {abs_path}")
            return False

        with open(abs_path, "r") as f:
            content = f.read()

        count = content.count(old_string)
        if count == 0:
            print(f"  [Patch 错误] 未找到匹配的 old_string，可能文件已变更")
            return False
        if count > 1:
            print(f"  [Patch 错误] old_string 在文件中出现了 {count} 次，需更精确的上下文")
            return False

        new_content = content.replace(old_string, new_string)
        with open(abs_path, "w") as f:
            f.write(new_content)
        print(f"  [Patch] 已修改 {file_path}")
        return True

    # ── 解析 API 响应 ──

    def _parse_patch(self, response_text):
        """从 API 响应中解析补丁信息。返回 (file_path, old_string, new_string) 或 None。"""
        # 期望格式:
        # FILE: <path>
        # OLD: <<<OLD
        # <code>
        # OLD
        # NEW: <<<NEW
        # <code>
        # NEW
        file_match = re.search(r'FILE:\s*(\S+)', response_text)
        old_match = re.search(r'<<<OLD\s*\n(.*?)\nOLD', response_text, re.DOTALL)
        new_match = re.search(r'<<<NEW\s*\n(.*?)\nNEW', response_text, re.DOTALL)

        if file_match and old_match and new_match:
            return file_match.group(1), old_match.group(1), new_match.group(1)
        return None

    # ── Agent 主循环 ──

    def run(self, problem, retest_fn=None):
        """
        运行单个问题的完整修复循环。

        problem: dict with keys:
            case_name, ref_file, category, label, status, detail,
            module (optional), functions (optional, list), note (optional)

        retest_fn: callable(problem) -> new_report_or_None

        返回 True 表示修复成功。
        """
        self.messages = []
        self.extra_reads = 0

        case_name = problem.get("case_name", "?")
        label = problem.get("label", "?")
        detail = problem.get("detail", "")
        category = problem.get("category", "comparison")

        print(f"\n{'=' * 60}")
        print(f"修复 Agent: {case_name}")
        print(f"问题: {label}")
        print(f"{'=' * 60}")

        # ── Step 1: 问题确认 ──
        note = problem.get("note", "")
        mod_name = problem.get("module", "")
        funcs = problem.get("functions", [])

        if note:
            print(f"\n  Agent 认为根因是: {note}")
        if mod_name:
            print(f"  关联模块: {mod_name}")
        if funcs:
            print(f"  关联函数: {', '.join(funcs)}")

        while True:
            choice = input("\n  [Step 1/5] 根因假设是否正确？[y/继续找/跳过]: ").strip().lower()
            if choice == "y":
                break
            elif choice == "跳过":
                print("  已跳过此问题。")
                return False
            elif choice == "继续找" or choice == "c":
                if self.extra_reads >= self.max_extra_reads:
                    print(f"  已达分析上限 ({self.max_extra_reads} 次追加)，建议人工介入。")
                    return False
                self.extra_reads += 1
                # 让 API 建议下一个要读的函数
                self._step2_read_more(mod_name, funcs, detail)
            else:
                print("  请输入 y / 继续找 / 跳过")

        # ── Step 2: Agent 读代码定位 ──
        fix_plan = self._step2_locate_root_cause(problem)
        if not fix_plan:
            print("  [Step 2] 未能生成修复方案。")
            return False

        # ── Step 3: 修改确认 ──
        print(f"\n  [Step 3/5] 修复方案:")
        print(f"  {'─' * 50}")
        print(f"  {fix_plan}")
        print(f"  {'─' * 50}")

        # 判断是否需要强制确认
        needs_force_confirm = any(
            kw in fix_plan for kw in ["坐标变换", "拓扑", "lattice", "fractional", "Bstar", "alpha_matrix"]
        )
        if needs_force_confirm:
            print("  【⚠️需确认】涉及坐标变换/拓扑核心算法，必须确认后才执行。")
            while True:
                ans = input("  确认修改？[y/n]: ").strip().lower()
                if ans == "y":
                    break
                elif ans == "n":
                    print("  已取消修改。")
                    return False
        else:
            ans = input("  确认修改？[y/n]: ").strip().lower()
            if ans != "y":
                print("  已取消修改。")
                return False

        # ── Step 4: 执行修改 ──
        print(f"\n  [Step 4/5] 执行修改...")
        success = self._step4_execute_fix(fix_plan, problem)
        if not success:
            print("  [Step 4] 修改未成功应用。")
            return False

        # ── Step 5: 重新测试 ──
        print(f"\n  [Step 5/5] 重新测试...")
        if retest_fn:
            new_report = retest_fn(problem)
            if new_report:
                print("  新的对比结果:")
                print(new_report.summary())
                if not new_report.has_errors:
                    print("  ✅ 修复成功！所有严重项已消除。")
                    return True
                else:
                    # 检查原始问题是否解决
                    original_label = problem.get("label", "")
                    still_fails = any(
                        label == original_label
                        for label, status, _ in new_report.checks
                        if status == new_report.STATUS_BAD
                    )
                    if not still_fails:
                        print("  ⚠️ 原始问题已修复，但仍有其他严重项。请人工检查。")
                        return True
                    else:
                        print("  ❌ 原始问题仍然存在，修复可能未生效。")
                        return False
            else:
                print("  ❌ 重测试失败（组装报错或其他异常）。")
                return False
        else:
            print("  (无重测试回调，跳过验证)")

        return True

    # ── Step 2 子步骤 ──

    def _step2_read_more(self, mod_name, funcs, detail):
        """'继续找' 时：让 API 建议下一个函数，读取并追加分析。"""
        sigs_text = ""
        mod_path = self._resolve_module_path(mod_name)
        if mod_path:
            sigs = self._get_function_signatures(mod_path)
            sigs_text = "\n".join(sigs)

        prompt = f"""\
当前问题: {detail}

关联模块 {mod_name} 的函数签名:
{sigs_text}

我们已经查看了函数: {', '.join(funcs)} (推测的关联函数)。

请建议下一个需要查看的具体函数名（只回复函数名，不要加解释）。"""

        resp = self._call_api(prompt)
        if not resp:
            return

        # 提取函数名
        next_func = resp.strip().split("\n")[0].strip()
        print(f"  Agent 建议查看: {next_func}")

        code, err = self._read_function(mod_path or "", next_func, max_lines=60)
        if err:
            print(f"  [读取失败] {err}")
            return
        print(f"  已读取 {next_func} ({len(code.splitlines())} 行)")

        funcs.append(next_func)

    def _step2_locate_root_cause(self, problem):
        """Step 2 核心：通过 API 定位根因并生成修复方案。"""
        label = problem.get("label", "")
        detail = problem.get("detail", "")
        category = problem.get("category", "comparison")
        mod_name = problem.get("module", "")
        funcs = problem.get("functions", [])
        note = problem.get("note", "")

        # 1. 构建初始上下文
        mod_path = self._resolve_module_path(mod_name)
        if not mod_path:
            print(f"  [Step 2] 无法定位模块: {mod_name}")
            return None

        sigs = self._get_function_signatures(mod_path)
        sigs_text = "\n".join(sigs) if sigs else "(无签名) "

        # 首次请求：给定问题+签名，让 API 选函数
        prompt = f"""\
## 问题描述
检查项: {label}
详情: {detail}
分类: {category}
根因假设: {note}

## 模块: {mod_name}
函数签名:
{sigs_text}

## 任务
根据问题描述和函数签名，判断需要查看哪个具体函数来确认根因。
只回复一个函数名（不要加解释）。"""

        resp = self._call_api(prompt)
        if not resp:
            return None

        target_func = resp.strip().split("\n")[0].strip()
        print(f"  Agent 选择查看: {target_func}")

        # 2. 读取该函数实现
        code, err = self._read_function(mod_path, target_func, max_lines=60)
        if err:
            print(f"  [读取失败] {err}")
            return None
        print(f"  已读取 {target_func} ({len(code.splitlines())} 行)")

        # 3. 传入函数实现，让 API 分析根因+修复方案
        prompt2 = f"""\
## 函数: {target_func} (来自 {mod_name})
```python
{code}
```

## 问题
检查项: {label}
详情: {detail}
根因假设: {note}

## 任务
1. 确认根因是否在此函数中
2. 如果是，给出完整修复方案（包括具体代码改动）
3. 如果不是，说明需要查看哪个其他函数

回复格式:
根因确认: <是/否，及解释>
修复方案: <具体步骤>
需要修改的文件: <文件路径>
修改内容: <old snippet → new snippet>"""

        resp2 = self._call_api(prompt2)
        if not resp2:
            return None

        return resp2

    def _step4_execute_fix(self, fix_plan, problem):
        """Step 4: 让 API 生成精确补丁并应用。"""
        mod_name = problem.get("module", "")

        prompt = f"""\
## 修复方案（已确认）
{fix_plan}

## 任务
生成精确的代码补丁。使用以下格式:

FILE: <文件相对路径>
<<<OLD
<需要替换的精确原代码>
OLD
<<<NEW
<替换后的新代码>
NEW

注意:
- OLD 代码必须与源文件完全一致（包括缩进和空行）
- 只输出一个补丁
- 确保 NEW 代码修复问题且不引入新错误"""

        resp = self._call_api(prompt)
        if not resp:
            return False

        patch = self._parse_patch(resp)
        if not patch:
            print("  [Step 4] 无法从 API 响应中解析补丁。响应内容:")
            print(f"  {resp[:500]}")
            return False

        file_path, old_str, new_str = patch
        print(f"  目标文件: {file_path}")
        return self._apply_patch(file_path, old_str, new_str)

    def _resolve_module_path(self, mod_name):
        """将模块名解析为绝对路径。"""
        if not mod_name:
            return None
        # 处理 "assembler.py / main.py" 这种多模块情况——取第一个
        parts = [p.strip() for p in mod_name.replace("/", ",").split(",")]
        for part in parts:
            if part.endswith(".py"):
                path = os.path.join(PROJECT_ROOT, part)
                if os.path.exists(path):
                    return path
        return None
