# 多机协作

## 新机器首次使用

```sh
git clone https://github.com/ttttxx/Symmetry-Aware-Enumeration-of-mof.git
cd Symmetry-Aware-Enumeration-of-mof
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Windows 激活命令为 `.venv\Scripts\activate`。Python 依赖在每台机器单独安装；密钥通过本机环境变量或本机 `.env` 配置。

macOS/Linux 已安装 C++17 编译器时可构建原生核：

```sh
python tools/build_enumeration_core.py
```

原生库按机器重新构建。默认后端为 `auto`，原生库不可用时回退优化 Python；Windows 可直接使用 Python 后端。

## 每次开始工作

先提交或暂存本机未完成的修改，再拉取更新：

```sh
git switch master
git pull --ff-only
```

Codex 接续工作前读取 `AGENTS.md` 及 `codex/` 四个主文档。

## 每次完成工作

```sh
git status
git add <本次修改的文件或目录>
git commit -m "描述本次修改"
git push origin master
```

GitHub 同步由 commit/push 和 pull 完成，不是实时文件同步。切换机器前先推送；到另一台机器先拉取。两台机器同时修改时，各自创建不同分支：

```sh
git switch -c codex/my-change
# 修改、提交后
git push -u origin codex/my-change
```

通过 Pull Request 合并。遇到 `pull --ff-only` 失败时检查两端提交并处理分叉，避免强制推送。

## 同步范围

仓库包含源码、输入模板、论文、项目记忆、验证材料和当前结果数据。本机 `.env`、虚拟环境、缓存、编译日志、临时渲染目录和原生二进制由 `.gitignore` 排除。克隆后重新安装依赖并构建原生库。

此开发仓库与 `data_availability/` 中的投稿数据包具有不同用途；开发同步不代表已完成投稿数据仓库公开核验。
