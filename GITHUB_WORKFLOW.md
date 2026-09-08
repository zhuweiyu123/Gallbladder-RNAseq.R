# 在远程 VS Code 中使用 GitHub

工作目录：`/home/zhuweiyu/codex-r`  
仓库：https://github.com/zhuweiyu123/Gallbladder-RNAseq.R  
分支：`main`

## 日常操作
1. 在远程 VS Code 中打开上述工作目录。
2. 工作区干净时，先在“源代码管理”的“…”菜单中选择“拉取”。
3. 编辑脚本，查看文件差异，点击文件旁的“+”暂存需要提交的文件。
4. 输入提交说明并提交，再选择“推送”。

若还有未提交的修改，先提交或暂存（stash），再拉取。
拉取采用 fast-forward only；如分支分叉会停止，需先检查差异再合并。

## 跟踪范围
只跟踪根目录 Markdown/TXT 说明、Git 配置文件及 scripts/ 内的 R、Python、Shell、PowerShell、Markdown、TXT、Rmd、Qmd 文本文件。
数据、results/、logs/、tmp/、backups/ 和其他目录默认忽略。
不要使用 git add -f 强制添加被忽略的文件。
.gitattributes 将文本统一保存为 LF，减少 Windows 与 Linux 之间的换行差异。

服务器本地提交钩子拒绝新增或修改的单个文件超过 5 MiB。
该钩子只安装于当前服务器，不会自动传播到其他克隆。

## 认证
服务器使用此仓库专用的可写部署密钥；私钥留在服务器 ~/.ssh/ 下。
可在 GitHub 仓库 Settings → Deploy keys 撤销名为 codex-r server 192.168.1.6 的密钥。
