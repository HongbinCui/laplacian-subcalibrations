# GitHub 上传说明

此包用于仓库 HongbinCui/laplacian-subcalibrations。

1. 解压 ZIP，把 `laplacian-subcalibrations` 文件夹**里面的内容**放在仓库根目录。
   根目录应直接显示 README.md、anc、paper，而不是再套一层文件夹。
2. 保留完整 anc 目录，尤其是 manuscript 中的冻结参考 PDF、labels.aux
   以及 results 中的本次运行记录。它们共同构成已验证的发布材料。
3. 更新已有仓库时，保留原有 LICENSE、仓库配置和 Git 历史。
   本包没有检查或授权清理远端旧文件，不建议把仓库整体删除重建。
   如果已有旧版入口，请让首页明确指向本包的当前入口。
4. 完成提交后，建议创建固定 Release，标签使用 verification-20261012-v15，
   并把本 ZIP 作为附件。RELEASE_NOTES.md 可作为发布说明的基础。
5. 发布后确认 README 中的 PDF、运行摘要和 notebook 链接可以访问。
   然后在 arXiv 第二版中使用该固定 Release 的链接。

当前论文是 V15；计算套件标识 20261011_memo_v7 是有意保留的，不能简单批量改名。
包中旧准备日期不等于运行日期；真实运行日期是 2026-10-12，根目录 README 已说明。
本次只整理上传包，并未重新运行 Mathematica，也未自动上传到 GitHub。

数学验证全部使用 Mathematica。SHA256SUMS 仅用于文件完整性检查。
