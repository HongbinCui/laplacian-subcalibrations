# GitHub 链接同步说明（2026-10-13）

论文完整题目保持为：
Laplacian Algebras and Subcalibrations for Minimal Hypercones

本次将 GitHub 完整包中的论文 PDF、LaTeX 源码和 README 的引用同步为
verification-20261013-v2，并更新发布说明和校验清单。
PDF 与本次修正后的 arXiv PDF 完全相同，仍为 166 页。
计算代码、运行记录和 anc 内全部 72 个文件保持不变，无需重跑。

## 1. 更新仓库 main 中的文件

1. 解压本次新下载的 Laplacian_GitHub_update_20261013.zip。
2. 在仓库根目录选择 Add file > Upload files。
3. 拖入解压后根目录的四个文件和整个 paper 文件夹。
   上传列表应显示 README.md、RELEASE_NOTES.md、UPLOAD_GUIDE_ZH.md、
   SHA256SUMS、paper 下的 PDF，以及
   paper/source/sections/computer_verification.tex。
4. 共 6 个文件；提交后保留已有 anc 和其它未改动文件。

## 2. 替换现有 Release 的完整附件

1. 打开 verification-20261013-v2 Release，点击编辑按钮。
2. 删除该 Release 中旧的同名 Laplacian_GitHub_complete_20261013.zip 附件。
3. 上传本次新下载的 Laplacian_GitHub_complete_20261013.zip。
4. 可将发布说明更新为 RELEASE_NOTES.md 的内容，点击 Update release。

不必创建新 Release，也不必改论文题目、Release 名称或标签。
同名附件替换后，原下载地址可以继续使用。
arXiv 使用上一轮已经提供的、指向 verification-20261013-v2 的源码包。

GitHub 自动生成的 Source code (zip/tar.gz) 取自标签所指的提交，
不会随手工附件或 main 的更新而自动更新。本次修正版请下载明确上传的
Laplacian_GitHub_complete_20261013.zip；本次不移动已发布标签。
paper 内的 PDF 文件名沿用旧路径以保持链接有效，文件名中的 v15
不代表 PDF 正文仍引用旧 Release。
