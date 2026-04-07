
python change.py                                
dir=/home/lipeixin27/data/index_genome/cisTarget_databases/  # 改成你自己下载数据库的位置
tfs=$dir/hs_hgnc_tfs.txt   # 包含所有人类 TF 的基因名（~1839 个）
feather=$dir/hg19-tss-centered-10kb-10species.mc9nr.genes_vs_motifs.rankings.feather  # motif 排名文件
tbl=$dir/motifs-v9-nr.hgnc-m0.001-o0.0.tbl  # motif 注释信息表
input_loom=./sample.loom  # 你的输入表达矩阵，需为 loom 格式
ls $tfs $feather $tbl
pyscenic grn \
--num_workers 10 \
--output adj.sample.tsv \
--method grnboost2 \
sample.loom \
$tfs
pyscenic ctx \
adj.sample.tsv $feather \
--annotations_fname $tbl \
--expression_mtx_fname $input_loom \
--mode "dask_multiprocessing" \
--output reg.csv \
--num_workers 20 \
--mask_dropouts
pyscenic aucell \
$input_loom \
reg.csv \
--output out_SCENIC.loom \
--num_workers 10
