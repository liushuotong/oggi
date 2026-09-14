# tree 聚类与无标准答案的 auto 评分

本版本给出内部聚类质量评分，不需要标准分类答案。auto 已停用原 B/R/A/Q
加权几何平均；所有方法在固定的共同基因集合、同一距离矩阵上评价。

## 本次 bHLH 的服务器命令

先把更新后的整个 `cluster_engine/` 目录同步到服务器 `~/oggi_v1/cluster_engine/`。
将 Windows 上的 `bHLH.nwk` 上传到服务器
`~/oggi_v1/ath_genome_and_annotation/bHLH_try/bHLH.nwk`。
原始家族蛋白和映射文件继续使用 result 目录中的文件。

```bash
conda activate oggi
python -c "import numpy, pandas, scipy, Bio"
cd ~/oggi_v1/ath_genome_and_annotation/bHLH_try/result

python ~/oggi_v1/oggi.py cluster \
  -i gene_family.fa \
  --gene-map gene_to_assembly.tsv \
  --gene-tree ../bHLH.nwk \
  --collinear-pairs ../subcoli.collinear_pairs.tsv \
  -M auto --tree-threshold 0.1 \
  --ranking-metric silhouette \
  -t 32 -o ../cluster_auto_tree_v1
```

若仅运行树聚类，将 `-M auto` 换成 `-M tree`，并使用新的输出目录。
`--tree` 仍表示物种/材料树；基因树必须使用 `--gene-tree`。
没有完整蛋白组或明确复用的完整蛋白组 OrthoFinder 结果时，自动跳过
OrthoFinder 及两种 OrthoFinder 组合方法。加权 MCL 需要共线性等加权证据。

新评分版本与旧结果的指纹不同，请使用新目录，不对旧运行添加 `--resume`。
同一代码、参数和输入的后续恢复可以使用 `--resume`。

## tree 的算法

1. 用 Biopython 解析一个 Newick 文件，严格检查叶子 ID 和分支长度。
2. 计算叶子间树路径长度之和，不要求树已定根；内部支持率标签不作为距离。
3. 使用 SciPy complete-linkage 层次聚类，在指定最大簇内距离处切割。
4. 树中缺失的输入基因保留为未解决单例，同时报告覆盖率。

该方法只使用基因树的拓扑与枝长。不是 TreeCluster 软件的实现，不保证每簇
一定为严格单系群，也不执行复制/物种分化判别。0.1 的单位是输入树的枝长单位，
通常为每位点替换数；它是初始参数，不是通用 OGG 生物学阈值。

默认各方法一个候选参数。若要检查 tree 参数敏感性，可另存配置：

```json
{
  "ranking_metric": "silhouette",
  "grid": {
    "tree": [
      {"threshold": 0.025}, {"threshold": 0.05}, {"threshold": 0.1},
      {"threshold": 0.2}, {"threshold": 0.4}
    ]
  }
}
```

在命令中加入 `--config 配置文件.json`。论文需要报告各方法参数搜索范围和数量，
不可把某方法充分调参后的最高值直接当作相同调参条件下的公平比较。

## 分数和输出

| 输出指标 | 含义 | 如何使用 |
| --- | --- | --- |
| silhouette_mean | 平均轮廓系数，范围 [-1,1] | 默认主指标，越大越好 |
| total_score | `50 × (silhouette_mean + 1)` | 默认 auto 的 0–100 分；50 对应轮廓系数 0 |
| dunn_index | 最小簇间距离 / 最大簇内直径 | 辅助指标，越大越好，对极端基因对敏感 |
| dunn_bounded | `Dunn / (1 + Dunn)` | Dunn 的有界形式，处理无限大情况 |
| evaluation_coverage | 共同可评价基因数 / 输入基因数 | 所有方法使用同一集合 |
| singleton_gene_fraction | 评价集合中处于单例簇的基因比例 | 辅助判断过度拆分 |
| assignment_coverage | 未被方法适配器标记为未解决的基因比例 | 运行诊断，不代表准确率 |

通过 `--ranking-metric dunn` 可以明确改用 `100 × dunn_bounded` 排序；
Silhouette 和 Dunn 始终同时报告，不将两者随意加权。

单例基因的 Silhouette 为 0；全部合并或全部单例的候选不参与排名。
某个候选不可评分，不影响其他候选。Dunn 分母为 0 且分子为正时，
原始值显示 NA 并标记为无界，有界形式为 1；0/0 则保持不可评价。

默认只报告覆盖率，不强行设置通用最低覆盖标准。可在配置中预先设置
`min_evaluation_coverage`；低于该标准仍报告数值，但不给推荐。

- `scores.tsv`：每个方法/参数的分数、两个指数、覆盖率及簇数。
- `selection.json`：推荐结果、主指标、并列结果和评价范围。
- `selected_clusters.tsv`：最佳分组；不同完整划分并列或无可用分数时只含表头。
- `method_summary.tsv`：各方法参数候选数、最高分、分数范围及簇数范围。
- `evaluation_genes.tsv`：每个输入基因是否进入固定评价集合。
- `candidates/*/clusters.tsv`：每个候选的完整分组。
- `candidates/*/cluster_statistics.tsv`：每簇平均树距离、直径、平均轮廓系数。
- `candidates/*/silhouette.tsv`：逐基因轮廓系数。

相同完整划分得到相同分数，并作为等价结果输出；不因方法名称不同给出不同分数。
在评价子集上相同、但完整输入上的划分不同，仍视为不同结果。

这些分数适用于比较同一数据、同一距离定义下的内部聚类质量。若树既用于
tree 聚类又用于评分，这是对该树距离的内部拟合评价，不是独立的生物学验证，
0–100 分也不是直系同源准确率。论文建议同时报告原始指数、覆盖率、簇数和参数；
仅凭最高分不能宣称差异具有统计显著性。公式、参考文献及其他参数见
[CLUSTER_AUTO.md](CLUSTER_AUTO.md)。
