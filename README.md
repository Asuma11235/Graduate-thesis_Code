# Graduate-thesis_Code
卒業論文で使用した解析パイプライン一式．CoIP（IP-MS）と AirID のタンパク質相互作用データに対する GO 解析・可視化・文献件数集計・セット比較，snAtlas（単細胞）を入力にした WGCNA や発現可視化を含む．

## ディレクトリ概要
- `3章_IP-MSによるPPI解析/`
  - `GO analysis.R` : IP-MS（CoIP）ヒットの enrichGO + 可視化一式（bar・dot・Chord・階層ツリー・DAG）．
  - `Euler Diagram.R` / `Euler list.R` : 条件間集合のオイラー図作成と交差リスト出力．
  - `related paper search.R` : 候補遺伝子の PubMed 件数（全体 / brassinosteroid 文脈）を取得しスコアリング．
- `4章 _AirIDによるPPI解析/`
  - `GO analysis.R` : AirID 結果の enrichGO + 可視化（上記と同等仕様）．
  - `PCA & Hierarchical Clustering Analysis.R` : Abundance 行列の PCA・階層クラスタリング・相関ヒートマップ．
  - `PCA Analysis(4 Conditions).R` / `PCA Analysis.R` : prot_acc の有無マトリクスを用いた PCA と dendrogram．
  - `Hierarchical Clustering Heatmap.R` : Jaccard 距離のヒートマップを描画．
  - `Pseudo-Volcano Plot*.R` : 条件別の擬似ボルケーノプロット生成（Log2FC × Scaled Abundance）．
  - `Euler Diagram.R` / `Euler list.R` : AirID/CoIP 条件セットのオイラー図と交差リスト．
  - `Fisher's Exact Test for Set Overlap.R` : 条件間重なりの有意性検定．
- `5章_snAtlas再解析/`（Seurat RDS を入力とした単細胞再解析）
  - `batch_WGCNA.R` : セルタイプ/ステージで擬似バルク化 → WGCNA → BIL7 との相関モジュールとハブ遺伝子抽出．
  - `batch_WGCNA_GO analysis.R` : 相関モジュールに enrichGO を実行し CSV とドットプロットを出力．
  - `WGCNA on Pseudobulk Bins (Single|Merged Datasets).R` : データセット別／統合の WGCNA 補助スクリプト．
  - 発現可視化: `Celltype_expression barplot.R` / `Celltype_expression dotmap.R` / `Non-zero Expression Violin Plotter.R` / `Gene Expression Detection.R` などで AGI 指定の棒・ドット・バイオリン図．
  - 時系列: `Time-Course Gene Expression Plot_leaf.R` / `_flower-silique.R`．
  - 発現分布: `Organ_expression barplot.R`．
  - `Binomial & Fisher's Exact Test.R` : snAtlas 擬似バルクでの発現頻度・割合に対する検定．
  - `Monocle3 Pseudotime Analysis.R` : Monocle3 を用いた疑似時間解析と主要遺伝子のトラジェクトリ可視化．

## 共通の使い方
- すべて `Rscript` で実行．各スクリプト冒頭の設定（例: `IN_XLSX`, `RDS_PATH`, `OUT_*`, 閾値パラメータ）を手元のファイル構成に合わせて修正する．
- 必要パッケージはスクリプト側で自動インストールを試行（`clusterProfiler`, `org.At.tair.db`, `Seurat`, `WGCNA` など Bioconductor 依存を含む）．権限やプロキシ環境では事前インストール推奨．
- 典型的な入力列名例
  - GO 系: `SYMBOL` or TAIR ID と `Log2FC`．
  - セット比較系: `prot_acc`（擬似ボルケーノは `Log2FC_*` と `Abundances (Scaled): *`）．
- 出力例
  - enrichGO の CSV / bar・dot・Chord・階層ツリー・DAG PNG
  - 文献集計 Excel（スコア付）
  - PCA / dendrogram / 相関ヒートマップ / オイラー図 / 交差リスト Excel
  - WGCNA のモジュール相関・Hub 遺伝子 CSV、GO 結果 CSV/PDF、擬似バルク発現の図（棒・ドット・バイオリン・時系列）
  - 検定系: Binomial/Fisher の検定結果 CSV、Monocle3 の擬似時間図

## 第５章で使用したデータセット
- `Guo et al.,2025` - `"An Arabidopsis single-nucleus atlas decodes leaf senescence and nutrient allocation"`
- DOI: `10.1016/j.cell.2025.03.024`
- `Cell 188, 2856–2871.e1–e6, May 29, 2025`
- Publication History:
  - Received `October 15, 2023`; 
  - Revised `July 30, 2024`; 
  - Accepted `March 12, 2025`; 
  - Published online `April 11, 2025`
- 本論文が提供している，Seuratによるアノテーション済みRDSファイルを用いて解析を行った．
  - データ元： `https://ftp.cngb.org/pub/CNSA/ data6/CNP0002614/Other/`


## コード一覧（簡易解説）
- IP-MS 系: `GO analysis.R`, `Euler Diagram.R`, `Euler list.R`, `related paper search.R`
- AirID 系: `GO analysis.R`, `PCA & Hierarchical Clustering Analysis.R`, `PCA Analysis(4 Conditions).R`, `PCA Analysis.R`, `Hierarchical Clustering Heatmap.R`, `Pseudo-Volcano Plot (Bikinin vs. Mock|DMSO vs. Mock|Bikinin vs. DMSO).R`, `Volcano Plot.R`, `Euler Diagram.R`, `Euler list.R`, `Fisher's Exact Test for Set Overlap.R`
- snAtlas 再解析: `batch_WGCNA.R`, `batch_WGCNA_GO analysis.R`, `WGCNA on Pseudobulk Bins (Single Datasets|Merged Datasets).R`, `Celltype_expression barplot.R`, `Celltype_expression dotmap.R`, `Organ_expression barplot.R`, `Non-zero Expression Violin Plotter.R`, `Gene Expression Detection.R`, `Time-Course Gene Expression Plot_leaf.R`, `Time-Course Gene Expression Plot_flower-silique.R`, `Binomial & Fisher's Exact Test.R`, `Monocle3 Pseudotime Analysis.R`
