# Graduate-thesis_Code
卒業論文で実施した PPI 解析まわりのスクリプト置き場です。主に CoIP（IP-MS）と AirID のデータを対象に、GO 解析・可視化・文献数集計・セット間比較を行います。さらに snAtlas 再解析（Seurat RDS を入力とする WGCNA / 発現量可視化）を収録しています。

## 構成
- `3章_IP-MSによるPPI解析/`
  - `GO analysis.R` : CoIP 結果の GO enrichment と可視化（bar, GOplot, GO 階層ツリー、DAG）。
  - `related paper search.R` : 候補遺伝子に対する PubMed 件数（一般・BR 文脈）を取得し、スコアリングして Excel 出力。
- `4章 _AirIDによるPPI解析/`
  - `GO analysis.R` : AirID 結果の GO enrichment と可視化（上と同等仕様）。
  - `PCA分析.R` : AirID 4 条件 + CoIP 4 条件の `prot_acc` を 0/1 行列にまとめ、PCA・階層クラスタリング・Jaccard ヒートマップを描画。
  - `ベン図.R` : AirID/CoIP 8 条件の集合を Euler 図（近似ベン図）として描画。
  - `ベン図_リスト.R` : AirID×CoIP の交差（16 通り）を `prot_acc` と SYMBOL 付きで Excel 出力し、件数と Jaccard 行列も併記。
- `5章 _snAtlas再解析/`（Seurat RDS を入力とした単細胞データ再解析）
  - `batch_WGCNA.R` : 擬似バルク化（セルタイプ×ステージなど）→ WGCNA → BIL7 との相関モジュール抽出と Hub 遺伝子抽出。
  - `batch_WGCNA_GO analysis.R` : WGCNA で得たモジュールに対し、相関しきい値でフィルタして GO enrichGO を実行（ドットプロットと CSV 出力）。
  - `器官別発現量解析.R` : 複数 Seurat RDS を読み込み、全セル合算の擬似バルク CPM を計算して棒グラフ出力。
  - `細胞型別発現量解析-ドット図.R` : セルタイプ（＋ステージ任意）ごとの平均発現・検出率を計算し、指定 AGI のドットプロットを PDF に保存。
  - `細胞型別発現量解析-棒グラフ.R` : 同じメトリクスを用い、指定 AGI の平均発現を棒グラフで比較（log10 変換オプションあり）。

## 使い方のメモ
- いずれも Rscript 想定。スクリプト冒頭のパス（例: `IN_XLSX`, `OUT_*`）を手元のファイルに合わせて書き換えてください。
- 必要パッケージはスクリプト内で自動インストールするようにしています（Bioconductor: `clusterProfiler`, `org.At.tair.db` など）。
- 入力 Excel の列名想定:
  - GO 解析: `SYMBOL` と `Log2FC`（TAIR ID でも可）。
  - 交差・可視化系: `prot_acc`（`ベン図.R` はデフォルトで `SYMBOL` 列）。
- Seurat RDS を使うスクリプト（snAtlas 再解析）は、`RDS_PATH` や `CELLTYPE_COL`（セルタイプ列）、`STAGE_COL`（任意）を環境に合わせて指定してください。`input_data/` や `results_*` は任意の出力ディレクトリに変更可能です。
- 主な出力:
  - GO 解析: CSV（enrichGO 結果・マッピング監査）、PNG（bar/dot/Chord/ツリー/DAG）。
  - 文献集計: `GenePriority_with_literature.xlsx`（文献件数と総合スコア）。
  - セット比較: PCA/Dendrogram/ヒートマップ PNG、Euler 図 PDF/PNG、交差リスト Excel。
  - snAtlas 再解析: WGCNA 結果 CSV（モジュール相関・Hub 遺伝子・全遺伝子割り当て）、GO 結果 CSV/PDF、擬似バルク CPM 棒グラフ PDF、ドット図/棒グラフ PDF。

## 実行例
```bash
Rscript "3章_IP-MSによるPPI解析/GO analysis.R"
Rscript "3章_IP-MSによるPPI解析/related paper search.R"
Rscript "4章 _AirIDによるPPI解析/PCA分析.R"
Rscript "4章 _AirIDによるPPI解析/ベン図.R"
Rscript "4章 _AirIDによるPPI解析/ベン図_リスト.R"
Rscript "5章_snAtlas再解析/batch_WGCNA.R"
Rscript "5章_snAtlas再解析/batch_WGCNA_GO analysis.R"
Rscript "5章_snAtlas再解析/器官別発現量解析.R"
Rscript "5章_snAtlas再解析/細胞型別発現量解析-ドット図.R"
Rscript "5章_snAtlas再解析/細胞型別発現量解析-棒グラフ.R"
```
必要に応じて各スクリプト内の閾値（例: `PV_CUTOFF`, `TOP_PCT`, 背景モード `BG_MODE`）を調整してください。
