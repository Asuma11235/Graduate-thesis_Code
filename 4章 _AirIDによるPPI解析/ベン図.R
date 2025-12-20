

# ============================================================
# 8シートの prot_acc を読み込み -> 8集合 -> Euler図（近似ベン図）を描く
#   - AirID: 暖色系
#   - CoIP : 寒色系
# ============================================================

IN_XLSX_AIRID <- "/Users/ryo/Downloads/251220_AirID.xlsx"
SHEETS_AIRID  <- c("DMSO","bikinin","mock","BAK1")

IN_XLSX_COIP  <- "/Users/ryo/Downloads/251108_CoIP.xlsx"
SHEETS_COIP   <- c("DMSO","bikinin","BR-up","BR-down")

OUT_PDF <- "/Users/ryo/Downloads/AirID_CoIP_Euler8setssym.pdf"
OUT_PNG <- "/Users/ryo/Downloads/AirID_CoIP_Euler8setssym.png"

# --- packages ---
pkgs <- c("readxl", "dplyr", "stringr", "eulerr")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

library(readxl)
library(dplyr)
library(stringr)
library(eulerr)

# --- helper: read prot_acc from a sheet ---
read_prot_acc_set <- function(xlsx, sheet, col = "SYMBOL") {
  df <- readxl::read_excel(xlsx, sheet = sheet)
  if (!col %in% colnames(df)) {
    stop(sprintf("Column '%s' not found in sheet '%s' of '%s'", col, sheet, xlsx))
  }
  x <- df[[col]] |>
    as.character() |>
    str_trim() |>
    na.omit()
  
  # 空文字除去
  x <- x[nzchar(x)]
  
  # 末尾の余計な注釈などが混じる場合の保険（必要なら調整）
  # x <- sub("\\s+.*$", "", x)  # 例: "P12345 something" -> "P12345"
  
  unique(x)
}

# --- make sets (named list) ---
sets_airid <- setNames(
  lapply(SHEETS_AIRID, function(sh) read_prot_acc_set(IN_XLSX_AIRID, sh)),
  paste0("AirID_", SHEETS_AIRID)
)

sets_coip <- setNames(
  lapply(SHEETS_COIP,  function(sh) read_prot_acc_set(IN_XLSX_COIP,  sh)),
  paste0("CoIP_", SHEETS_COIP)
)

sets_all <- c(sets_airid, sets_coip)

# --- quick sanity check ---
cat("=== set sizes ===\n")
print(sapply(sets_all, length))

# --- Fit Euler diagram (area is approximate) ---
# shape="ellipse" は見た目が安定しやすい
fit <- eulerr::euler(sets_all, shape = "ellipse")

# --- colors: warm for AirID, cool for CoIP ---
fill_cols <- c(
  # AirID (warm)
  "AirID_DMSO"    = "#F4A261",
  "AirID_bikinin" = "#E76F51",
  "AirID_mock"    = "#F6BD60",
  "AirID_BAK1"    = "#E9C46A",
  # CoIP (cool)
  "CoIP_DMSO"     = "#4D96FF",
  "CoIP_bikinin"  = "#1F77B4",
  "CoIP_BR-up"    = "#3FB6B2",
  "CoIP_BR-down"  = "#2A9D8F"
)

# 念のため、名前が合ってるかチェック（足りないとエラーになりうる）
missing_cols <- setdiff(names(sets_all), names(fill_cols))
if (length(missing_cols) > 0) {
  stop("fill_cols に色が定義されていない集合があります: ",
       paste(missing_cols, collapse = ", "))
}

# --- plot ---
# quantities = TRUE で各領域に個数表示
# legend は TRUE にすると右に出ます
p <- function() {
  plot(
    fit,
    fills = list(fill = fill_cols, alpha = 0.45),
    edges = list(col = "grey30", lwd = 1),
    labels = list(col = "grey10", cex = 0.9),
    quantities = list(col = "grey10", cex = 0.75),
    legend = list(cex = 0.8),
    main = "AirID (warm) vs CoIP (cool): prot_acc overlap (Euler diagram)"
  )
}


# --- save ---
pdf(OUT_PDF, width = 10, height = 7)
p()
dev.off()

png(OUT_PNG, width = 1800, height = 1200, res = 200)
p()
dev.off()

cat("\nSaved:\n", OUT_PDF, "\n", OUT_PNG, "\n")

