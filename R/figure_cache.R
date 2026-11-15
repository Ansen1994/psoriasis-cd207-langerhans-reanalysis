## ============================================================================
## figure_cache.R - 图件输入缓存
##
## 目的：06_figures.R 单次 76 s，其中 33 s 花在 readRDS 一个 84,440 细胞的
## Seurat 对象上，而绘图只用到它的四样东西（UMAP、12 个标记基因、metadata
## 三列、细胞数）。缓存那个原对象毫无意义（读同样大的副本一样慢），
## 所以这里缓存的是「投影」：从大对象里抽出来的小对象，约 20 MB，读取 < 1 s。
##
## 用法：
##   s1 <- fig_cached("fig1_s1_projection", deps = <源文件路径>, build = function() {...})
##
## 失效规则（三重）：源文件路径 + 大小 + mtime 任一变化即重建。
## 手动清缓存：rm -rf data/cache
## ============================================================================

CACHE_DIR <- "data/cache"

fig_cache_path <- function(key) file.path(CACHE_DIR, paste0(key, ".rds"))

## 依赖签名：文件名 + 字节数 + mtime；任一变化就视为过期。
fig_sig <- function(deps) {
  deps <- deps[file.exists(deps)]
  if (!length(deps)) return("NO_DEPS")
  paste(basename(deps), file.size(deps),
        format(file.mtime(deps), "%Y%m%d%H%M%S"), collapse = "|")
}

fig_cached <- function(key, deps, build, verbose = TRUE) {
  path <- fig_cache_path(key)
  meta <- paste0(path, ".sig")
  sig  <- fig_sig(deps)

  if (file.exists(path) && file.exists(meta) &&
      identical(readLines(meta, warn = FALSE), sig)) {
    t0 <- Sys.time()
    v  <- readRDS(path)
    if (verbose)
      cat(sprintf("  CACHE HIT  %-26s (%.2f s)\n", key,
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    return(v)
  }

  if (verbose) cat(sprintf("  CACHE MISS %-26s building...\n", key))
  t0 <- Sys.time()
  v  <- build()
  dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
  saveRDS(v, path)
  writeLines(sig, meta)
  if (verbose)
    cat(sprintf("  CACHE SAVED %-25s (%.1f s, %.1f MB on disk)\n", key,
                as.numeric(difftime(Sys.time(), t0, units = "secs")),
                file.size(path) / 1024^2))
  v
}

## 只重绘某一张图，用于改版式时快速迭代：
##   FIG_ONLY=Fig3 Rscript R/06_figures.R
## 注意：跳过的那几张不会重写 PDF/PNG，最终定稿必须不带 FIG_ONLY 跑一次。
fig_should_render <- function(stem) {
  only <- Sys.getenv("FIG_ONLY", "")
  if (!nzchar(only)) return(TRUE)
  hit <- grepl(only, stem, ignore.case = TRUE)
  if (!hit) cat("  SKIP  ", stem, " (FIG_ONLY=", only, ")\n", sep = "")
  hit
}
