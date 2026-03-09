suppressPackageStartupMessages(library(ggplot2))

a <- commandArgs(trailingOnly = TRUE)
xT <- a[1]
mT <- a[2]
out <- a[3]

dir.create(file.path(out, "plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out, "tables"), showWarnings = FALSE, recursive = TRUE)

x <- read.table(xT, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
m <- read.table(mT, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

sid_x <- x$id[1]
sid_m <- m$id[1]

snv <- data.frame(
  id = x$id,
  FP = x$ct_snv_fp / x$cb,
  pseudoTP = x$ct_snv_tp / x$cb
)
write.table(
  snv,
  file.path(out, "tables", paste0(sid_x, ".snv_xy.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
p1 <- ggplot(snv, aes(FP, pseudoTP)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", formula = y ~ 0 + x, se = FALSE) +
  theme_classic()
ggsave(
  file.path(out, "plots", paste0(sid_x, ".snv_tp_vs_fp.pdf")),
  p1,
  width = 6,
  height = 4,
  bg = "transparent"
)

ind <- data.frame(
  id = x$id,
  FP = x$ct_ind_fp / x$cb,
  pseudoTP = x$ct_ind_tp / x$cb
)
write.table(
  ind,
  file.path(out, "tables", paste0(sid_x, ".ind_xy.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
p2 <- ggplot(ind, aes(FP, pseudoTP)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", formula = y ~ 0 + x, se = FALSE) +
  theme_classic()
ggsave(
  file.path(out, "plots", paste0(sid_x, ".ind_tp_vs_fp.pdf")),
  p2,
  width = 6,
  height = 4,
  bg = "transparent"
)

mix <- data.frame(
  id = m$id,
  hsnp = 2 * m$gHetSNV / m$cb,
  mixAdj = (2 * m$gHetSNV / m$cb) * m$refLen / m$pureS,
  pureS = m$pureS
)
write.table(
  mix,
  file.path(out, "tables", paste0(sid_m, ".mix.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
p3 <- ggplot(mix, aes(hsnp, mixAdj)) +
  geom_point(size = 3) +
  theme_classic()
ggsave(
  file.path(out, "plots", paste0(sid_m, ".mix_inferred.pdf")),
  p3,
  width = 6,
  height = 4,
  bg = "transparent"
)

prec_snv <- (x$ct_snv_tot - x$ct_snv_fp) / x$ct_snv_tot
prec_ind <- (x$ct_ind_tot - x$ct_ind_fp) / x$ct_ind_tot

exp_snv <- x$gHetSNV / (2 * x$refLen)
exp_ind <- x$gHetIND / (2 * x$refLen)

obs_snv <- x$ct_snv_tp / x$cb
obs_ind <- x$ct_ind_tp / x$cb

c_snv <- prec_snv / (obs_snv / exp_snv)
c_ind <- prec_ind / (obs_ind / exp_ind)

bur <- data.frame(
  id = x$id,
  snv = (x$somSNV / x$cb) * c_snv,
  ind = (x$somIND / x$cb) * c_ind,
  c_snv = c_snv,
  c_ind = c_ind
)
write.table(
  bur,
  file.path(out, "tables", paste0(sid_x, ".burden.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

p4 <- ggplot(bur, aes(id, snv)) +
  geom_point(size = 3) +
  theme_classic()
ggsave(
  file.path(out, "plots", paste0(sid_x, ".burden_snv.pdf")),
  p4,
  width = 5,
  height = 3.5,
  bg = "transparent"
)

p5 <- ggplot(bur, aes(id, ind)) +
  geom_point(size = 3) +
  theme_classic()
ggsave(
  file.path(out, "plots", paste0(sid_x, ".burden_ind.pdf")),
  p5,
  width = 5,
  height = 3.5,
  bg = "transparent"
)

cat("OK\n")
