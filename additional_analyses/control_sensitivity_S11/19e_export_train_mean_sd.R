cat("PDK2,PDK4 rows from GSE96804 training cohort\n")
e <- readRDS("data/GSE96804_expr_symbol.rds")
stopifnot(all(c("PDK2","PDK4") %in% rownames(e)))
m <- e[c("PDK2","PDK4"), , drop = FALSE]
mu <- apply(m, 1, mean); sd <- apply(m, 1, sd)
out <- data.frame(gene = c("PDK2","PDK4"), mean_training = round(mu, 6), sd_training = round(sd, 6),
                  n_samples = ncol(m))
if (!dir.exists("results/revision")) dir.create("results/revision", recursive = TRUE)
write.csv(out, "results/revision/train_pdk2_pdk4_mean_sd.csv", row.names = FALSE)
print(out)
cat("TRAINSD_DONE\n")
