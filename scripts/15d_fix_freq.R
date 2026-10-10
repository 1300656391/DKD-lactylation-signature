# 15d_fix_freq.R ----
# 修正M4选择频率分母（原/168可>1；改为计数与/200两列，定义明确）
library(data.table)
sel <- fread("results/revision/M4_selection_frequency.csv")
sel$count <- round(sel$selection_freq * 168)
sel$freq_of_200 <- sel$count / 200
setorder(sel, -count)
fwrite(sel, "results/revision/M4_selection_frequency.csv")
file.copy("results/revision/M4_selection_frequency.csv",
          "manuscript/supplement/TableS9b_selection_frequency.csv", overwrite = TRUE)
print(sel[, .(gene, count, freq_of_200)][1:8])
