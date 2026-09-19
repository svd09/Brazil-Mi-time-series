# --- Data (as provided) ---
df <- data.frame(
  Subgroup = rep(c("Males","Females","Young","Old"), each=6),
  Event = rep(c("EHE_90","EHE_90_2","EHE_90_3","EHE_99","EHE_99_2","EHE_99_3"), times=4),
  est = c(1.007,1.025,1.039,1.117,1.153,1.160,
          1.051,1.065,1.069,1.181,1.195,1.191,
          1.009,1.018,1.023,1.115,1.129,1.122,
          1.036,1.060,1.073,1.175,1.213,1.219),
  low = c(1.001,1.018,1.032,1.103,1.135,1.139,
          1.043,1.056,1.060,1.163,1.172,1.164,
          1.003,1.011,1.016,1.099,1.109,1.099,
          1.030,1.053,1.066,1.159,1.193,1.195),
  high = c(1.013,1.031,1.046,1.131,1.172,1.182,
           1.059,1.073,1.078,1.200,1.219,1.218,
           1.016,1.025,1.031,1.132,1.149,1.146,
           1.043,1.068,1.081,1.191,1.233,1.243)
)

# --- Plot settings ---
events <- levels(factor(df$Event))
cols   <- c("blue","red","darkgreen","purple")
pad    <- 0.01  # vertical padding above each CI line for the label

# Expand ylim to make space for vertical labels above the highest CI
ymax <- max(df$high) + 0.05
tiff("graphs/EHE_subgroup_plot.tiff", 
     width=10, height=9, 
     units="in", res=600)
par(mar=c(6,5,4,2))
plot(NA, xlim=c(0.5, length(events) + 0.5), ylim=c(1.0, ymax),
     xaxt="n", xlab="", ylab="Odds Ratio (95% CI)",
     frame = F)

axis(1, at=1:length(events), labels=events, las=2)
abline(h=1, lty=2, col="gray")

# --- Draw points, CI lines (no whiskers), and vertical labels ---

for (i in 1:4) {
  subname <- unique(df$Subgroup)[i]
  dsub <- df[df$Subgroup == subname, ]
  x <- 1:length(events) + (i - 2.5) * 0.15  # small horizontal offset by subgroup
  
  # points + CI lines
  points(x, dsub$est, pch=19, col=cols[i])
  segments(x, dsub$low, x, dsub$high, col=cols[i], lwd=1.5)
  
  # vertical labels: "est (low, high)" with 2 decimals, placed above the CI line
  labs <- sprintf("%.2f (%.2f, %.2f)", dsub$est, dsub$low, dsub$high)
  text(x, dsub$high + pad, labels=labs, srt=90, col=cols[i],
       cex=0.7, adj=c(0, 0.5))  # adj so text starts just above the line
}

legend("topleft", legend=unique(df$Subgroup), col=cols, pch=19, bty="n")
dev.off()