# --- Create data frame ---
df <- data.frame(
  Subgroup = rep(c("Males","Females",
                   "Young","Old"), each=6),
  Event = rep(c("ECE_10","ECE_10_2","ECE_10_3","ECE_01","ECE_01_2","ECE_02_3"), times=4),
  est = c(1.123,1.125,1.119,1.213,1.214,1.173,
          1.114,1.129,1.126,1.206,1.210,1.176,
          1.087,1.096,1.091,1.158,1.150,1.108,
          1.156,1.164,1.159,1.264,1.272,1.233),
  low = c(1.116,1.118,1.111,1.198,1.194,1.148,
          1.106,1.120,1.117,1.188,1.185,1.146,
          1.080,1.088,1.082,1.141,1.128,1.082,
          1.149,1.156,1.150,1.247,1.250,1.206),
  high = c(1.129,1.132,1.127,1.229,1.235,1.199,
           1.122,1.138,1.136,1.225,1.234,1.207,
           1.095,1.104,1.100,1.175,1.172,1.136,
           1.163,1.172,1.168,1.281,1.295,1.262)
)

# --- Plot setup ---
events <- levels(factor(df$Event))
cols   <- c("blue","red","darkgreen","purple")
pad    <- 0.01
ymax   <- max(df$high) + 0.07

tiff("graphs/ECE_effects_by_subgroup.tiff",
     width=8, height=8, res=600,units="in")
par(mar=c(6,5,4,2))
plot(NA, xlim=c(0.5, length(events)+0.5), ylim=c(1.0, ymax),
     xaxt="n", xlab="", ylab="Odds Ratio (95% CI)",
     frame = F)

axis(1, at=1:length(events), labels=events, las=2)
abline(h=1, lty=2, col="gray")

# --- Draw points, CI lines, and rotated labels ---
for (i in 1:4) {
  dsub <- df[df$Subgroup == unique(df$Subgroup)[i], ]
  x <- 1:length(events) + (i - 2.5) * 0.15
  
  # point estimates + CI lines
  points(x, dsub$est, pch=19, col=cols[i])
  segments(x, dsub$low, x, dsub$high, col=cols[i], lwd=1.5)
  
  # labels as "est (low, high)" vertically above each line
  labs <- sprintf("%.3f (%.3f, %.3f)", dsub$est, dsub$low, dsub$high)
  text(x, dsub$high + pad, labels=labs, srt=90, col=cols[i],
       cex=0.7, adj=c(0,0.5))
}

legend("topleft", legend=unique(df$Subgroup), col=cols, pch=19, bty="n")
dev.off()