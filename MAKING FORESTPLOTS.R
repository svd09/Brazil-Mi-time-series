#---------------------------#
# FORESTPLOTS FOR THE PAPER #
#---------------------------#
# salil v deo
# 2025-09-14
# ckbplotr package to make forest plots for the paper.
rm(list = ls())
library(easypackages)
libraries(c("tidyverse","ckbplotr",
            "data.table","janitor",
            "stringr"))
#---------------#
# PLOTS FOR EHE #
#---------------# 
data = read_csv("results/ehe_by_location_grouped.csv")
glimpse(data)
data$` ` <- paste(rep(" ", 20), collapse = " ")
data$location <- ifelse(is.na(data$location), "", data$location)
data$variable <- ifelse(is.na(data$variable), "", data$variable)
data$OR = paste0(data$est," (",data$low," -",data$high," )")
data$OR <- ifelse(is.na(data$OR), "", data$OR)
res <- c(data$OR)

num = c(1,5,9,13,17)
for(k in num){
  res[k] <- " "
}
data$OR = res
# format the data.
locations = c(rep("Brazil",3),rep("North",3),rep("Northeast",3),rep("Southeast",3),rep("South",3))
events = c("EHE_90","EHE_90_2","EHE_90_3")
colors = c("#0072B2","#CC79A7","#D55E00")
df_values = data.frame(
  est = data$est,
  low = data$low,
  high = data$high
)

df_values = df_values %>%
  filter(!is.na(est)) %>%
  mutate(
    key = 1:15)
df_values = df_values %>%
  mutate(color = rep(colors,5)) 

df_labels = data.frame(
  group = locations,
  subgroup = rep(events,5),
  key = 1:nrow(df_values))


# make table to save this in the folder.
df_table = cbind(df_values,df_labels)
write_csv(df_table,"results/forest_plot_ehe_table.csv")
# make plot
p1 <- forest_plot(panels = list(df_values),
            row.labels = df_labels,
            col.estimate = 'est',
            col.lci = 'low',
            col.uci = 'high',
            digits = 3,
            colour = "color",
            exponentiate = F,
            xlim = c(0.95,1.2),
            col.right.heading = "Odds Ratio (95%CI)",
            xlab = "Odds Ratio (95%CI)",
            row.labels.heading = "Regions") 

# add vertical line.
p2 <- p1$plot + geom_vline(xintercept = 1, linetype="dashed")
# save this plot
ggsave(p2,
       filename = "graphs/forest_plot_ehe_all_locations.tiff",
       width = 6,
       height = 7,
       dpi = 600,device = "tiff")

# plot for the EHE_99 events.
# get the data and then reformat to get into the correct format.
#------------------#
# PLOT FOR EHE_99. #
#------------------#

df = read_csv("results/ehe_estimates_long.csv")
glimpse(df)
df2 = df %>% filter(variable %in% c("EHE_99","EHE_99_2","EHE_99_3"))
df2 = df2 %>% arrange(location,variable)
colors = c("#0072B2","#CC79A7","#D55E00")
df2$color = rep(colors,5)
df2$key = 1:nrow(df2)
locations = c(rep("Brazil",3),rep("North",3),rep("Northeast",3),rep("Southeast",3),rep("South",3))
events = c("EHE_99","EHE_99_2","EHE_99_3")
df_labels = data.frame(
  group = locations,
  subgroup = rep(events,5),
  key = 1:nrow(df2))

p99_1 <- forest_plot(panels = list(df2),
                  row.labels = df_labels,
                  col.estimate = 'est',
                  col.lci = 'low',
                  col.uci = 'high',
                  digits = 3,
                  colour = "color",
                  exponentiate = F,
                  xlim = c(0.95,1.4),
                  col.right.heading = "Odds Ratio (95%CI)",
                  xlab = "Odds Ratio (95%CI)",
                  row.labels.heading = "Regions") 

# add vertical line.
p99_2 <- p99_1$plot + geom_vline(xintercept = 1, linetype="dashed")
p99_2
# use patchwork to get the panel of both plots together.
library(cowplot)
combined_plot <- cowplot::plot_grid(p2,p99_2,
                   labels = c("A","B"),
                   ncol = 2,
                   align = "h")
# save this plot
ggsave(combined_plot,
       filename = "graphs/forest_plot_ehe.tiff",
       width = 12,
       height = 7,
       dpi = 600,device = "tiff")

#-----------#
# ECE PLOTS #
#-----------#

data = read_csv('results/ece_rr_table.csv')
glimpse(data)
data = data %>% rename(
  Location = location,
  Event = variable
)
# limit to ECE_10
data1 = data %>% filter(Event %in% c("ECE_10","ECE_10_2","ECE_10_3"))
# now to make key data.
data1$key = 1:nrow(data1)
colors = c("#0072B2","#CC79A7","#D55E00")
data1$color = rep(colors,5)
# naming.
names = data.frame(
  Location = c(rep("Brazil",3),rep("North",3),rep("Northeast",3),rep("Southeast",3),rep("South",3)),
  Event = rep(c("ECE_10","ECE_10_2","ECE_10_3"),5),
  key = 1:nrow(data1)
)
p1 <- forest_plot(panels = list(data1),
                  row.labels = names,
                  col.estimate = 'est',
                  col.lci = 'low',
                  col.uci = 'high',
                  digits = 3,
                  colour = "color",
                  exponentiate = F,
                  xlim = c(0.95,1.4),
                  col.right.heading = "Odds Ratio (95%CI)",
                  xlab = "Odds Ratio (95%CI)",
                  row.labels.heading = "Regions") 

# add vertical line.
p2 <- p1$plot + geom_vline(xintercept = 1, linetype="dashed")
p2

#- ECE_01

data2 = data %>% filter(Event %in% c("ECE_01","ECE_01_2","ECE_01_3"))
# now to make key data.
data2$key = 1:nrow(data1)
colors = c("#0072B2","#CC79A7","#D55E00")
data2$color = rep(colors,5)
# naming.
names = data.frame(
  Location = c(rep("Brazil",3),rep("North",3),rep("Northeast",3),rep("Southeast",3),rep("South",3)),
  Event = rep(c("ECE_10","ECE_10_2","ECE_10_3"),5),
  key = 1:nrow(data1)
)
p_1 <- forest_plot(panels = list(data2),
                  row.labels = names,
                  col.estimate = 'est',
                  col.lci = 'low',
                  col.uci = 'high',
                  digits = 3,
                  colour = "color",
                  exponentiate = F,
                  xlim = c(0.95,2.5),
                  col.right.heading = "Odds Ratio (95%CI)",
                  xlab = "Odds Ratio (95%CI)",
                  row.labels.heading = "Regions") 

# add vertical line.
p_2 <- p_1$plot + geom_vline(xintercept = 1, linetype="dashed")
p_2
library(cowplot)
p_comb = cowplot::plot_grid(p2,p_2,
                   labels = c("A","B"),
                   ncol = 2,
                   align = "h")
ggsave(p_comb,
       filename = "graphs/forest_plot_ece.tiff",
       width = 12,
       height = 7,
       dpi = 600,device = "tiff")

#--- END SCRIPT ---=============================================================