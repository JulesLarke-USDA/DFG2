## Title: 08_top_food_contributors
## Author: Jules Larke
## Date: 031524
## Purpose: visualize top food/ingredient contributors by amount of raffinose and xyloglucan
## intake for total and non-alcohol consumers respectively

library(tidyverse)
library(RColorBrewer)
library(cowplot)

data = read.csv('../output/02/ingred_source_top_glycans.csv')
outcome <- read_csv("../output/02/outcome_homa.csv")
data = left_join(data, outcome, by='SEQN')

######### Get sum total of each monosaccharide for unique foods/ingredients. Simple name is the glycopedia food name
data2 <- data %>%
  group_by(SEQN, simple_name) %>%
  summarise(raffinose_sum = sum(Free_Raffinose)
  )

data3 <- data2 %>%
  group_by(simple_name) %>%
  summarise(
    raffinose_mean = mean(raffinose_sum),
    raffinose_sem = (sd(raffinose_sum) / sqrt(13661))
  )

raff = data3 %>% arrange(desc(raffinose_mean))
raff = raff[1:5,]

raff$simple_name[1] = "Beer, lager"
raff$simple_name[2] = "Lentils"
raff$simple_name[3] = "Refried beans"
raff$simple_name[4] = "Cashews"
raff$simple_name[5] = "Pine nuts"

#labs <- c("IR", "non-IR")
lab1 <- expression(bgroup("",frac('Raffinose (g)','day • 1000 kcal'),""))


raff_plot <- ggplot(data = raff,
       aes(
         x = reorder(simple_name, -raffinose_mean),
         y = raffinose_mean,
         fill = factor(simple_name),
         #color = factor(simple_name)
       )) +  geom_errorbar(aes(ymin=raffinose_mean - raffinose_sem, ymax=raffinose_mean + raffinose_sem),
                position=position_dodge(0.4), size=0.2, width=0.1, color='black') +
  stat_summary(
    fun = mean,
    color = "black",
    geom = "bar",
    position = 'dodge',
    width = 0.4,
    size = 0.2
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "", y = lab1) +
  scale_fill_manual(values = c("darkorange", "darkorange", 'darkorange', 'darkorange', 'darkorange')) +
  #scale_fill_manual(values = c("seagreen3", "darkorchid3")) +
  theme_bw() +
  theme(
    title = element_text(size = 10),
    axis.line.x = element_line(color = "black", size = .5),
    axis.line.y = element_line(color = "black", size = .5),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 10, color = "black", vjust = 0.2),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.position = 'none'
  )
raff_plot
ggsave("../plots/raffinose_top_5.png", width = 6, height = 4, units = "in", dpi=1000)

# create boxplots for figure
labs <- c("IR", "non-IR")
lab1 <- expression(bgroup("",frac('Raffinose (g)','day • 1000 kcal'),""))

# compare mean intake on IR vs Non-IR
per_pid = data2 %>% group_by(SEQN, ir) %>% summarise(per_pid = sum(raffinose_sum))
per_pid$ir_name = ifelse(per_pid$ir > 0, 'IR', 'non-IR')
t.test(per_pid~ir_name, per_pid)

raff_by_ir <- ggplot(data = per_pid,
       aes(x = ir_name,
           y = per_pid,
           fill = ir_name)) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.1,
    size = 0.2
  ) +  stat_summary(
    fun = mean,
    color = "black",
    geom = "bar",
    width = 0.4,
    size = 0.2
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "", y = lab1) +
  scale_fill_manual(values = c("seagreen3", "darkorchid3")) +
  annotate(
    "segment",
    x = 1,
    xend = 2,
    y = 0.04,
    yend = 0.04,
    size = 0.5
  ) +
  annotate(
    "text",
    x = 1.5,
    y = 0.041,
    label = "***",
    size = 6
  ) +
  theme_bw() +
  theme(
    title = element_text(size = 8),
    axis.line.x = element_line(color = "black", size = .5),
    axis.line.y = element_line(color = "black", size = .5),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.ticks.x = element_blank(),
    legend.title = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #legend.position = "none",
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.position = 'none'
  )
raff_by_ir
ggsave("../plots/raffinose_by_IR.png", width = 2, height = 3, units = "in", dpi=1000)

# load in data for participants with zero alcohol intake
no_alc = read.csv('../output/02/homa_no_alc_ingred.csv')
no_alc_outcome <- read_csv("../output/02/homa_outcome_no_alc.csv")

no_alc = left_join(no_alc, no_alc_outcome, by='SEQN')

######### Get sum total of each monosaccharide for unique foods/ingredients. Simple name is the glycopedia food name
data2 <- no_alc %>%
  group_by(SEQN, simple_name) %>%
  summarise(xylo_sum = sum(Xyloglucan)
  )

data3 <- data2 %>%
  group_by(simple_name) %>%
  summarise(
    xylo_mean = mean(xylo_sum),
    xylo_sem = (sd(xylo_sum) / sqrt(10260))
  )

xylo = data3 %>% arrange(desc(xylo_mean))
xylo = xylo[1:5,]

xylo$simple_name[1] = "Cannellini beans"
xylo$simple_name[2] = "Apple, raw"
xylo$simple_name[3] = "Banana"
xylo$simple_name[4] = "Mango"
xylo$simple_name[5] = "Wheat flour"


#labs <- c("IR", "non-IR")
lab1 <- expression(bgroup("",frac('Xyloglucan (g)','day • 1000 kcal'),""))


xylo_plot <- ggplot(data = xylo,
                    aes(
                      x = reorder(simple_name, -xylo_mean),
                      y = xylo_mean,
                      fill = factor(simple_name),
                      #color = factor(simple_name)
                    )) +  geom_errorbar(aes(ymin=xylo_mean - xylo_sem, ymax=xylo_mean + xylo_sem),
                                        position=position_dodge(0.4), size=0.2, width=0.1, color='black') +
  stat_summary(
    fun = mean,
    color = "black",
    geom = "bar",
    position = 'dodge',
    width = 0.4,
    size = 0.2
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "", y = lab1) +
  scale_fill_manual(values = c("darkorange", "darkorange", 'darkorange', 'darkorange', 'darkorange')) +
  #scale_fill_manual(values = c("seagreen3", "darkorchid3")) +
  theme_bw() +
  theme(
    title = element_text(size = 10),
    axis.line.x = element_line(color = "black", size = .5),
    axis.line.y = element_line(color = "black", size = .5),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 10, color = "black", vjust = 0.2),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.position = 'none'
  )
xylo_plot
ggsave("../plots/xyloglucan_top_5.png", width = 6, height = 4, units = "in", dpi=1000)

per_pid = data2 %>% group_by(SEQN, ir) %>% summarise(per_pid = sum(xylo_sum))
per_pid$ir_name = ifelse(per_pid$ir > 0, 'IR', 'non-IR')
t.test(per_pid~ir_name, per_pid)

ggplot(data = per_pid,
       aes(x = ir_name,
           y = per_pid,
           fill = ir_name)) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.1,
    size = 0.2
  ) +  stat_summary(
    fun = mean,
    color = "black",
    geom = "bar",
    width = 0.4,
    size = 0.2
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "", y = lab1) +
  scale_fill_manual(values = c("seagreen3", "darkorchid3")) +
  annotate(
    "segment",
    x = 1,
    xend = 2,
    y = 0.45,
    yend = 0.45,
    size = 0.5
  ) +
  annotate(
    "text",
    x = 1.5,
    y = 0.46,
    label = "***",
    size = 6
  ) +
  theme_bw() +
  theme(
    title = element_text(size = 8),
    axis.line.x = element_line(color = "black", size = .5),
    axis.line.y = element_line(color = "black", size = .5),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.ticks.x = element_blank(),
    legend.title = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #legend.position = "none",
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.position = 'none'
  )

ggsave("../plots/xyloglucan_by_IR.png", width = 2, height = 3, units = "in", dpi=1000)
