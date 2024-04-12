## Title: 09_PCA
## Author: Jules Larke
## Date: 032024
## Purpose: Perform PCA and visualize glycan features in euclidean space

##Load libraries
library(ggplot2)
library(dplyr)
library(cowplot)
library(vegan)

##Load data sets
free = read.csv('../data/glycan/free_saccharide_v2.csv')
poly = read.csv('../data/glycan/polysaccharide_v2.csv')
mono = read.csv('../data/glycan/monosaccharide_v2.csv')
link = read.csv('../data/glycan/glycosidic_linkage_v2.csv')

free = free %>% arrange(sample_id)
poly = poly %>% arrange(sample_id)
mono = mono %>% arrange(sample_id)
link = link %>% arrange(sample_id)

free$sample_id = NULL
free$simple_name = NULL
poly$sample_id = NULL
poly$simple_name = NULL
poly$food_group = NULL
mono$sample_id = NULL
mono$simple_name = NULL
mono$food_group = NULL
link$sample_id = NULL
link$simple_name = NULL
link$food_group = NULL

colnames(free)[2:13] <- paste0('free_',colnames(free)[2:13])
colnames(mono)[1:10] <- paste0('mono_',colnames(mono)[1:10])
colnames(link) <- c('t-Glucose',	'4-Glucose',	'6-Glucose',	'3-Glucose/3-Galactose',	'2-Glucose',	'4,6-Glucose',	'3,4-Glucose',	'2,4-Glucose',	'3,4,6-Glucose',	't-Galactose',	'6-Galactose',	'4-Galactose',	'2-Galactose',	'4,6-Galactose',	'3,6-Galactose',	'3,4-Galactose',	't-p-Xylose',	'4-p-Xylose',	'3-Xylose', 	'2-Xylose',	'3,4-P-Xylose/3,5-Arabinose',	'2,4-p-Xylose',	't-f-Arabinose',	't-p-Arabinose',	'5-f-Arabinose',	'3-Arabinose',	'2-f-Arabinose',	'2,3-f-Arabinose',	't-Fucose',	't-Rhamnose',	'4-Rhamnose', 	'2-Rhamnose',	'2,4-Rhamnose',	'4-Glucosamine/GlcNac',	'3-Glucosamine/GlcNac',	't-Mannose',	'4-Mannose',	'3-Mannose',	'2-Mannose',	'4,6-Mannose',	'2,X-Mannose',	'3,4,6-Mannose',	'X-Hexose',	'2,X,X-Hexose (I)',	'2,X,X-Hexose (II)',	'2,X,X-Hexose (III)',	't-Deoxyhexose')

#names.change <- grep("^[A-Z]+$", names(df)) # remove preceding X for colnames with numeric starting chacacter

all = cbind(free, poly, mono, link)
all$food_group <- as.factor(all$food_group)

#############################################################################
## PCAs
#############################################################################

sum(is.na(all))

pca <- prcomp(all[c(2:ncol(all))], scale. = T, center = T)
summary(pca)

## Create a scree plot to decide how many components to keep.
## screeplot(pca)
screeplot(pca, type = "line", main = "Scree Plot")

## Create a new dataframe with the metadata, along with the values for the first 2 PCs in each.
PC12 <- data.frame(all$food_group, pc1 = pca$x[,1], pc2 = pca$x[,2])
##View(PC1234)

## Variation explained by each PC. This is the same as we see in the summary(pca)
percent.var <- (pca$sdev)^2/sum(pca$sdev^2)

## Create a data frame with loadings for PCs 1-2
loadings <- data.frame(L1= pca$rotation[,1], L2= pca$rotation[,2])
loadings

##############################################################################
## PC1 vs PC2

## Loadings Plot  
plot(loadings$L1, loadings$L2, xlab= "PC1", cex.lab = 1.5, ylab="PC2", pch = 16, cex = 1.5)+
  text (x=loadings[,1], y=loadings[,2], 
        labels=row.names (loadings), 
        pos=3,
        cex=1)

##Centroids

cent <- aggregate(cbind(pc1, pc2) ~ all.food_group, data = PC12, FUN = mean)
g1 <- merge(PC12,aggregate(cbind(mean.x=pc1,mean.y=pc2)~all.food_group,PC12,mean),by="all.food_group")


## Scores: all features legend
scores.PC12 <- ggplot(g1, aes(x=pc1, y=pc2, color=factor(all.food_group))) +
  geom_point(size=3) +
  #geom_point(aes(x=mean.x,y=mean.y),size=5) +
  #geom_segment(aes(x=mean.x, y=mean.y, xend=pc1, yend=pc2), size = 1) +
  xlab(paste("PC1 (",round(percent.var[1]*100, digits=1), "%)", sep="")) +
  ylab(paste("PC2 (",round(percent.var[2]*100, digits=1), "%)", sep="")) + 
  scale_colour_manual(values = c('#ff7f0e', '#7f7f7f', '#bcbd22', '#d62728', "#1f77b4","#9467bd", "#e377c2", "#8c564b", "#2ca02c")) +
  ggtitle("")+
  #annotate("text", x = -4, y = 5, label = "italic(R) ^ 2 == 0.177",
  #         parse = TRUE, size = 6) +
  #annotate("text", x = -4, y = 4.5, label = "italic(p) < 0.001", parse = TRUE, size = 6) +
  stat_ellipse(geom = "polygon", size = 1.2, aes(fill = all.food_group), alpha = 0.25, show.legend = FALSE, type = "t", level = 0.95) +
  theme_bw() +
  theme(axis.line.x = element_line(color="black", size = .5),
        axis.line.y = element_line(color="black", size = .5),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 10),
        panel.border = element_blank(),
        panel.background = element_blank()) +
        #legend.position = "bottom")  +        
        scale_fill_manual(values = c('#ff7f0e', '#7f7f7f', '#bcbd22', '#d62728', "#1f77b4","#9467bd", "#e377c2", "#8c564b", "#2ca02c"))
scores.PC12
ggsave("dfg2_all_features_PCA.png",
       width = 9,
       height = 7,
       dpi=300)
#leg = get_legend(scores.PC12)

##################### TEST HOMOGENEITY AND EFFECT SIZE (ADONIS)
pcadisp <- vegdist(all[2:ncol(all)], method = "euclidean")

disper.all <- betadisper(pcadisp, all$food_group)

anova(disper.all)

plot(disper.all)

##reject null hypothesis; the groups have different dispersions; cannot test with permanova
#adtest <- adonis(free[2:ncol(free)]~free$food_group, data=free[2:ncol(free)], method = "euclidean", permutations = 999)
#adtest

# loadings plot 
library(ggrepel)
## Scores: all features
loading.PC12 <- ggplot(loadings, aes(x=L1, y=L2)) +
  geom_point(size=3) +
  xlab("PC1") +
  ylab("PC2") + 
  geom_text_repel(
    data = loadings,
    aes(x = L1, y = L2, label = rownames(loadings)),
    #force = 5,
    min.segment.length = 0.0001,
    max.overlaps = 20,
    cex = 3,
    color = "black",
  ) +
  theme_bw() +
  theme(axis.line.x = element_line(color="black", size = .5),
        axis.line.y = element_line(color="black", size = .5),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 10),
        panel.border = element_blank(),
        panel.background = element_blank())
loading.PC12
ggsave("dfg2_all_features_loadings_PCA.png",
       width = 7,
       height = 7,
       dpi=300)

all_dfg2 = cowplot::plot_grid(scores.PC12, loading.PC12, ncol = 2)
all_dfg2

