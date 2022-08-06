library(R.matlab)
library(tidyverse)
library(ggplot2)
library(ggpubr)

figPath = "Z:/projects/ANLS Modeling/Figures/"

setwd("C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Single-cell Modeling/NeuronModeling/")

############################################
# Fig Half-lives
############################################

d = readMat("data/HalfLives.mat")
glycData = as.numeric(d$s[,,1]$glycVals)
mitoData = as.numeric(d$s[,,1]$mitoVals)


names = c("Glyc.","Mito.")

y = c(glycData, mitoData);
x = factor(c(rep(1,length(glycData)), rep(2,length(mitoData))), 1:2, names)

color_palette = c('#B5D39B','#6B97BC')  # light green, light blue

ds = tibble(x=x, y=y)

pHL = ggplot(ds, aes(x = x, y = y, fill = x)) +
  geom_violin(trim=F, show.legend=F, scale='count') +
  scale_fill_manual(values=color_palette) +
  theme_classic() + 
  ylab('Enzyme half-life (h)') +
  xlab('') +
  theme(text = element_text(size=14),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5,
                                   color='black', size=14),
        axis.text.y = element_text(color='black', size=14),
        axis.line.x = element_blank()) #+
  #ylim(c(0.1,0.38)) # +
pHL


ggsave(
  paste0(figPath, "FigHalfLife.png"),
  plot = pHL,
  width = 2.5, height = 3, dpi = 300)

ggsave(
  paste0(figPath, "FigHalfLife.svg"),
  plot = pHL,
  width = 2.5, height = 3, dpi = 300)
