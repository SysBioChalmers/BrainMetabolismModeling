library(R.matlab)
library(tidyverse)
library(ggplot2)
library(ggpubr)

figPath = "Z:/projects/ANLS Modeling/Figures/"

setwd("C:/Work/MatlabCode/projects/BrainMetabolismModeling/BrainMetabolismModeling/")

#These were calculated with the model in MATLAB
mitoDivGlycWithMT = 13.9088
mitoDivGlycWithoutMT = 10.5520

############################
# Make a figure that shows how much more expensive it is to use
# mitochondria compared to glycolysis if they are static.
#################################

#values from the matlab modeling
#Just copied here from above for convenience
#mitoDivGlycWithMT = 13.9088
#mitoDivGlycWithoutMT = 10.5520


#scales are 1 for glyc, 13.9088 and 10.5520 for mito.
getEAMCA = function (scaleWithMT, scaleWithoutMT, trCost, util) {
  totalCost = (scaleWithMT + trCost*scaleWithoutMT)/util  - scaleWithMT
}

util = (1:100)*0.01
extraCostGlyc = rep(NA, length(util))
extraCostMito = rep(NA, length(util))
#extraCostMitoLim = rep(NA, 40)
extraCostMitoLim = rep(NA, length(util))
for (i in 1:length(util)) {
  #hardcode transport cost to 0%
  extraCostGlyc[i] = getEAMCA(1,1, 0, util[i])
  extraCostMito[i] = getEAMCA(mitoDivGlycWithMT, mitoDivGlycWithoutMT, 0, util[i])
}

for (i in 1:length(extraCostMitoLim)) {
  #hardcode transport cost to 0%
  #extraCostMitoLim[i] = getEAMCA(mitoDivGlycWithMT, mitoDivGlycWithoutMT, 0, 0.4)
  extraCostMitoLim[i] = getEAMCA(mitoDivGlycWithMT, mitoDivGlycWithoutMT, 0, util[i] + 0.4*(1-util[i]))
}

names = c("Glyc.", "Mit.", "Mit. mob.       " );
#allX = c(util,util,util[1:40]) 
allX = c(util,util,util) 
allY = c(extraCostGlyc, extraCostMito, extraCostMitoLim) 
curve = factor(c(rep(1,length(util)), rep(2,length(util)), rep(3,length(extraCostMitoLim))), 1:3, names)

ds = tibble(x=allX, y=log2(allY), curve = curve)

#cols = c('#95B37B','#E0B064','#6B97BC')
cols = c('#95B37B','#6B97BC','#6B97BC')
#y=expression(Log[2]*" fold change maint. cost")
pY = ggplot(ds, aes(x = x, y = y, colour = curve, linetype = curve)) +
  geom_line(size=1.3) +
  scale_linetype_manual(values = c(1,1,2), labels = names) +
  scale_color_manual(values = c(cols[1],cols[2],cols[3]), labels = names) +
  ggplot2::labs(y=expression("Log"[2]*"(EAMCA)"), x="Static utilization", colour="Pathway", linetype="Pathway") +
#  ylim(0,1) +
  ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=14, color="black"), axis.text = element_text(size=14, color="black"), legend.text = element_text(size=14, color="black"), legend.position = "bottom", legend.title=element_blank())
#  guides(colour = guide_legend(title.position = element_blank()),
#         linetype = guide_legend(title.position = element_blank()))
  #guides(colour = guide_legend(title.position = "top", title.hjust = 0.5),
  #       linetype = guide_legend(title.position = "top", title.hjust = 0.5))
  pY


ggsave(
  paste0(figPath, "MaintPerUtil.png"),
  plot = pY,
  width = 3, height = 3.2, dpi = 300)

ggsave(
  paste0(figPath, "MaintPerUtil.eps"),
  plot = pY,
  width = 3, height = 3.2, dpi = 300)

###################################
# Similar plot for transport cost
###################################

util2 = 1
transp = (1:200)*0.01
extraCostGlyc = rep(NA, length(transp))
extraCostMito = rep(NA, length(transp))
for (i in 1:length(transp)) {
  extraCostGlyc[i] = getEAMCA(1,1, transp[i], util2)
  extraCostMito[i] = getEAMCA(mitoDivGlycWithMT, mitoDivGlycWithoutMT, transp[i], util2)
}

names = c("Glyc.", "Mito. resp.");
allX = c(transp,transp) 
allY = c(extraCostGlyc, extraCostMito) 
curve = factor(c(rep(1,length(transp)), rep(2,length(transp))), 1:2, names)

ds = tibble(x=allX, y=log2(allY), curve = curve)

cols = c('#95B37B','#6B97BC')
pY = ggplot(ds, aes(x = x, y = y, colour = curve)) +
  geom_line(size=1.3) +
  scale_color_manual(values = c(cols[1],cols[2],cols[3]), labels = names) +
  ggplot2::labs(y=expression("Log"[2]*"(EAMCA)"), x="Transportation cost", colour="Pathway") +
  ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=14, color="black"), axis.text = element_text(size=14, color="black"), legend.text = element_text(size=14, color="black"), legend.position = "bottom", legend.title = element_blank())
pY


ggsave(
  paste0(figPath, "MaintPerTransp.png"),
  plot = pY,
  width = 3.1, height = 3.2, dpi = 300)

ggsave(
  paste0(figPath, "MaintPerTransp.eps"),
  plot = pY,
  width = 3.1, height = 3.2, dpi = 300)



###########################
#Now a  plot showing optimal usage (glyc)
#neurons and astrocytes when lactate export is limited
###########################

#calculate the cost for both glyc and mito if we assume mito util >= 0.3 and transport = 0.1

util = (1:1000)*0.001
mitoExtraCost = rep(NA, length(util))
cytoExtraCost = rep(NA, length(util))

for (i in 1:length(util)) {
  #Transportation cost set to 0.1
  mitoExtraCost[i] = getEAMCA(mitoDivGlycWithMT, mitoDivGlycWithoutMT, 0.1, util[i] + 0.4*(1-util[i]))
  cytoExtraCost[i] = getEAMCA(1,1, 0.1, util[i])
}

#check at which positions glycolysis is more optimal. The lower switch can be calculated as 
#the saved maintenance cost per lactate output, where lactate output is just represented by

#check at which positions glycolysis is more optimal. The lower switch can be calculated as 
#the saved maintenance cost per lactate output, where lactate output is just represented by
#the static utilization, since the enzyme is inactive the rest of the time.
diff = mitoExtraCost - cytoExtraCost
sel = diff >= diff[length(util)]
gainY = diff[sel]

firstEdge = 0.0529
secEdge = 0.0889

allX = c(firstEdge,util[sel]) #firstEdge is for making the curve go to the bottom, the x value - manual interpolation
allY = c(min(gainY),gainY) #min(gainY) is to find the bottom of the curve
allY[secEdge*1000]
ds = tibble(x=log2(allX + 0.10), y=log2(allY))

cols2 = c('#c0D0F8','#C0E9B0','#FBF0D0')
names2 = c('Mit.', 'Glyc.','Both')

maxY = 4;

d = tibble(x1=log2(c(0.01,firstEdge,secEdge,0.40) + 0.10), x2=log2(c(firstEdge, secEdge, 0.40, 1) + 0.10), y1=log2(rep(min(allY), 4)), y2=rep(maxY, 4), Optimal = factor(c(1,3,2,3),1:3,names2))
dLine1 = tibble(x=log2(c(secEdge,secEdge) + 0.10), y=c(min(ds$y), maxY))
dLine2 = tibble(x=log2(c(0.40,0.40) + 0.10), y=c(min(ds$y), maxY))
dLine3 = tibble(x=log2(c(firstEdge,firstEdge) + 0.10), y=c(min(ds$y), maxY))

pZ = ggplot() +
  geom_rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=Optimal), color=NA, size=0,alpha=1) +
  geom_line(data = ds,size=1.3, mapping=aes(x = x, y = y), color = "black") +
  geom_line(data=dLine1,aes(x = x, y=y), color = "black", linetype = "dashed") +
  geom_line(data=dLine2,aes(x = x, y=y), color = "black", linetype = "dashed") +
  geom_line(data=dLine3,aes(x = x, y=y), color = "black", linetype = "dashed") +
  scale_fill_manual(values = c(cols2[1],cols2[2],cols2[3]), labels = names2) +
  scale_x_continuous(breaks = log2(c(0,0.10,0.25,0.5,1) + 0.10), labels = c('0','0.10','0.25','0.5','1' )) +
  ggplot2::labs(y=expression("Log"[2]*"(EAMCA"[m]*" - EAMCA"[g]*")      "), x="Static utilization", colour="Pathway") +
  ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=14, color="black"), axis.text = element_text(size=14, color="black"), legend.text = element_text(size=14, color="black"), legend.position = "bottom", legend.box = "horizontal", plot.margin = margin(r = 13), legend.title.align = 0.5) +
  guides(fill = guide_legend(title.position = "top", order=1, direction="horizontal", title = "Optimal pathway:")) + #, title.justification='center')) + 
  ylim(log2(min(allY)),maxY)
pZ

ggsave(
  paste0(figPath, "OptimalUse.png"),
  plot = pZ,
  width = 5.25, height = 2.54, dpi = 300)

ggsave(
  paste0(figPath, "OptimalUse.eps"),
  plot = pZ,
  width = 2.8, height = 3.35, dpi = 300)




#################################
# Modeling of the full model - Utilization
#################################


library(R.matlab)

d = readMat("data/FullModelOutput.mat")$d[,,1] 


ATPProdGlycNeur = as.numeric(d$ATPProdGlycNeur)
ATPProdGlycAstr = as.numeric(d$ATPProdGlycAstr)
ATPProdMitoNeur = as.numeric(d$ATPProdMitoNeur)
ATPProdMitoAstr = as.numeric(d$ATPProdMitoAstr)

#Neurons
##############


NPureGlyc = ATPProdGlycNeur - ATPProdMitoNeur/14.5
NPureGlyc[NPureGlyc < 0] = 0
NPureMito = ATPProdMitoNeur - ATPProdGlycNeur*14.5
NPureMito[NPureMito < 0] = 0
NMix = ATPProdGlycNeur + ATPProdMitoNeur - NPureGlyc - NPureMito

allY = NPureGlyc + NPureMito + NMix
NPureGlyc = NPureGlyc / allY
NPureMito = NPureMito / allY
NMix = NMix / allY

#avoid slopes between two points by repeating all points 10 times and spreading out the X
NPGExp = rep(NPureGlyc, each=10)
NPMExp = rep(NPureMito, each=10)
NPBExp = rep(NMix, each=10)


df = tibble(x=log2(((1:1000) + 5)/1000 + 0.10) , Mito = NPMExp, Glyc = NPGExp, Both = NPBExp, zeros = rep(0,1000))


cols3 = c('#c0D0F8','#C0E9B0','#FBF0D0')

names3 = c("Mito.", "Glyc.", "Both")

pZ = ggplot(df, aes(x=x)) +
  geom_area(aes(y=Mito+Glyc+Both, fill="Mix")) + 
  geom_area(aes(y=Mito+Glyc, fill="Glyc")) + 
  geom_area(aes(y=Mito, fill="Mito")) +
  geom_line(aes(y=Mito+Glyc+Both)) + 
  geom_line(aes(y=Mito+Glyc)) + 
  geom_line(aes(y=Mito)) +
  geom_line(aes(y=zeros)) +
  scale_fill_manual(values = cols3, labels = names3, breaks=c('Mito', 'Glyc', 'Mix')) +
  scale_x_continuous(breaks = log2(c(0,0.10,0.25,0.5,1) + 0.10), labels = c('0','0.10','0.25','0.5','1' )) +
  scale_y_continuous(breaks = c(0,1), labels = c('0','1' )) +
  ggplot2::labs(y="ATP Frac.", x="Static utilization") +
  ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=14, color="black"), axis.text = element_text(size=14, color="black"), legend.text = element_text(size=14, color="black"), legend.position = "none", legend.box = "vertical", plot.margin = margin(t=20), axis.text.x = element_blank(), axis.title.x = element_blank()) +
  #guides(fill = guide_legend(order=1, direction="horizontal", title = element_blank()))
  annotate(geom="text", x=-1.54, y=0.5, label="Neurons", color="black", size=5)  #fontface=2
pZ


ggsave(
  paste0(figPath, "CombModel1.eps"),
  plot = pZ,
  width = 3, height = 1.13, dpi = 300)


#Astrocytes
##############

APureGlyc = ATPProdGlycAstr - ATPProdMitoAstr/14.5
APureGlyc[APureGlyc < 0] = 0
APureMito = ATPProdMitoAstr - ATPProdGlycAstr*14.5
APureMito[APureMito < 0] = 0
AMix = ATPProdGlycAstr + ATPProdMitoAstr - APureGlyc - APureMito

allY = APureGlyc + APureMito + AMix
APureGlyc = APureGlyc / allY
APureMito = APureMito / allY
AMix = AMix / allY

#avoid slopes between two points by repeating all points 10 times and spreading out the X
APGExp = rep(APureGlyc, each=10)
APMExp = rep(APureMito, each=10)
APBExp = rep(AMix, each=10)


df = tibble(x=log2(((1:1000) + 5)/1000 + 0.10) , Mito = APMExp, Glyc = APGExp, Both = APBExp, zeros = rep(0,1000))


cols3 = c('#c0D0F8','#C0E9B0','#FBF0D0')
names3 = c("Mit.", "Glyc.", "Both")

pZ = ggplot(df, aes(x=x)) +
  geom_area(aes(y=Mito+Glyc+Both, fill="Mix")) + 
  geom_area(aes(y=Mito+Glyc, fill="Glyc")) + 
  geom_area(aes(y=Mito, fill="Mito")) +
  geom_line(aes(y=Mito+Glyc+Both)) + 
  geom_line(aes(y=Mito+Glyc)) + 
  geom_line(aes(y=Mito)) +
  geom_line(aes(y=zeros)) +
  scale_fill_manual(values = cols3, labels = names3, breaks=c('Mito', 'Glyc', 'Mix')) +
  scale_x_continuous(breaks = log2(c(0,0.10,0.25,0.5,1) + 0.10), labels = c('0','0.10','0.25','0.5','1' )) +
  scale_y_continuous(breaks = c(0,1), labels = c('0','1' )) +
  #scale_fill_discrete(breaks=c('Mito', 'Glyc', 'Mix')) +
  ggplot2::labs(y="ATP Frac.", x="Static utilization") +
  ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=14, color="black"), axis.text = element_text(size=14, color="black"), legend.text = element_text(size=14, color="black"), legend.position = "bottom", legend.box = "vertical", plot.margin = margin(t=20)) +
  #guides(colour = guide_legend(title="Curve", direction="vertical", title.position = "top", title.hjust = 0.5, order=1), fill = guide_legend(title.position = "top", title.hjust = 0.5, order=2, title = "Optimal pathway"))
  guides(fill = guide_legend(order=1, direction="horizontal", title = element_blank())) +
  annotate(geom="text", x=-1.57, y=0.5, label="Astrocytes", color="black", size=5)  #fontface=2
pZ


ggsave(
  paste0(figPath, "CombModel2.eps"),
  plot = pZ,
  width = 3, height = 2.1, dpi = 300)


#################################
# Modeling of the full model - Transportation
#################################


d = readMat("data/FullModelOutputTransport.mat")$d[,,1] 


ATPProdGlycNeur = as.numeric(d$ATPProdGlycNeur)
ATPProdGlycAstr = as.numeric(d$ATPProdGlycAstr)
ATPProdMitoNeur = as.numeric(d$ATPProdMitoNeur)
ATPProdMitoAstr = as.numeric(d$ATPProdMitoAstr)

#Neurons
##############


NPureGlyc = ATPProdGlycNeur - ATPProdMitoNeur/14.5
NPureGlyc[NPureGlyc < 0] = 0
NPureMito = ATPProdMitoNeur - ATPProdGlycNeur*14.5
NPureMito[NPureMito < 0] = 0
NMix = ATPProdGlycNeur + ATPProdMitoNeur - NPureGlyc - NPureMito

allY = NPureGlyc + NPureMito + NMix
NPureGlyc = NPureGlyc / allY
NPureMito = NPureMito / allY
NMix = NMix / allY

#avoid slopes between two points by repeating all points 10 times and spreading out the X
NPGExp = rep(NPureGlyc, each=10)
NPMExp = rep(NPureMito, each=10)
NPBExp = rep(NMix, each=10)


df = tibble(x=log2(((1:1000) + 5)/1000 + 0.10) , Mito = NPMExp, Glyc = NPGExp, Both = NPBExp, zeros = rep(0,1000))


cols3 = c('#c0D0F8','#C0E9B0','#FBF0D0')

names3 = c("Mito.", "Glyc.", "Both")

pZ = ggplot(df, aes(x=x)) +
  geom_area(aes(y=Mito+Glyc+Both, fill="Mix")) + 
  geom_area(aes(y=Mito+Glyc, fill="Glyc")) + 
  geom_area(aes(y=Mito, fill="Mito")) +
  geom_line(aes(y=Mito+Glyc+Both)) + 
  geom_line(aes(y=Mito+Glyc)) + 
  geom_line(aes(y=Mito)) +
  geom_line(aes(y=zeros)) +
  scale_fill_manual(values = cols3, labels = names3, breaks=c('Mito', 'Glyc', 'Mix')) +
  scale_x_continuous(breaks = log2(c(0,0.10,0.25,0.5,1) + 0.10), labels = c('0','0.10','0.25','0.5','1' )) +
  scale_y_continuous(breaks = c(0,1), labels = c('0','1' )) +
  ggplot2::labs(y="ATP Frac.", x="Static utilization") +
  ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=14, color="black"), axis.text = element_text(size=14, color="black"), legend.text = element_text(size=14, color="black"), legend.position = "none", legend.box = "vertical", plot.margin = margin(t=20), axis.text.x = element_blank(), axis.title.x = element_blank()) +
  #guides(fill = guide_legend(order=1, direction="horizontal", title = element_blank()))
  annotate(geom="text", x=-1.54, y=0.5, label="Neurons", color="black", size=5)  #fontface=2
pZ


ggsave(
  paste0(figPath, "CombModel1Transp.eps"),
  plot = pZ,
  width = 3, height = 1.13, dpi = 300)


#Astrocytes
##############

APureGlyc = ATPProdGlycAstr - ATPProdMitoAstr/14.5
APureGlyc[APureGlyc < 0] = 0
APureMito = ATPProdMitoAstr - ATPProdGlycAstr*14.5
APureMito[APureMito < 0] = 0
AMix = ATPProdGlycAstr + ATPProdMitoAstr - APureGlyc - APureMito

allY = APureGlyc + APureMito + AMix
APureGlyc = APureGlyc / allY
APureMito = APureMito / allY
AMix = AMix / allY

#avoid slopes between two points by repeating all points 10 times and spreading out the X
APGExp = rep(APureGlyc, each=10)
APMExp = rep(APureMito, each=10)
APBExp = rep(AMix, each=10)


df = tibble(x=log2(((1:1000) + 5)/1000 + 0.10) , Mito = APMExp, Glyc = APGExp, Both = APBExp, zeros = rep(0,1000))


cols3 = c('#c0D0F8','#C0E9B0','#FBF0D0')
names3 = c("Mit.", "Glyc.", "Both")

pZ = ggplot(df, aes(x=x)) +
  geom_area(aes(y=Mito+Glyc+Both, fill="Mix")) + 
  geom_area(aes(y=Mito+Glyc, fill="Glyc")) + 
  geom_area(aes(y=Mito, fill="Mito")) +
  geom_line(aes(y=Mito+Glyc+Both)) + 
  geom_line(aes(y=Mito+Glyc)) + 
  geom_line(aes(y=Mito)) +
  geom_line(aes(y=zeros)) +
  scale_fill_manual(values = cols3, labels = names3, breaks=c('Mito', 'Glyc', 'Mix')) +
  scale_x_continuous(breaks = log2(c(0,0.10,0.25,0.5,1) + 0.10), labels = c('0','0.10','0.25','0.5','1' )) +
  scale_y_continuous(breaks = c(0,1), labels = c('0','1' )) +
  #scale_fill_discrete(breaks=c('Mito', 'Glyc', 'Mix')) +
  ggplot2::labs(y="ATP Frac.", x="Static utilization") +
  ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=14, color="black"), axis.text = element_text(size=14, color="black"), legend.text = element_text(size=14, color="black"), legend.position = "bottom", legend.box = "vertical", plot.margin = margin(t=20)) +
  #guides(colour = guide_legend(title="Curve", direction="vertical", title.position = "top", title.hjust = 0.5, order=1), fill = guide_legend(title.position = "top", title.hjust = 0.5, order=2, title = "Optimal pathway"))
  guides(fill = guide_legend(order=1, direction="horizontal", title = element_blank())) +
  annotate(geom="text", x=-1.57, y=0.5, label="Astrocytes", color="black", size=5)  #fontface=2
pZ


ggsave(
  paste0(figPath, "CombModel2Transp.eps"),
  plot = pZ,
  width = 3, height = 2.1, dpi = 300)




