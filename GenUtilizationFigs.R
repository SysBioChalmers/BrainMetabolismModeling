library(R.matlab)
library(tidyverse)
library(ggplot2)
library(ggpubr)

figPath = "Z:/projects/ANLS Modeling/Figures/"

setwd("C:/Work/MatlabCode/projects/BrainMetabolismModeling/BrainMetabolismModeling/")

#These were calculated with the model in MATLAB
mitoDivGlycWithMT = 13.9088
mitoDivGlycWithoutMT = 10.5520
glycProtPerATP = 0.00312



#################################
#4B
#################################
s = readMat("data/2B-DData.mat")$s[,,1] 
transp = as.numeric(s$transp)
EAMCAGlyc = as.numeric(s$b.EAMCAGlyc)
EAMCAMito = as.numeric(s$b.EAMCAMito)

names = c("Glyc.", "Mito. resp.");
allX = c(transp,transp) 
allY = c(EAMCAGlyc, EAMCAMito) 
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
  width = 3, height = 3.2, dpi = 300)


#################################
#4C
#################################
s = readMat("data/2B-DData.mat")$s[,,1] 
util = as.numeric(s$util)
EAMCAGlyc = as.numeric(s$c.EAMCAGlyc)
EAMCAMito = as.numeric(s$c.EAMCAMito)
EAMCAMitoMob = as.numeric(s$c.EAMCAMitoMob)

names = c("Glyc.", "Mit.", "Mit. mob.       " );
#allX = c(util,util,util[1:40]) 
allX = c(util,util,util) 
allY = c(EAMCAGlyc, EAMCAMito, EAMCAMitoMob) 
curve = factor(c(rep(1,length(util)), rep(2,length(util)), rep(3,length(EAMCAMitoMob))), 1:3, names)

ds = tibble(x=allX, y=log2(allY), curve = curve)

cols = c('#95B37B','#6B97BC','#6B97BC')
pY = ggplot(ds, aes(x = x, y = y, colour = curve, linetype = curve)) +
  geom_line(size=1.3) +
  scale_linetype_manual(values = c(1,1,2), labels = names) +
  scale_color_manual(values = c(cols[1],cols[2],cols[3]), labels = names) +
  ggplot2::labs(y=expression("Log"[2]*"(EAMCA)"), x="Static utilization", colour="Pathway", linetype="Pathway") +
  ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=14, color="black"), axis.text = element_text(size=14, color="black"), legend.text = element_text(size=14, color="black"), legend.position = "bottom", legend.title=element_blank())
pY


ggsave(
  paste0(figPath, "MaintPerUtil.png"),
  plot = pY,
  width = 3, height = 3.2, dpi = 300)

ggsave(
  paste0(figPath, "MaintPerUtil.eps"),
  plot = pY,
  width = 3, height = 3.2, dpi = 300)


#################################
#4D
#################################
s = readMat("data/2B-DData.mat")$s[,,1] 
util = as.numeric(s$util2)
EAMCAGlyc = as.numeric(s$d.EAMCAGlyc)
EAMCAMitoMob = as.numeric(s$d.EAMCAMitoMob)

#check at which positions glycolysis is more optimal. The lower switch can be calculated as 
#the saved maintenance cost per lactate output, where lactate output is just represented by

#check at which positions glycolysis is more optimal. The lower switch can be calculated as 
#the saved maintenance cost per lactate output, where lactate output is just represented by
#the static utilization, since the enzyme is inactive the rest of the time.
diff = EAMCAMitoMob - EAMCAGlyc
sel = diff > 0
gainY = diff[sel]

firstEdge = 0.0473 #0.0529
secEdge = 0.0791# 0.0889  # y = 2.155420e-03 at 0.4

allX = c(firstEdge,util[sel]) #firstEdge is for making the curve go to the bottom, the x value - manual interpolation
allY = c(min(gainY),gainY) #min(gainY) is to find the bottom of the curve
allY[secEdge*1000]
ds = tibble(x=log2(allX + 0.10), y=log2(allY))

#cols2 = c('#c0D0F8','#C0E9B0','#FBF0D0')
#names2 = c('Mit.', 'Glyc.','Both')
cols2 = c('#c0D0F8','#C0E9B0')
names2 = c('Mito. Resp.', 'Glyc.')

maxY = -1;

#d = tibble(x1=log2(c(0.01,firstEdge,secEdge,0.40) + 0.10), x2=log2(c(firstEdge, secEdge, 0.40, 1) + 0.10), y1=log2(rep(min(allY), 4)), y2=rep(maxY, 4), Optimal = factor(c(1,3,2,3),1:3,names2))
d = tibble(x1=log2(c(0.01,firstEdge) + 0.10), x2=log2(c(firstEdge, 1) + 0.10), y1=log2(rep(min(allY), 2)), y2=rep(maxY, 2), Optimal = factor(c(1,2),1:2,names2))
#dLine1 = tibble(x=log2(c(secEdge,secEdge) + 0.10), y=c(min(ds$y), maxY))
#dLine2 = tibble(x=log2(c(0.40,0.40) + 0.10), y=c(min(ds$y), maxY))
#dLine3 = tibble(x=log2(c(firstEdge,firstEdge) + 0.10), y=c(min(ds$y), maxY))

pZ = ggplot() +
  geom_rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=Optimal), color=NA, size=0,alpha=1) +
  geom_line(data = ds,size=1.3, mapping=aes(x = x, y = y), color = "black") +
#  geom_line(data=dLine1,aes(x = x, y=y), color = "black", linetype = "dashed") +
#  geom_line(data=dLine2,aes(x = x, y=y), color = "black", linetype = "dashed") +
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
  width = 3.05, height = 3.35, dpi = 300)


#################################
# The car transport analogy
#################################
day = 1:20
need = c(2,1,2,3,1,4,2,2,3,6,5,7,8,4,2,1,2,4,3,3)
ds = tibble(day=as.factor(day), need=need)

pY = ggplot(ds, aes(x = day, y = need)) +
  geom_bar(stat="identity", size=1.3) +
  labs(y="Cars needed", x="Day") +
  ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=14, color="black"), axis.text = element_text(size=14, color="black"), legend.text = element_text(size=14, color="black"), legend.position = "bottom", legend.title=element_blank())
pY

ggsave(
  paste0(figPath, "Car demand.png"),
  plot = pY,
  width = 5.9, height = 2.5, dpi = 300)

pY = ggplot(ds, aes(x = day, y = need)) +
  geom_bar(stat="identity", size=1.3) +
  labs(y="ATP demand", x="Hour") +
  ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=14, color="black"), axis.text = element_text(size=14, color="black"), legend.text = element_text(size=14, color="black"), legend.position = "bottom", legend.title=element_blank())
pY

ggsave(
  paste0(figPath, "ATP demand.png"),
  plot = pY,
  width = 5.9, height = 2.5, dpi = 300)



hist = tabulate(need)
cs = rev(cumsum(rev(hist)))
ds = tibble(cars=as.factor(1:8), util=cs/20)

pY = ggplot(ds, aes(x = cars, y = util)) +
  geom_bar(stat="identity", size=1.3) +
  labs(y="Utilization", x="Car slot") +
  ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=14, color="black"), axis.text = element_text(size=14, color="black"), legend.text = element_text(size=14, color="black"), legend.position = "bottom", legend.title=element_blank())
pY

ggsave(
  paste0(figPath, "Car util.png"),
  plot = pY,
  width = 5.5, height = 2.5, dpi = 300)

util= cs/20
util2 = util + (1- util)*0.7
ds = tibble(cars=rep(as.factor(1:8),2), pathway=factor(c(rep(1,8), rep(2,8)), c(1,2), c("Gly.", "Mit.")), util=c(util, util2))


bp = ggplot(data=dfPlot, aes(x=x, y=y, fill=pathway)) +
  geom_bar(stat="identity",position=position_dodge()) +
  geom_errorbar(aes(ymin=lb, ymax=ub), width=.2,
                position=position_dodge(.9)) +
  coord_cartesian(ylim=c(0.8, 0.89)) +
  labs( y="Correlation, 10x vs bulk", x="Covariates regressed out", fill="Fit: ") +
  scale_fill_manual(values=c("#82E182", "#229A22", "#E47060", "#AC210E")) + 
  theme_bw() +
  theme(axis.text = element_text(size = 8), legend.position="bottom")


pY = ggplot(ds, aes(x = cars, y = util, fill=pathway)) +
  geom_bar(stat="identity", size=1.3,position=position_dodge()) +
  labs(y="Utilization", x="ATP production capacity slice") +
  ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_fill_manual(values=c("#444444", "#AAAAAA")) + 
  theme(text = element_text(size=14, color="black"), axis.text = element_text(size=14, color="black"), 
        legend.text = element_text(size=14, color="black"), legend.position = c(0.8,0.85), legend.title=element_blank(), 
        legend.direction="horizontal")
pY

ggsave(
  paste0(figPath, "Cat cap util.png"),
  plot = pY,
  width = 5.5, height = 2.5, dpi = 300)

#################################
# Modeling of the full model - Utilization (Fig. 4G)
#################################


library(R.matlab)

d = readMat("data/FullModelOutputRedMitMob.mat")$d[,,1] 


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
  ggplot2::labs(y="ATP Frac.", x="Static utilization") +
  ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=14, color="black"), axis.text = element_text(size=14, color="black"), legend.text = element_text(size=14, color="black"), legend.position = "bottom", legend.box = "vertical", plot.margin = margin(t=20)) +
  guides(fill = guide_legend(order=1, direction="horizontal", title = element_blank())) +
  annotate(geom="text", x=-1.57, y=0.5, label="Astrocytes", color="black", size=5)  #fontface=2
pZ


ggsave(
  paste0(figPath, "CombModel2.eps"),
  plot = pZ,
  width = 3, height = 2.1, dpi = 300)


#################################
# Modeling of the full model - Transportation (Fig. S3)
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
  ggplot2::labs(y="ATP Frac.", x="Static utilization") +
  ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=14, color="black"), axis.text = element_text(size=14, color="black"), legend.text = element_text(size=14, color="black"), legend.position = "bottom", legend.box = "vertical", plot.margin = margin(t=20)) +
  guides(fill = guide_legend(order=1, direction="horizontal", title = element_blank())) +
  annotate(geom="text", x=-1.57, y=0.5, label="Astrocytes", color="black", size=5)  #fontface=2
pZ


ggsave(
  paste0(figPath, "CombModel2Transp.eps"),
  plot = pZ,
  width = 3, height = 2.1, dpi = 300)



#################################
# Modeling of the full model - Check importance of transport cost of MT genes (Fig. S3)
#################################


library(R.matlab)

d = readMat("data/FullModelOutputRedMitMobTMT.mat")$d[,,1] 


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
  annotate(geom="text", x=-1.54, y=0.5, label="Neurons", color="black", size=5)  #fontface=2
pZ


ggsave(
  paste0(figPath, "CombModelTMT1.eps"),
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
  ggplot2::labs(y="ATP Frac.", x="Static utilization") +
  ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=14, color="black"), axis.text = element_text(size=14, color="black"), legend.text = element_text(size=14, color="black"), legend.position = "bottom", legend.box = "vertical", plot.margin = margin(t=20)) +
  guides(fill = guide_legend(order=1, direction="horizontal", title = element_blank())) +
  annotate(geom="text", x=-1.57, y=0.5, label="Astrocytes", color="black", size=5)  #fontface=2
pZ


ggsave(
  paste0(figPath, "CombModelTMT2.eps"),
  plot = pZ,
  width = 3, height = 2.1, dpi = 300)



#############################
# Fig. 4E
#############################

d1 = readMat("data/simpleModelOutputRedMitMob.mat")$d[,,1] 
d2 = readMat("data/simpleModelOutputTransport.mat")$d[,,1] 

glycATPGenNeur1 = as.numeric(d1$glycATPGenNeur)
glycATPGenAstr1 = as.numeric(d1$glycATPGenAstr)
mitoATPGenNeur1 = as.numeric(d1$mitoATPGenNeur)
mitoATPGenAstr1 = as.numeric(d1$mitoATPGenAstr)

glycATPGenNeur2 = as.numeric(d2$glycATPGenNeur)
glycATPGenAstr2 = as.numeric(d2$glycATPGenAstr)
mitoATPGenNeur2 = as.numeric(d2$mitoATPGenNeur)
mitoATPGenAstr2 = as.numeric(d2$mitoATPGenAstr)

sumNeur = glycATPGenNeur1 + mitoATPGenNeur1
sumAstr = glycATPGenAstr1 + mitoATPGenAstr1

y = c(glycATPGenAstr1/sumAstr, mitoATPGenAstr1/sumAstr, glycATPGenNeur1/sumNeur, mitoATPGenNeur1/sumNeur)
x = factor(1:4,1:4,c("Astr. Glyc.","Astr. Mit.", "Neur. Glyc.","Neur. Mit."))
ct = as.factor(c(1,1,2,2))
df = tibble(x=x,y=y,ct=ct)

p1<-ggplot(data=df, aes(x=x, y=y, fill=ct)) +
  #  geom_bar(stat="identity", width = 0.8) +
  geom_bar(stat="identity", color = "NA") +
  scale_fill_manual(values=c("#AAAAAA", "#444444")) +
  coord_flip() +
  labs(x="",y="ATP frac.") +
  scale_y_continuous(breaks = c(0,1), labels = c('0','1' )) +
  ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.ticks.y = element_blank()) +
  theme(text = element_text(size=14, color="black"), axis.text.x = element_text(size=14, color="black"), 
        axis.text.y = element_text(size=14, color="black", vjust=0.3), legend.text = element_text(size=14, color="black"), 
        panel.background = element_blank(), plot.margin = margin(r = 13), legend.position = "none")
p1

sumNeur = glycATPGenNeur2 + mitoATPGenNeur2
sumAstr = glycATPGenAstr2 + mitoATPGenAstr2

y = c(glycATPGenAstr2/sumAstr, mitoATPGenAstr2/sumAstr, glycATPGenNeur2/sumNeur, mitoATPGenNeur2/sumNeur)
x = factor(1:4,1:4,c("Astr. glyc.","Astr. mit.", "Neur. glyc.","Neur. mit."))
ct = as.factor(c(1,1,2,2))
df = tibble(x=x,y=y,ct=ct)

p2<-ggplot(data=df, aes(x=x, y=y, fill=ct)) +
  geom_bar(stat="identity", color = "NA") +
  scale_fill_manual(values=c("#AAAAAA", "#444444")) +
  coord_flip() +
  labs(x="",y="ATP frac.") +
  scale_y_continuous(breaks = c(0,1), labels = c('0','1' )) +
  ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.ticks.y = element_blank()) +
  theme(text = element_text(size=14, color="black"), axis.text.x = element_text(size=14, color="black"), 
        axis.text.y = element_text(size=14, color="black", vjust=0.3), legend.text = element_text(size=14, color="black"), 
        panel.background = element_blank(), plot.margin = margin(r = 13), legend.position = "none")
p2


ggsave(
  paste0(figPath, "Fig S2B.svg"),
  plot = p1,
  width = 3, height = 1.5, dpi = 300)

ggsave(
  paste0(figPath, "Fig S2C.svg"),
  plot = p2,
  width = 3, height = 1.5, dpi = 300)




################################
################################
# Below are tests for verifying
# Fig. 4B-D
################################
################################



############################
# Make a figure that shows how much more expensive it is to use
# mitochondria compared to glycolysis if they are static.
#################################

#values from the matlab modeling
#Just copied here from above for convenience
#mitoDivGlycWithMT = 13.9088
#mitoDivGlycWithoutMT = 10.5520


#scales are 1 for glyc, 13.9088 and 10.5520 for mito.
getEAMCA = function (scaleWithMT, scaleWithoutMT, trCost, util, glycProtPerATP) {
  totalCost = ((scaleWithMT + trCost*scaleWithoutMT)/util  - scaleWithMT)*glycProtPerATP
}

util = (1:100)*0.01
extraCostGlyc = rep(NA, length(util))
extraCostMito = rep(NA, length(util))
#extraCostMitoLim = rep(NA, 40)
extraCostMitoLim = rep(NA, length(util))
for (i in 1:length(util)) {
  #hardcode transport cost to 0%
  extraCostGlyc[i] = getEAMCA(1, 1, 0, util[i], glycProtPerATP)
  extraCostMito[i] = getEAMCA(mitoDivGlycWithMT, mitoDivGlycWithoutMT, 0, util[i], glycProtPerATP)
}

for (i in 1:length(extraCostMitoLim)) {
  #hardcode transport cost to 0%
  #extraCostMitoLim[i] = getEAMCA(mitoDivGlycWithMT, mitoDivGlycWithoutMT, 0, 0.4, glycProtPerATP)
  extraCostMitoLim[i] = getEAMCA(mitoDivGlycWithMT, mitoDivGlycWithoutMT, 0, util[i] + 0.4*(1-util[i]), glycProtPerATP)
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

###################################
# Similar plot for transport cost
###################################

util2 = 1
transp = (1:200)*0.01
extraCostGlyc = rep(NA, length(transp))
extraCostMito = rep(NA, length(transp))
for (i in 1:length(transp)) {
  extraCostGlyc[i] = getEAMCA(1,1, transp[i], util2, glycProtPerATP)
  extraCostMito[i] = getEAMCA(mitoDivGlycWithMT, mitoDivGlycWithoutMT, transp[i], util2, glycProtPerATP)
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


###########################
#Now a  plot showing optimal usage (glyc)
#neurons and astrocytes when lactate export is limited
###########################

#calculate the cost for both glyc and mito if we assume mito mobility = 0.4 and transport = 0.1

util = (1:1000)*0.001
mitoExtraCost = rep(NA, length(util))
cytoExtraCost = rep(NA, length(util))

for (i in 1:length(util)) {
  #Transportation cost set to 0.1
  mitoExtraCost[i] = getEAMCA(mitoDivGlycWithMT, mitoDivGlycWithoutMT, 0.1, util[i] + 0.4*(1-util[i]), glycProtPerATP)
  cytoExtraCost[i] = getEAMCA(1,1, 0.1, util[i], glycProtPerATP)
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






