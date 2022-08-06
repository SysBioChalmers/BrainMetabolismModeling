library(R.matlab)
library(tidyverse)
library(ggplot2)
library(ggpubr)

figPath = "Z:/projects/ANLS Modeling/Figures/"

setwd("C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Single-cell Modeling/NeuronModeling/")

#These were calculated with the model in MATLAB
mitoDivGlycWithMT = 13.9088
mitoDivGlycWithoutMT = 10.5520

########################
# Make the range of transport cost from 0.01 to 2
# trCost = Abs extra transport const / abs cost prot maintenance in "normal" cells
# Balance: utMito = ((1 + trCost)*mitoDivGlyc) /((1 + trCost)/Utc - 1 + mitoDivGlyc)
# for different curves
########################

getUtMito = function (mitoDivGlycWithMT, mitoDivGlycWithoutMT, trCost, utCyto) {
  utMito = (mitoDivGlycWithMT + trCost*mitoDivGlycWithoutMT) /((1 + trCost)/utCyto - 1 + mitoDivGlycWithMT)
}



trCost = seq(0.01, 2, by=0.01)
ut01 = rep(NA,length(trCost))
for (i in 1:length(trCost)) {
  ut01[i] = getUtMito(mitoDivGlycWithMT, mitoDivGlycWithoutMT, trCost[i], 0.01)
}

ut05 = rep(NA,length(trCost))
for (i in 1:length(trCost)) {
  ut05[i] = getUtMito(mitoDivGlycWithMT, mitoDivGlycWithoutMT, trCost[i], 0.05)
}

ut10 = rep(NA,length(trCost))
for (i in 1:length(trCost)) {
  ut10[i] = getUtMito(mitoDivGlycWithMT, mitoDivGlycWithoutMT, trCost[i], 0.1)
}

names = c("0.01", "0.05","0.10");
allX = c(trCost,trCost,trCost) 
allY = c(ut01, ut05, ut10) 
curve = factor(c(rep(1,length(trCost)), rep(2,length(trCost)), rep(3,length(trCost))), 
               1:3, names)

ds = tibble(x=allX, y=allY, curve = curve)

cols = c('#95B37B','#E0B064','#6B97BC')

pX = ggplot(ds, aes(x = x, y = y, colour = curve)) +
  geom_line(size=1.3) +
#  scale_linetype_manual(values = c(1,2), labels = typeNames) +
  scale_color_manual(values = c(cols[1],cols[2],cols[3]), labels = names) +
  ggplot2::labs(y="Req. mito. utilization", x="Transportation cost", colour=expression("Static utilization (U"[s]*")")) +
  ylim(0,1) +
  ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=14, color="black"), axis.text = element_text(size=14, color="black"), legend.text = element_text(size=14, color="black"), legend.position = "bottom") +
  guides(colour = guide_legend(title.position = "top", title.hjust = 0.45))
pX


ggsave(
  paste0(figPath, "FigReqMitoUtil.png"),
  plot = pX,
  width = 3, height = 3.5, dpi = 300)

ggsave(
  paste0(figPath, "FigReqMitoUtil.eps"),
  plot = pX,
  width = 3, height = 3.5, dpi = 300)

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
extraCostMitoLim = rep(NA, 40)
for (i in 1:length(util)) {
  #hardcode transport cost to 0%
  extraCostGlyc[i] = getEAMCA(1,1, 0, util[i])
  extraCostMito[i] = getEAMCA(mitoDivGlycWithMT, mitoDivGlycWithoutMT, 0, util[i])
}

for (i in 1:length(extraCostMitoLim)) {
  #hardcode transport cost to 0.1%
  extraCostMitoLim[i] = getEAMCA(mitoDivGlycWithMT, mitoDivGlycWithoutMT, 0.1, 0.4)
}

names = c("Glyc.", "Mit.", "Mit. mob.       " );
allX = c(util,util,util[1:40]) 
allY = c(extraCostGlyc, extraCostMito, extraCostMitoLim) 
curve = factor(c(rep(1,length(util)), rep(2,length(util)), rep(3,40)), 1:3, names)

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

util2 = 1#0.5
transp = (1:200)*0.01
extraCostGlyc = rep(NA, length(transp))
extraCostMito = rep(NA, length(transp))
for (i in 1:length(transp)) {
  #hardcode transport cost to 10%
  extraCostGlyc[i] = getEAMCA(1,1, transp[i], util2)
  extraCostMito[i] = getEAMCA(mitoDivGlycWithMT, mitoDivGlycWithoutMT, transp[i], util2)
}

names = c("Glyc.", "Mito. resp.");
allX = c(transp,transp) 
allY = c(extraCostGlyc, extraCostMito) 
curve = factor(c(rep(1,length(transp)), rep(2,length(transp))), 1:2, names)

ds = tibble(x=allX, y=log2(allY), curve = curve)

#cols = c('#95B37B','#E0B064','#6B97BC')
cols = c('#95B37B','#6B97BC')
#y=expression(Log[2]*" fold change maint. cost")
pY = ggplot(ds, aes(x = x, y = y, colour = curve)) +
  geom_line(size=1.3) +
  #  scale_linetype_manual(values = c(1,2), labels = typeNames) +
  scale_color_manual(values = c(cols[1],cols[2],cols[3]), labels = names) +
  ggplot2::labs(y=expression("Log"[2]*"(EAMCA)"), x="Transportation cost", colour="Pathway") +
  #  ylim(0,1) +
  ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=14, color="black"), axis.text = element_text(size=14, color="black"), legend.text = element_text(size=14, color="black"), legend.position = "bottom", legend.title = element_blank())
#  guides(colour = guide_legend(title.position = "top", title.hjust = 0.5))
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

#getExtraCost = function (scaleWithMT, scaleWithoutMT, trCost, util) {
#  extraCost = (scaleWithMT + trCost*scaleWithoutMT)/util - scaleWithMT
#}

range = (1:100)*0.01
mitoExtraCost = rep(NA, length(range))
cytoExtraCost = rep(NA, length(range))

for (i in 1:length(range)) {
  mitoExtraCost[i] = getEAMCA(mitoDivGlycWithMT, mitoDivGlycWithoutMT, 0.1, max(range[i], 0.4))
  cytoExtraCost[i] = getEAMCA(1,1, 0.1, range[i])
}

#check at which positions glycolysis is more optimal. The lower switch can be calculated as 
#the saved maintenance cost per lactate output, where lactate output is just represented by

#plot 2:
#check at which positions glycolysis is more optimal. The lower switch can be calculated as 
#the saved maintenance cost per lactate output, where lactate output is just represented by
#the static utilization, since the enzyme is inactive the rest of the time.
gainPerLactateOutput = (mitoExtraCost - cytoExtraCost)/range
gainYTmp = gainPerLactateOutput[gainPerLactateOutput > 0]
#just add one value at the bottom of the fig to make the curve go all the way down
gainY = c(min(gainYTmp),gainYTmp)



#allX = c(0.028,range[3:100]) 
#allY = c(gainY) 
names = c("Glycolysis", "Mitochondrial resp.", "Gain per lactate exp.");
allX = c(0.045,range[5:100]) #0.045 is for making the curve go to the bottom, the x value
allY = gainY 

ds = tibble(x=log2(allX + 0.10), y=log2(allY))

#cols2 = c('#D5E5F8','#D8E9CB','#FBF0D0')
cols2 = c('#c0D0F8','#C0E9B0','#FBF0D0')
names2 = c('Mit.', 'Glyc.','Both')
#d = tibble(x1=log2(c(0.01,0.028,0.036,0.40) + 0.10), x2=log2(c(0.028, 0.036, 0.40, 1) + 0.10), y1=log2(rep(min(allY), 4)), y2=rep(8, 4), Optimal = factor(c(1,3,2,3),1:3,names2))
#dLine1 = tibble(x=log2(c(0.036,0.036) + 0.10), y=c(min(ds$y), 8))
#dLine2 = tibble(x=log2(c(0.40,0.40) + 0.10), y=c(min(ds$y), 8))
#dLine3 = tibble(x=log2(c(0.028,0.028) + 0.10), y=c(min(ds$y), 8))

d = tibble(x1=log2(c(0.01,0.044,0.050,0.40) + 0.10), x2=log2(c(0.044, 0.050, 0.40, 1) + 0.10), y1=log2(rep(min(allY), 4)), y2=rep(8, 4), Optimal = factor(c(1,3,2,3),1:3,names2))
dLine1 = tibble(x=log2(c(0.050,0.050) + 0.10), y=c(min(ds$y), 8))
dLine2 = tibble(x=log2(c(0.40,0.40) + 0.10), y=c(min(ds$y), 8))
dLine3 = tibble(x=log2(c(0.044,0.044) + 0.10), y=c(min(ds$y), 8))

pZ = ggplot() +
  geom_rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=Optimal), color=NA, size=0,alpha=1) +
  geom_line(data = ds,size=1.3, mapping=aes(x = x, y = y), color = "black") +
  geom_line(data=dLine1,aes(x = x, y=y), color = "black", linetype = "dashed") +
  geom_line(data=dLine2,aes(x = x, y=y), color = "black", linetype = "dashed") +
  geom_line(data=dLine3,aes(x = x, y=y), color = "black", linetype = "dashed") +
  scale_fill_manual(values = c(cols2[1],cols2[2],cols2[3]), labels = names2) +
  scale_x_continuous(breaks = log2(c(0,0.10,0.25,0.5,1) + 0.10), labels = c('0','0.10','0.25','0.5','1' )) +
  ggplot2::labs(y=expression("Log"[2]*"(EAMCA gain per lact.)        "), x="Static utilization", colour="Pathway") +
  ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=14, color="black"), axis.text = element_text(size=14, color="black"), legend.text = element_text(size=14, color="black"), legend.position = "bottom", legend.box = "horizontal", plot.margin = margin(r = 13), legend.title.align = 0.5) +
  #guides(colour = guide_legend(title="Curve", direction="vertical", title.position = "top", title.hjust = 0.5, order=1), fill = guide_legend(title.position = "top", title.hjust = 0.5, order=2, title = "Optimal pathway"))
  guides(fill = guide_legend(title.position = "top", order=1, direction="horizontal", title = "Optimal pathway:")) + #, title.justification='center')) + 
  ylim(log2(min(allY)),8)
pZ

#axis.title.y = element_text(hjust=-0.5)
ggsave(
  paste0(figPath, "OptimalUse2B.png"),
  plot = pZ,
  width = 5.25, height = 2.54, dpi = 300)

ggsave(
  paste0(figPath, "OptimalUse2B.eps"),
  plot = pZ,
  width = 2.8, height = 3.35, dpi = 300)












#names = c("Glyc.", "Mito. resp.");
names = c("Glycolysis", "Mitochondrial resp.", "Gain per lactate exp.");
allX = c(range,range,c(0.045,range[5:100])) #0.045 is for making the curve go to the bottom, the x value
allY = c(cytoExtraCost, mitoExtraCost, gainY) 
curve = factor(c(rep(1,length(range)), rep(2,length(range)), rep(3,length(gainY))), 1:3, names)

ds = tibble(x=log2(allX + 0.10), y=log2(allY), curve = curve)

#cols = c('#95B37B','#E0B064','#6B97BC')
cols = c('#95B37B','#6B97BC','#E07064')
#cols2 = c('#C8D9BB','#B5C7DC','#ECD8B0')
#cols2 = c('#D8E9CB','#C5D7EC','#F6E6C0')
#cols2 = c('#E7F6D9','#D5E5F8','#FBF0D0')
cols2 = c('#D5E5F8','#D8E9CB','#FBF0D0')
#y=expression(Log[2]*" fold change maint. cost")
names2 = c('Mit.', 'Glyc.','Both')
d = tibble(x1=log2(c(0.01,0.036,0.40) + 0.10), x2=log2(c(0.036, 0.40, 1) + 0.10), y1=log2(rep(min(allY), 3)), y2=log2(rep(max(allY), 3)), Optimal = factor(1:3,1:3,names2))
dLine1 = tibble(x=log2(c(0.036,0.036) + 0.10), y=c(min(ds$y), max(ds$y)))
dLine2 = tibble(x=log2(c(0.40,0.40) + 0.10), y=c(min(ds$y), max(ds$y)))


pZ = ggplot() +
  geom_rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=Optimal), color=NA, size=0,alpha=1) +
  geom_line(data = ds,size=1.3, mapping=aes(x = x, y = y, colour = curve, linetype=curve)) +
  geom_line(data=dLine1,aes(x = x, y=y), color = "black", linetype = "dashed") +
  geom_line(data=dLine2,aes(x = x, y=y), color = "black", linetype = "dashed") +
  scale_color_manual(values = c(cols[1],cols[2],cols[3]), labels = names) +
  scale_linetype_manual(values = c('solid','solid','solid'), labels = names) +
  scale_fill_manual(values = c(cols2[1],cols2[2],cols2[3]), labels = names2) +
  scale_x_continuous(breaks = log2(c(0,0.10,0.25,0.5,1) + 0.10), labels = c('0','0.10','0.25','0.5','1' )) +
  ggplot2::labs(y=expression("Log"[2]*"(ext. cost per ATP)"), x="Static utilization", colour="Pathway") +
  ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=14, color="black"), axis.text = element_text(size=14, color="black"), legend.text = element_text(size=14, color="black"), legend.position = "bottom", legend.box = "vertical", plot.margin = margin(r = 13)) +
  #guides(colour = guide_legend(title="Curve", direction="vertical", title.position = "top", title.hjust = 0.5, order=1), fill = guide_legend(title.position = "top", title.hjust = 0.5, order=2, title = "Optimal pathway"))
  guides(colour = guide_legend(title=element_blank(), direction="vertical", order=1), 
         linetype = guide_legend(title=element_blank(), direction="vertical", order=1), 
         fill = guide_legend(title.position = "top", title.hjust = 0.5, order=2, title = "Optimal pathway"))
pZ


ggsave(
  paste0(figPath, "OptimalUse.png"),
  plot = pZ,
  width = 3, height = 4.3, dpi = 300)

ggsave(
  paste0(figPath, "OptimalUse.eps"),
  plot = pZ,
  width = 3, height = 4.54, dpi = 300)

#split into two plots instead:

#plot 1:

names = c("Glyc.", "Mito. resp.");
allX = c(range,range) 
allY = c(cytoExtraCost, mitoExtraCost) 
curve = factor(c(rep(1,length(range)), rep(2,length(range))),1:2, names)

ds = tibble(x=log2(allX + 0.10), y=log2(allY), curve = curve)

cols = c('#95B37B','#6B97BC')

pZ = ggplot() +
  geom_line(data = ds,size=1.3, mapping=aes(x = x, y = y, colour = curve, linetype=curve)) +
  scale_color_manual(values = c(cols[1],cols[2]), labels = names) +
  scale_linetype_manual(values = c('solid','solid'), labels = names) +
  scale_x_continuous(breaks = log2(c(0,0.10,0.25,0.5,1) + 0.10), labels = c('0','0.10','0.25','0.5','1' )) +
  ggplot2::labs(y=expression("Log"[2]*"(ext. cost per ATP)"), x="Static utilization", colour="Pathway") +
  ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=14, color="black"), axis.text = element_text(size=14, color="black"), legend.text = element_text(size=14, color="black"), legend.position = "bottom", legend.box = "vertical", plot.margin = margin(r = 13)) +
  #guides(colour = guide_legend(title="Curve", direction="vertical", title.position = "top", title.hjust = 0.5, order=1), fill = guide_legend(title.position = "top", title.hjust = 0.5, order=2, title = "Optimal pathway"))
  guides(colour = guide_legend(title=element_blank(), direction="horizontal", order=1), 
         linetype = guide_legend(title=element_blank(), direction="horizontal", order=1))
pZ

ggsave(
  paste0(figPath, "OptimalUse2a.png"),
  plot = pZ,
  width = 3, height = 4.3, dpi = 300)

ggsave(
  paste0(figPath, "OptimalUse2a.eps"),
  plot = pZ,
  width = 3, height = 4.54, dpi = 300)


#plot 2:
#check at which positions glycolysis is more optimal. The lower switch can be calculated as 
#the saved maintenance cost per lactate output, where lactate output is just represented by
#the static utilization, since the enzyme is inactive the rest of the time.
gainPerLactateOutput = (mitoExtraCost - cytoExtraCost)/range
gainYTmp = gainPerLactateOutput[gainPerLactateOutput > 0]
#just add one value at the bottom of the fig to make the curve go all the way down
gainY = c(min(c(gainYTmp, cytoExtraCost,mitoExtraCost)),gainYTmp)



allX = c(0.028,range[3:100]) 
allY = c(gainY) 

ds = tibble(x=log2(allX + 0.10), y=log2(allY))

#cols2 = c('#D5E5F8','#D8E9CB','#FBF0D0')
cols2 = c('#c0D0F8','#C0E9B0','#FBF0D0')
names2 = c('Mitochondrial. resp.', 'Glycolysis','Both')
d = tibble(x1=log2(c(0.01,0.028,0.036,0.40) + 0.10), x2=log2(c(0.028, 0.036, 0.40, 1) + 0.10), y1=log2(rep(min(allY), 4)), y2=log2(rep(max(allY), 4)), Optimal = factor(c(1,3,2,3),1:3,names2))
dLine1 = tibble(x=log2(c(0.036,0.036) + 0.10), y=c(min(ds$y), max(ds$y)))
dLine2 = tibble(x=log2(c(0.40,0.40) + 0.10), y=c(min(ds$y), max(ds$y)))
dLine3 = tibble(x=log2(c(0.028,0.028) + 0.10), y=c(min(ds$y), max(ds$y)))


pZ = ggplot() +
  geom_rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=Optimal), color=NA, size=0,alpha=1) +
  geom_line(data = ds,size=1.3, mapping=aes(x = x, y = y), color = "black") +
  geom_line(data=dLine1,aes(x = x, y=y), color = "black", linetype = "dashed") +
  geom_line(data=dLine2,aes(x = x, y=y), color = "black", linetype = "dashed") +
  geom_line(data=dLine3,aes(x = x, y=y), color = "black", linetype = "dashed") +
  scale_fill_manual(values = c(cols2[1],cols2[2],cols2[3]), labels = names2) +
  scale_x_continuous(breaks = log2(c(0,0.10,0.25,0.5,1) + 0.10), labels = c('0','0.10','0.25','0.5','1' )) +
  ggplot2::labs(y=expression("Log"[2]*"(ext. cost per ATP)"), x="Static utilization", colour="Pathway") +
  ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=14, color="black"), axis.text = element_text(size=14, color="black"), legend.text = element_text(size=14, color="black"), legend.position = "right", legend.box = "vertical", plot.margin = margin(r = 13)) +
  #guides(colour = guide_legend(title="Curve", direction="vertical", title.position = "top", title.hjust = 0.5, order=1), fill = guide_legend(title.position = "top", title.hjust = 0.5, order=2, title = "Optimal pathway"))
  guides(fill = guide_legend(title.position = "top", order=1, direction="vertical", title = "Optimal pathway:"))
pZ


ggsave(
  paste0(figPath, "OptimalUse2B.png"),
  plot = pZ,
  width = 5.25, height = 2.54, dpi = 300)

ggsave(
  paste0(figPath, "OptimalUse2B.eps"),
  plot = pZ,
  width = 5.25, height = 2.54, dpi = 300)


###########################
#Now a  plot showing optimal usage (glyc)
#neurons and astrocytes when lactate export is limited
#This is a different vairant with the legends placed differently
###########################

#calculate the cost for both glyc and mito if we assume mito util >= 0.3 and transport = 0.1

getExtraCost = function (scaleWithMT, scaleWithoutMT, trCost, util) {
  extraCost = (scaleWithMT + trCost*scaleWithoutMT)/util - scaleWithMT
}

range = (1:100)*0.01
mitoExtraCost = rep(NA, length(range))
cytoExtraCost = rep(NA, length(range))

for (i in 1:length(range)) {
  mitoExtraCost[i] = getExtraCost(mitoDivGlycWithMT, mitoDivGlycWithoutMT, 0.1, max(range[i], 0.3))
  cytoExtraCost[i] = getExtraCost(1,1, 0.1, range[i])
}
#plot(1:100,cytoExtraCost)
#plot(1:100,mitoExtraCost)

#check at which positions glycolysis is more optimal. The lower switch can be calculated as 
#the saved maintenance cost per lactate output, where lactate output is just represented by
#the static utilization, since the enzyme is inactive the rest of the time.
gainPerLactateOutput = (mitoExtraCost - cytoExtraCost)/range/5
gainYTmp = gainPerLactateOutput[gainPerLactateOutput > 0]
#just add one value at the bottom of the fig to make the curve go all the way down
gainY = c(min(c(gainYTmp, cytoExtraCost,mitoExtraCost)),gainYTmp)



#names = c("Glyc.", "Mito. resp.");
names = c("Glycolysis", "Mitochondrial resp.", "Gain per lactate exp.");
allX = c(range,range,c(0.028,range[3:100])) 
allY = c(cytoExtraCost, mitoExtraCost, gainY) 
curve = factor(c(rep(1,length(range)), rep(2,length(range)), rep(3,length(gainY))), 1:3, names)

ds = tibble(x=log2(allX + 0.10), y=log2(allY), curve = curve)

#cols = c('#95B37B','#E0B064','#6B97BC')
cols = c('#95B37B','#6B97BC','#E07064')
#cols2 = c('#C8D9BB','#B5C7DC','#ECD8B0')
#cols2 = c('#D8E9CB','#C5D7EC','#F6E6C0')
#cols2 = c('#E7F6D9','#D5E5F8','#FBF0D0')
cols2 = c('#D5E5F8','#D8E9CB','#FBF0D0')
#y=expression(Log[2]*" fold change maint. cost")
names2 = c('Mitochondrial resp.', 'Glycolysis','Both')
d = tibble(x1=log2(c(0.01,0.036,0.40) + 0.10), x2=log2(c(0.036, 0.40, 1) + 0.10), y1=log2(rep(min(allY), 3)), y2=log2(rep(max(allY), 3)), Optimal = factor(1:3,1:3,names2))
dLine1 = tibble(x=log2(c(0.036,0.036) + 0.10), y=c(min(ds$y), max(ds$y)))
dLine2 = tibble(x=log2(c(0.40,0.40) + 0.10), y=c(min(ds$y), max(ds$y)))


pZ = ggplot() +
  geom_rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=Optimal), color=NA, size=0,alpha=1) +
  geom_line(data = ds,size=1.3, mapping=aes(x = x, y = y, colour = curve, linetype=curve)) +
  geom_line(data=dLine1,aes(x = x, y=y), color = "black", linetype = "dashed") +
  geom_line(data=dLine2,aes(x = x, y=y), color = "black", linetype = "dashed") +
  scale_color_manual(values = c(cols[1],cols[2],cols[3]), labels = names) +
  scale_linetype_manual(values = c('solid','solid','solid'), labels = names) +
  scale_fill_manual(values = c(cols2[1],cols2[2],cols2[3]), labels = names2) +
  scale_x_continuous(breaks = log2(c(0,0.10,0.25,0.5,1) + 0.10), labels = c('0','0.10','0.25','0.5','1' )) +
  ggplot2::labs(y=expression("Log"[2]*"(ext. cost per ATP)"), x="Static utilization", colour="Pathway") +
  ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(text = element_text(size=14, color="black"), axis.text = element_text(size=14, color="black"), legend.text = element_text(size=14, color="black"), legend.position = "right", legend.box = "vertical", plot.margin = margin(t=20)) +
  #guides(colour = guide_legend(title="Curve", direction="vertical", title.position = "top", title.hjust = 0.5, order=1), fill = guide_legend(title.position = "top", title.hjust = 0.5, order=2, title = "Optimal pathway"))
  guides(colour = guide_legend(title="Extra maint. cost per ATP:", direction="vertical", order=1), 
         linetype = guide_legend(title="Extra maint. cost per ATP:", direction="vertical", order=1), 
         fill = guide_legend(title.position = "top", order=2, direction="vertical", title = "Optimal pathway:"))
pZ


ggsave(
  paste0(figPath, "OptimalUse2.png"),
  plot = pZ,
  width = 6, height = 3, dpi = 300)

ggsave(
  paste0(figPath, "OptimalUse2.eps"),
  plot = pZ,
  width = 5.47, height = 2.78, dpi = 300)











###########################
#Now a bar plot for lactate export between 
#neurons and astrocytes when lactate export is limited
###########################

#these numbers were calculated in Matlab, RunFullModelSim.m
expNeur = 0.5380
expAstr = 1

y = c(expNeur,expAstr) * 100
x = factor(1:2,1:2,c("Neurons","Astrocytes"))
df = tibble(x=x,y=y)

p<-ggplot(data=df, aes(x=x, y=y)) +
#  geom_bar(stat="identity", width = 0.8) +
  geom_bar(stat="identity") +
  coord_flip() +
  labs(x="",y="Lactate exp. (%)") +
  #theme_minimal() + 
  ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.ticks.y = element_blank()) +
  theme(text = element_text(size=14, color="black"), axis.text.x = element_text(size=14, color="black"), axis.text.y = element_text(size=14, color="black", vjust=0.3), legend.text = element_text(size=14, color="black"), panel.background = element_blank(), plot.margin = margin(r = 13))
p

ggsave(
  paste0(figPath, "AstrNeurBar.png"),
  plot = p,
  width = 2.7, height = 1.1, dpi = 300)

ggsave(
  paste0(figPath, "AstrNeurBar.svg"),
  plot = p,
  width = 2.7, height = 1.1, dpi = 300)


#names = c("Glyc., Glyc = 0.1", "Mito., Glyc = 0.1","Glyc., Glyc = 0.3", "Mito., Glyc = 0.3", "Glyc., Glyc = 0.5", "Mito., Glyc = 0.5");
#allX = c(0.01,2,trCost, 0.01,2,trCost, 0.01,2,trCost) 
#allY = c(0.01,0.01,ut01, 0.05,0.05,ut05, 0.1,0.1,ut10) 
#curve = factor(c(1,1,rep(2,length(trCost)), 3,3,rep(4,length(trCost)), 5,5,rep(6,length(trCost))), 
#               1:6, names)
#clrNames = c("0.1","0.3","0.5")
#clr = factor(c(1,1,rep(1,length(trCost)), 2,2,rep(2,length(trCost)), 3,3,rep(3,length(trCost))), 
#             1:3, clrNames)
#typeNames = c("Glyc.", "Mito.")
#type = factor(c(1,1,rep(2,length(trCost)), 1,1,rep(2,length(trCost)), 1,1,rep(2,length(trCost))), 
#             1:2, typeNames)

#ds = tibble(x=allX, y=allY, Pathway=type, clr = clr)
#ds = ds[ds$y <= 1,]

#cols = c('#B5D39B','#E7B56C','#6B97BC')
#cols = c('#95B37B','#E0B064','#6B97BC')

#pX = ggplot(ds, aes(x = x, y = y, colour = clr, linetype = Pathway)) +
#  geom_line(size=1.3) +
#  scale_linetype_manual(values = c(1,2), labels = typeNames) +
#  scale_color_manual(values = c(cols[1],cols[2],cols[3]), labels = clrNames) +
#  ggplot2::labs(y="Utilization", x="Transport cost", colour="Glyc. Util.") +
#  ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#  theme(text = element_text(size=14, color="black"), axis.text = element_text(size=14, color="black"), legend.text = element_text(size=14, color="black"))
#pX
