# TL - 060218

# Ce script permet de verifier la qualité des sites sur les variants appelés chez Colivacea

library("ggplot2")
library("grid")

covqual11pc=read.table("~/Collab/Montpellier/Colivacea_CovQual_joinVcf/Colivacea_FiltrationQualityjointVCF_variantonly_allsites_060218.PROV.30pc.sed.withoutNA", header=TRUE)
covqualind=read.table("~/Collab/Montpellier/Colivacea_CovQual_joinVcf/Colivacea_bidon", header=FALSE)

pdf("C_olivacea_summary_calling_newFilteredGVCF_060218.pdf",12,8)

###### SNP QUALITY
### DP / MQ
# Nuage de points colorés par groupes
scatterPlot <- ggplot(covqual11pc,aes(x=DP, y=MQ)) + xlim(0,200) + ylim(0,58)+
  xlab("Coverage (DP)") + ylab("Mapping Quality (MQ)") +
  geom_density_2d(aes(colour=label),geom = "polygon") + 
  #scale_colour_manual(values = c('#999999','#E69F00',"red")) + 
  theme(legend.position=c(0.85,0.85))
# Courbe de densité marginale de x (panel du haut)
xdensity <- ggplot(covqual11pc, aes(DP, colour=label)) + xlim(0,200)+
  geom_density() + xlab("Coverage (DP)") #+
#scale_fill_manual(values = c('#999999','#E69F00')) + 
#theme(legend.position = "none")
# Courbe de densité marginale de y (panel de droite)
ydensity <- ggplot(covqual11pc, aes(MQ, colour=label)) + xlim(0,70)+
  geom_density() + xlab("Mapping Quality (MQ)")+
  #scale_fill_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position = "none")

# Nouvelle page
grid.newpage()
# Créer la mise en page : nrow = 2, ncol = 2
pushViewport(viewport(layout = grid.layout(2, 2)))
# Une fonction pour definir une region dans la mise en page
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
# Arranger les graphiques
print(scatterPlot, vp=define_region(1, 1:2))
print(xdensity, vp = define_region(2, 1))
print(ydensity, vp = define_region(2, 2))

### All Freq / DP
# Nuage de points colorés par groupes
scatterPlot <- ggplot(covqual11pc,aes(x=AF, y=DP)) + xlim(0,1) + ylim(0,200)+
  xlab("Allele Frequency (AF)") + ylab("Coverage (DP)") +
  geom_density_2d(aes(colour=label),geom = "polygon") + 
  #scale_colour_manual(values = c('#999999','#E69F00',"red")) + 
  theme(legend.position=c(0.85,0.85))
# Courbe de densité marginale de x (panel du haut)
xdensity <- ggplot(covqual11pc, aes(AF, colour=label)) + xlim(0,1)+
  geom_density() + xlab("Allele Frequency (AF)") #+
#scale_fill_manual(values = c('#999999','#E69F00')) + 
#theme(legend.position = "none")
# Courbe de densité marginale de y (panel de droite)
ydensity <- ggplot(covqual11pc, aes(DP, colour=label)) + xlim(0,200)+
  geom_density() + xlab("Coverage (DP)")+
  #scale_fill_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position = "none")

# Nouvelle page
grid.newpage()
# Créer la mise en page : nrow = 2, ncol = 2
pushViewport(viewport(layout = grid.layout(2, 2)))
# Une fonction pour definir une region dans la mise en page
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
# Arranger les graphiques
print(scatterPlot, vp=define_region(1, 1:2))
print(xdensity, vp = define_region(2, 1))
print(ydensity, vp = define_region(2, 2))


### QD / Mapping quality

# Nuage de points colorés par groupes
scatterPlot <- ggplot(covqual11pc,aes(x=QD, y=MQ)) + xlim(0,60) + ylim(0,58)+
  xlab("Quality by Depth (QD)") + ylab("Mapping Quality (MQ)") +
  geom_density_2d(aes(colour=label),geom = "polygon") + 
  #scale_colour_manual(values = c('#999999','#E69F00',"red")) + 
  theme(legend.position=c(0.85,0.85))
# Courbe de densité marginale de x (panel du haut)
xdensity <- ggplot(covqual11pc, aes(QD, colour=label)) + xlim(0,60)+
  geom_density() + xlab("Quality by Depth (DP)") #+
#scale_fill_manual(values = c('#999999','#E69F00')) + 
#theme(legend.position = "none")
# Courbe de densité marginale de y (panel de droite)
ydensity <- ggplot(covqual11pc, aes(MQ, colour=label)) + xlim(0,70)+
  geom_density() + xlab("Mapping Quality (MQ)")+
  #scale_fill_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position = "none")

# Nouvelle page
grid.newpage()
# Créer la mise en page : nrow = 2, ncol = 2
pushViewport(viewport(layout = grid.layout(2, 2)))
# Une fonction pour definir une region dans la mise en page
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
# Arranger les graphiques
print(scatterPlot, vp=define_region(1, 1:2))
print(xdensity, vp = define_region(2, 1))
print(ydensity, vp = define_region(2, 2))

###### COVERAGE PER IND
ggplot(covqualind, aes(V4, colour=V1)) + xlim(0,100)+
    geom_density() + xlab("Coverage (DP) per individual")+
    #scale_fill_manual(values = c('#999999','#E69F00')) + 
    theme(legend.position = c(0.8,0.8))

dev.off()

#summary(lm(formula=covqual$QUAL ~ covqual$QD))
#summary(lm(formula=covqual$QD ~ I(covqual$QUAL/covqual$DP)))


quantiles_99 <- function(x) {
  r <- quantile(x, probs=c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99))
  names(r) <- c("0.01", "0.05", "0.1","0.25", "median", "0.75", "0.9","0.95","0.99")
  r
}


z <- aggregate(DP ~ label, covqual11pc, function(x) c(quant = round(quantiles_99(x),0)))
zind <- aggregate(V4 ~ V1, covqualind, function(x) c(quant = round(quantiles_99(x),0)))

z <- data.frame(z)
zind <- data.frame(zind)

write.table(noquote(z), file = "summary_statistics_coverage.txt", sep="\t")
write.table(noquote(zind), file = "summary_statistics_coverage_ind.txt", sep= "\t")

