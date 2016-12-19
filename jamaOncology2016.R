################################################################################################
# set working directory to jamaOncology_2016
userDirectory <- "~/Desktop/jamaOncology_2016-master"
setwd(userDirectory)

# prerequisite packages, be sure these are already installed using installRequiredPackages.R
requiredPackages.cran <- c(
  "devtools", 
  "reshape",  
  "RColorBrewer", 
  "ggplot2", 
  "gplots", 
  "cowplot")

requiredPackages.bioconductor <- c(
  "genefu",
  "ComplexHeatmap")

lapply(c(requiredPackages.bioconductor,requiredPackages.cran), require, character.only = TRUE)
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

################################################################################################
# functions

# http://stackoverflow.com/questions/15253954/replace-multiple-arguments-with-gsub
mgsub <- function(pattern, replacement, x, ...) {
  if (length(pattern)!=length(replacement)) {
    stop("pattern and replacement do not have the same length.")
  }
  result <- x
  for (i in 1:length(pattern)) {
    result <- gsub(pattern[i], replacement[i], result, ...)
  }
  result
}

# for tilePlot coloring, genesUp
fcColScaleUp <- function(fcValue) {
  if (fcValue < 1) {
    return(0)
  }
  if (fcValue > 1 & fcValue < 2) {
    return(1)
  }
  if (fcValue > 2) {
    return(2)
  }
}

# for tilePlot coloring, genesDown
fcColScaleDown <- function(fcValue) {
  if (fcValue > -1) {
    return(0)
  }
  if (fcValue < -1 & fcValue > -2) {
    return(1)
  }
  if (fcValue < -2) {
    return(2)
  }
}

################################################################################################
# import data

brain.met.pairs.ns.normalized <- read.table("dataInput/2016-03-10_normalized-counts.txt")
brain.met.pairs.ns.normalized.log2 <- log2(brain.met.pairs.ns.normalized)

sample.info.table <- read.table("dataInput/2016-03-10_sampleTable.txt",  sep = "\t", row.names = 1, header = TRUE, 
  stringsAsFactors = FALSE)
sample.info.table <- replace(sample.info.table, is.na(sample.info.table), "NA") 

################################################################################################
# FIGURE 1A

mydist <- function(x) {as.dist(1 - cor(t(x), use = "pa"))}
myclust <- function(x) {hclust(x,method="average")}

ER.colors <- mgsub(c("Pos","Neg", "NA"), c("darkgreen", "black", "gray80"), 
  sample.info.table$ER.Status)
HER2.colors <- mgsub(c("Pos","Neg", "NA"), c("darkgreen", "black", "gray80"), 
  sample.info.table$Her2.Status)
Institution.colors <- mgsub(c("RCS","Pitt"), c("goldenrod2","darkviolet"), 
  sample.info.table$Tissue.Source)
PrimaryMet.colors <- mgsub(c("Primary","Metastasis"), c("darkblue", "firebrick3"), 
  sample.info.table$Tumor.Type)

clab <- cbind(PrimaryMet.colors,Institution.colors,HER2.colors,ER.colors)
colnames(clab)=c("Tumor","Institution","HER2","ER")

breaks=seq(-2, 2, by=0.1)
mycol <- colorpanel(n=length(breaks)-1,low="blue3",mid="white",high="red3")

data.hm <- as.matrix(brain.met.pairs.ns.normalized.log2)

pdf(file = "dataOutput/Figure_1A.pdf", height = 12, width = 8)
  heatmap.3(data.hm, distfun = mydist, hclustfun = myclust, 
    scale = "row", dendrogram = "column", ColSideColors=clab,
    col=mycol, cexRow = 0.5, cexCol = 1.0, trace="none", key=FALSE,
    margins=c(12,12), ColSideColorsSize=3,
    keysize=1, key.title = NA, offsetRow = 0, offsetCol = 0, 
    adjCol = 1, breaks = breaks, symkey = TRUE)
dev.off()

pdf(file = "dataOutput/Figure_1A-legend.pdf", width = 5, height = 5 )
  plot(c(1,1))
  legend("top", legend=c("Positive","Negative","","RCS","Pitt","", "Primary", "Metastasis"),
    fill=c("darkgreen","black", "white","goldenrod2","darkviolet","white","darkblue","firebrick3"), 
    border=FALSE, bty="n", y.intersp = 0.8, x.intersp =0.6, cex=1.5)
dev.off()

################################################################################################
# Figure 1B

pam50.gl <- scan("dataInput/PAM50-gl.txt", what = "character")
pam50.ns <- t(brain.met.pairs.ns.normalized.log2[pam50.gl,])
PAM50.subtype = rep(NA,nrow(pam50.ns))

esr1.pos.patients <- rownames(subset(sample.info.table, ER.Status == "Pos" & Tumor.Type == "Primary"))
esr1.neg.patients <- rownames(subset(sample.info.table, ER.Status == "Neg" & Tumor.Type == "Primary"))

pam50.output.df <- data.frame(PAM50.subtype = PAM50.subtype, row.names = row.names(pam50.ns))
pam50.probability.df <- data.frame(Basal = rep(NA,40), Her2 = rep(NA,40), 
  LumA = rep(NA,40), LumB = rep(NA,40), Normal = rep(NA,40),  row.names = row.names(pam50.ns))

# PAM50 assignments
for (i in row.names(pam50.ns)) {
  print(i)
  # 20 samples of ER balanced
  test.set <- c(esr1.pos.patients,esr1.neg.patients)
  # add 1 query sample
  sample.set <- c(test.set,i)
  # generate PAM50 matrix for sample set
  pam50.sample.set.final <- pam50.ns[sample.set,]
  # call PAM50
  PAM50.subtype <- intrinsic.cluster.predict(sbt.model=pam50.robust, 
    data=pam50.sample.set.final, annot=annot.nkis, do.mapping=FALSE, 
    do.prediction.strength=FALSE, verbose=TRUE)
  # store query sample PAM50 in pam50.output.df, probability in pam50.probability.df
  pam50.output.df[i,1] <- as.character(PAM50.subtype$subtype[i])
  pam50.probability.df[i,] <- PAM50.subtype$subtype.proba[21,]
}

write.table(pam50.output.df, file = "dataOutput/genefu_pam50.brain.met.pairs-ns-samples.txt", 
  quote = FALSE, col.names = NA, sep = "\t")
write.table(signif(pam50.probability.df,3), 
  file = "dataOutput/genefu_pam50.probabilities.brain.met.pairs-ns-samples.txt", 
  quote = FALSE, col.names = NA, sep = "\t")

# PAM50 plot
pam50.subtypes <- read.table("dataOutput/genefu_pam50.brain.met.pairs-ns-samples.txt",
  sep = "\t", row.names = 1, header = TRUE, stringsAsFactors = FALSE)
pam50.subtypes$Case <- substring(rownames(pam50.subtypes),4,100)
pam50.subtypes$Tumor <- sample.info.table[rownames(pam50.subtypes),]$Tumor.Type
pam50.subtypes$PAM50.plot <- mgsub(c("Basal","Her2", "LumA", "LumB"), 
  c("1", "2", "3", "4"), pam50.subtypes$PAM50.subtype)

pam50.subtypes$Tumor2 <- factor(pam50.subtypes$Tumor, levels = c("Primary", "Metastasis"))
pam50.subtypes$Case2 <- factor(pam50.subtypes$Case, levels = rev(unique(pam50.subtypes$Case)) )

pdf(file = "dataOutput/Figure_1B.pdf", width = 1.3, height = 5 )
  qplot(data=pam50.subtypes,x=Tumor2,y=Case2,fill=factor(PAM50.plot),geom="tile") +
    scale_fill_manual(values=c("1"="firebrick3", "2"="darkgreen", "3"="darkblue", "4" = "goldenrod2")) +
    theme_classic() +
    theme(legend.position="none", 
      panel.border=element_blank(),
      plot.title=element_text(lineheight=.8, face="bold"),
      axis.text=element_text(size=5,face="bold"),
      axis.title=element_text(size=4,face="bold"),
      axis.text.x=element_text(size=5,face="bold"),
      axis.title.x=element_text(size=0,face="bold"),
      axis.title.y=element_text(size=0,face="bold")) + 
    scale_y_discrete(expand = c(0,0)) + 
    scale_x_discrete(expand = c(0,0), position = "top") + 
    ggtitle("PAM50")
dev.off()

pdf(file = "dataOutput/Figure_1B-legend.pdf", width = 5, height = 5 )
  plot(c(1,1))
  legend("top", legend=c("Basal","Her2","LumA","LumB"),
    fill=c("firebrick3","darkgreen","darkblue", "goldenrod2"),
    border=FALSE, bty="n", y.intersp = 0.8, x.intersp =0.6, cex=2)
dev.off()

################################################################################################
# eFigure 1

nanostring.log2.fc <- data.frame(matrix(NA, nrow = nrow(brain.met.pairs.ns.normalized.log2), 
  ncol = ncol(brain.met.pairs.ns.normalized.log2)/2))

caseIDs <- substr(rownames(sample.info.table),4,100)
colnames(nanostring.log2.fc) <- unique(caseIDs)
rownames(nanostring.log2.fc) <- rownames(brain.met.pairs.ns.normalized.log2)

for (i in seq(2,ncol(brain.met.pairs.ns.normalized.log2),2)) {
  caseID <- caseIDs[i]
  nanostring.log2.fc[,caseID] <-  brain.met.pairs.ns.normalized.log2[,i] - brain.met.pairs.ns.normalized.log2[,(i-1)] 
}

write.table(nanostring.log2.fc, file = "dataOutput/nanostring.log2.fc.txt", 
  sep = "\t", quote = FALSE, col.names = NA)
nanostring.log2.fc.melt <- melt(nanostring.log2.fc)
colnames(nanostring.log2.fc.melt) <- c("Case", "FoldChange")

mean(nanostring.log2.fc.melt$FoldChange)
sd(nanostring.log2.fc.melt$FoldChange)

pdf("dataOutput/eFigure_1.pdf", height = 10, width = 10)
  ggplot(nanostring.log2.fc.melt, aes(FoldChange, fill = Case)) + 
    geom_density(alpha = 0.3) + 
    ggtitle("Fold Change Value Distribution in Patient-Matched Pairs") +
    theme(plot.title = element_text(size = 15, face = "bold"),
      axis.title.x = element_text(size = 12, face = "bold"), 
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size=12, face = "bold"),
      panel.background = element_rect(fill = "gray92")) +
    scale_y_continuous(expand = c(0,0)) + 
    scale_x_continuous(expand = c(0,0)) +
    geom_vline(xintercept=-1, col = "firebrick3", lwd = 2) + 
    geom_vline(xintercept=1,col = "firebrick3", lwd = 2) + 
    xlab("Log2 Fold Change") + 
    ylab("Density")
dev.off()

################################################################################################
# Figure 2A & eFigure 2A

upregulated.matrix <- nanostring.log2.fc > 1
upregulated.matrix.ordered <- upregulated.matrix[order(-rowSums(upregulated.matrix)),]
upregulated.inAtLeast1.index <- grep("TRUE", rowSums(upregulated.matrix.ordered) >= 1)
upregulated.matrix.final <- upregulated.matrix.ordered[upregulated.inAtLeast1.index,]

downregulated.matrix <- nanostring.log2.fc < -1
downregulated.matrix.ordered <- downregulated.matrix[order(-rowSums(downregulated.matrix)),]
downregulated.inAtLeast1.index <- grep("TRUE", rowSums(downregulated.matrix.ordered) >= 1)
downregulated.matrix.final <- downregulated.matrix.ordered[downregulated.inAtLeast1.index,]

class(upregulated.matrix.final) <- "numeric"
class(downregulated.matrix.final) <- "numeric"

data.up <- upregulated.matrix.final
data.down <- downregulated.matrix.final

# format data for oncoPrint
mat.list <- list(Increased = data.up, Decreased = data.down )
mat.list.unified <- unify_mat_list(mat.list)

alter_fun_list = list(
    background = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
    },
    Decreased = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "blue3", col = NA))
    },
    Increased = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "red3", col = NA))
    }
)

col = c("Increased" = "red3", "Decreased" = "blue3")

increased.decreased <- cbind(mat.list.unified$Increased, mat.list.unified$Decreased)
recurrent.changes.genes.4 <- names(which(rowSums(increased.decreased) > 4))
recurrent.changes.genes.1 <- names(which(rowSums(increased.decreased) > 1))
mat.list.recurrent.4 <- list(Increased = mat.list.unified$Increased[recurrent.changes.genes.4,], 
  Decreased = mat.list.unified$Decreased[recurrent.changes.genes.4,])
mat.list.recurrent.1 <- list(Increased = mat.list.unified$Increased[recurrent.changes.genes.1,], 
  Decreased = mat.list.unified$Decreased[recurrent.changes.genes.1,])

pdf("dataOutput/Figure_2A.pdf", height = 10, width = 6)
  oncoPrint(mat.list.recurrent.4, get_type = function(x) strsplit(x, ";")[[1]],
    alter_fun_list = alter_fun_list, col = col, 
    column_title = "Recurrent Expression Alterations 
    (>4 pairs) in Brain Metastases",
    heatmap_legend_param = list(title = "Expression", at = c("Increased", "Decreased"), 
      labels = c("Increased", "Decreased")),
    row_names_gp = gpar(fontsize = 7,fontface = "bold"), 
    pct_gp = gpar(fontsize =5, fontface = "bold"))
dev.off()

pdf("dataOutput/eFigure_2.pdf", height = 10, width = 6)
  oncoPrint(mat.list.recurrent.1, get_type = function(x) strsplit(x, ";")[[1]],
    alter_fun_list = alter_fun_list, col = col, 
    column_title = "Recurrent Expression Alterations 
    (>1 pair) in Brain Metastases",
    heatmap_legend_param = list(title = "Expression", at = c("Increased", "Decreased"), 
      labels = c("Increased", "Decreased")),
    row_names_gp = gpar(fontsize = 7,fontface = "bold"), 
    pct_gp = gpar(fontsize =5, fontface = "bold"))
dev.off()

################################################################################################
# Figure 2B

genesOfInterest <- scan("dataInput/clinicallyActionable-gl.txt", what = "character")

nanostring.log2.fc.cA <- nanostring.log2.fc[genesOfInterest,]

upregulated.matrix.cA <- nanostring.log2.fc.cA > 1
upregulated.matrix.cA.final <- upregulated.matrix.cA[order(-rowSums(upregulated.matrix.cA)),]
geneOrderUp <- rownames(upregulated.matrix.cA.final)

downregulated.matrix.cA <- nanostring.log2.fc.cA < -1
downregulated.matrix.cA.final <- downregulated.matrix.cA[order(-rowSums(downregulated.matrix.cA)),]
geneOrderDown <- rownames(downregulated.matrix.cA.final)

# Figure 2B Top

nanostring.log2.fc.cA.ordered.up <- as.matrix(nanostring.log2.fc.cA[geneOrderUp,])
data.up.cA <- melt(nanostring.log2.fc.cA.ordered.up)
data.up.cA$X1 <- factor(data.up.cA$X1, levels = rev(data.up.cA$X1))
colnames(data.up.cA) <- c("Gene", "Case", "Expression")
data.up.cA$color.translation <- unlist(sapply(data.up.cA$Expression, FUN = fcColScaleUp))

pdf(file = "dataOutput/Figure_2B-top.pdf", width = 6, height = 4 )
  qplot(data=data.up.cA, x=Case, y=Gene, geom="tile", fill=factor(color.translation)) + 
    scale_fill_manual(values=c("0"="gray92", "1"="brown2", "2"="darkred")) + 
    theme(legend.position="none", 
      panel.border=element_blank(),
      plot.title = element_text(lineheight=.8, face="bold"),
      axis.text=element_text(size=10,face="bold"),
      axis.title=element_text(size=9,face="bold"),
      axis.text.x=element_text(size=10,face="bold", angle = -90),
      axis.title.x=element_text(size=0,face="bold"),
      axis.title.y=element_text(size=0,face="bold")) + 
    scale_y_discrete(expand = c(0,0)) + 
    scale_x_discrete(expand = c(0,0), position = "top")
dev.off()

pdf(file = "dataOutput/Figure_2B-top-legend.pdf", width = 6, height = 6 )
  plot(1,1)
  legend("top", legend=c(">2-Fold Increase",">4-Fold Increase"),
    fill=c("brown2","darkred"), 
    border=FALSE, bty="n", y.intersp = 0.8, x.intersp =0.6, cex=1.5)
dev.off()

# Figure 2B Bottom

nanostring.log2.fc.cA.ordered.down <- as.matrix(nanostring.log2.fc.cA[geneOrderDown,])
data.down.cA <- melt(nanostring.log2.fc.cA.ordered.down)
data.down.cA$X1 <- factor(data.down.cA$X1, levels = rev(data.down.cA$X1))
colnames(data.down.cA) <- c("Gene", "Case", "Expression")

data.down.cA$color.translation <- unlist(sapply(data.down.cA$Expression, FUN = fcColScaleDown))

pdf(file = "dataOutput/Figure_2B-bottom.pdf", width = 6, height = 4 )
  qplot(data=data.down.cA, x=Case, y=Gene, geom="tile", fill=factor(color.translation)) + 
    scale_fill_manual(values=c("0"="gray92", "1"="dodgerblue2", "2"="darkblue")) + 
    theme(legend.position="none", 
      panel.border=element_blank(),
      plot.title = element_text(lineheight=.8, face="bold"),
      axis.text=element_text(size=10,face="bold"),
      axis.title=element_text(size=9,face="bold"),
      axis.text.x=element_text(size=10,face="bold", angle = -90),
      axis.title.x=element_text(size=0,face="bold"),
      axis.title.y=element_text(size=0,face="bold")) + 
    scale_y_discrete(expand = c(0,0)) + 
    scale_x_discrete(expand = c(0,0), position = "top")
dev.off()

pdf(file = "dataOutput/Figure_2B-bottom-legend.pdf", width = 6, height = 6 )
  plot(1,1)
  legend("top", legend=c(">2-Fold Decrease",">4-Fold Decrease"),
    fill=c("dodgerblue2","darkblue"), 
    border=FALSE, bty="n", y.intersp = 0.8, x.intersp =0.6, cex=1.5)
dev.off()

################################################################################################
# Figure 3A

primary.tumors <- rownames(subset(sample.info.table, Tumor.Type == "Primary"))
metastatic.tumors <- rownames(subset(sample.info.table, Tumor.Type == "Metastasis"))

geneOfInterest = "ERBB2"

paired.plot.data <- data.frame(Patient = caseIDs, 
  Gene = t(brain.met.pairs.ns.normalized.log2[geneOfInterest,]), 
  Tumor.Type = sample.info.table$Tumor.Type)

paired.plot.data$Tumor.Type <- gsub("Primary", "Primary Tumors", paired.plot.data$Tumor.Type)
paired.plot.data$Tumor.Type <- gsub("Metastasis", "Brain Metastases", paired.plot.data$Tumor.Type)

paired.plot.data$Tumor.Type <- factor(paired.plot.data$Tumor.Type, as.character(paired.plot.data$Tumor.Type))
paired.plot.data$Tumor.Type <- factor(paired.plot.data$Tumor.Type, c("Primary Tumors", "Brain Metastases"))

wilcox.result.signedRank <- wilcox.test(paired.plot.data[primary.tumors,2], 
  paired.plot.data[metastatic.tumors,2], paired = TRUE)
wilcoxSignedRank.table.p.value <- wilcox.result.signedRank$p.value
wilcoxSignedRank.table.p.value.text <- as.character(signif(wilcoxSignedRank.table.p.value,1))

her2.neg.inc <- subset(paired.plot.data, Patient == "RCS_4" | Patient == "Pitt_62" | Patient == "RCS_6")
other.her2 <- subset(paired.plot.data, Patient != "RCS_4" & Patient != "Pitt_62" & Patient != "RCS_6")

pdf(file = "dataOutput/Figure_3A.pdf", height = 6, width = 6)
  ggplot(other.her2, aes_string(x="Tumor.Type", y=geneOfInterest, group="Patient", label="Patient")) +
  ggtitle(paste(geneOfInterest, sep = " ")) +
  theme(plot.title = element_text(size = 25, face = "bold"),
    axis.title.x = element_text(size = 0, face = "bold"), 
    axis.text.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size=18, face = "bold"),
    panel.background = element_rect(fill = "gray92")) + 
  geom_point(data=her2.neg.inc, col = "forestgreen",size=8.0, alpha = 0.8, position=position_dodge(width=0.2)) +
  geom_line(data=her2.neg.inc,size=0.5, alpha=1.0, position=position_dodge(width=0.2)) +
  geom_point(aes(colour=Tumor.Type), size=7.0, alpha = 0.8, position=position_dodge(width=0.15)) +
  geom_line(size=0.2, alpha=1.0, position=position_dodge(width=0.15)) +
  scale_colour_manual(values=c("blue3", "red3"), guide=FALSE) +
  xlab('Tumor Type') +
  ylab(paste("Log2 Normalized Counts")) + 
  scale_shape_manual(values = c(17, 15,19)) +
  geom_text(x = Inf, y = Inf, label = paste("P=",wilcoxSignedRank.table.p.value.text, sep = " "), 
    hjust = 1, vjust = 1.1, size = 7)
dev.off()

########################################################################################################
# eFigure 3

geneOfInterest = "ESR1"

paired.plot.data <- data.frame(Patient = caseIDs, 
  Gene = t(brain.met.pairs.ns.normalized.log2[geneOfInterest,]), 
  Tumor.Type = sample.info.table$Tumor.Type)

paired.plot.data$Tumor.Type <- gsub("Primary", "Primary Tumors", paired.plot.data$Tumor.Type)
paired.plot.data$Tumor.Type <- gsub("Metastasis", "Brain Metastases", paired.plot.data$Tumor.Type)
paired.plot.data$Tumor.Type <- factor(paired.plot.data$Tumor.Type, as.character(paired.plot.data$Tumor.Type))
paired.plot.data$Tumor.Type <- factor(paired.plot.data$Tumor.Type, c("Primary Tumors", "Brain Metastases"))

wilcox.result.signedRank <- wilcox.test(paired.plot.data[primary.tumors,2], 
  paired.plot.data[metastatic.tumors,2], paired = TRUE)
wilcoxSignedRank.table.p.value <- wilcox.result.signedRank$p.value
wilcoxSignedRank.table.p.value.text <- as.character(signif(wilcoxSignedRank.table.p.value,1))

er.pos.dec <- subset(paired.plot.data, Patient == "RCS_2" | Patient == "Pitt_47")
other.er <- subset(paired.plot.data, Patient != "RCS_2" & Patient != "Pitt_47")

pdf(file = "dataOutput/eFigure_3.pdf", height = 6, width = 6)
  ggplot(other.er, aes_string(x="Tumor.Type", y=geneOfInterest, group="Patient", label="Patient")) +
  ggtitle(paste(geneOfInterest, sep = " ")) +
  theme(plot.title = element_text(size = 25, face = "bold"),
    axis.title.x = element_text(size = 0, face = "bold"), 
    axis.text.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size=18, face = "bold"),
    panel.background = element_rect(fill = "gray92")) + 
  geom_point(data=er.pos.dec, col = "forestgreen",size=8.0, alpha = 0.8, position=position_dodge(width=0.2)) +
  geom_line(data=er.pos.dec,size=0.5, alpha=1.0, position=position_dodge(width=0.2)) +
  geom_point(aes(colour=Tumor.Type), size=7.0, alpha = 0.8, position=position_dodge(width=0.15)) +
  geom_line(size=0.2, alpha=1.0, position=position_dodge(width=0.15)) +
  scale_colour_manual(values=c("blue3", "red3"), guide=FALSE) +
  xlab('Tumor Type') +
  ylab(paste("Log2 Normalized Counts")) + 
  scale_shape_manual(values = c(17, 15,19)) +
  geom_text(x = Inf, y = Inf, label = paste("P=",wilcoxSignedRank.table.p.value.text, sep = " "), 
    hjust = 1, vjust = 1.1, size = 7)
dev.off()


########################################################################################################
# eFigure 2B

geneList <- c("KRT14", "KRT5", "RAB6B", "GRB7")

pdf(file = "dataOutput/eFigure_2B.pdf", height = 6, width = 6) 
  for (geneOfInterest in geneList) {
    
    paired.plot.data <- data.frame(Patient = caseIDs, 
      Gene = t(brain.met.pairs.ns.normalized.log2[geneOfInterest,]), 
      Tumor.Type = sample.info.table$Tumor.Type)

    paired.plot.data$Tumor.Type <- gsub("Primary", "Primary Tumors", paired.plot.data$Tumor.Type)
    paired.plot.data$Tumor.Type <- gsub("Metastasis", "Brain Metastases", paired.plot.data$Tumor.Type)
    paired.plot.data$Tumor.Type <- factor(paired.plot.data$Tumor.Type, as.character(paired.plot.data$Tumor.Type))
    paired.plot.data$Tumor.Type <- factor(paired.plot.data$Tumor.Type, c("Primary Tumors", "Brain Metastases"))

    wilcox.result.signedRank <- wilcox.test(paired.plot.data[primary.tumors,2], 
      paired.plot.data[metastatic.tumors,2], paired = TRUE)
    wilcoxSignedRank.table.p.value <- wilcox.result.signedRank$p.value
    wilcoxSignedRank.table.p.value.text <- as.character(signif(wilcoxSignedRank.table.p.value,1))

    print(ggplot(paired.plot.data, aes_string(x="Tumor.Type", y=geneOfInterest, group="Patient", label="Patient")) +
      ggtitle(paste(geneOfInterest, sep = " ")) +
      theme(plot.title = element_text(size = 25, face = "bold"),
        axis.title.x = element_text(size = 0, face = "bold"), 
        axis.text.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size=18, face = "bold"),
        panel.background = element_rect(fill = "gray92")) + 
      geom_point(aes(colour=Tumor.Type), size=7.0, alpha = 0.8, position=position_dodge(width=0.15)) +
      geom_line(size=0.2, alpha=1.0, position=position_dodge(width=0.15)) +
      scale_colour_manual(values=c("blue3", "red3"), guide=FALSE) +
      xlab('Tumor Type') +
      ylab(paste("Log2 Normalized Counts")) + 
      scale_shape_manual(values = c(17, 15,19)) +
      geom_text(x = Inf, y = Inf, label = paste("P=",wilcoxSignedRank.table.p.value.text, sep = " "), 
        hjust = 1, vjust = 1.1, size = 7)

    )

  }
 
dev.off()

########################################################################################################
# end
