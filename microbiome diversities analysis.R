#environment
library(tidyverse)
library(tidytext)
library(vegan)
library(ape)
library(philr)
library(ggrepel)
library(phyloseq)
library(rbiom)
library(microbiome)
library(RColorBrewer)
library(VennDiagram)
library(Maaslin2)
library(plotly)
library(forestplot)
library(forcats)
library(factoextra)

#import metadata
metadata <- read.table("Metadata.txt",header=T)

####alpha diversity ####
  ##import taxa estimated absolute counts####
Abs.Filenames <- list.files("./Abs", full.names = T)

Abs <- lapply(Abs.Filenames, function(x) {
  df <- suppressMessages(read_tsv(x, col_names = T, skip=5) %>% dplyr::select(c(1,5)))
  col <- x %>% basename %>% gsub (".abs.txt","",.) %>% gsub ("337","W14",.) %>% gsub("352","W4",.)
  colnames(df) <- c("taxa",col)
  return(df) })

Abs.all <- purrr:: reduce (Abs, full_join, by="taxa")

#export all samples absolute counts
write.csv(Abs.all,"Abs.all.csv")

rm(Abs.all, Abs.Filenames, Abs)

  ##import abs SGB table####
abs.spe <- read.table("SGB.abs.txt",header=T,row.names = 1)
abs.spe[is.na(abs.spe)] <- 0


  ##rarefy to min counts####
abs.min <- min(colSums(abs.spe))
abs.spe.rare <-vegan::rrarefy(t(abs.spe), abs.min) %>% t()

  ##Rarefaction plot by ggplot2####

out <- rarecurve(t(abs.spe), step = 6000, sample = abs.min, col="blue",cex=0.6)
names(out) <-  colnames(abs.spe)

protox <- mapply(FUN = function(x, y) {
  mydf <- as.data.frame(x)
  colnames(mydf) <- "value"
  mydf$SampleID <- y
  mydf$subsample <- attr(x, "Subsample")
  mydf
}, x = out, y = as.list(names(out)), SIMPLIFY = F)

xy <- do.call(rbind, protox)
row.names(xy) <- NULL

xy <- data.frame(xy, metadata[match(xy$SampleID, metadata$ID),])

seq.no1 <-seq(1:5)
seq.no2 <-seq(1:6) 
seq.no3 <-seq(1:9) 
seq.no4 <-seq(1:7) 

xy <- xy %>% mutate(Sample_ID=factor(Sample_ID, 
                                     levels=c(paste("W4.WT",seq.no1,sep="."), 
                                              paste("W4.CF",seq.no2,sep="."),
                                              paste("W14.WT",seq.no3,sep="."),
                                              paste("W14.CF",seq.no4,sep=".")))) 

mycolors2 <- colorRampPalette(c("blue4","mediumpurple", "olivedrab",
                                "mediumspringgreen" ,"firebrick","orange","gold"))(27)

tiff("rarecurve1.tiff", unit="in", width=6.5, height =4.5, res=600, compression = 'lzw' )
ggplot(xy, aes(x = subsample, y = value, color = Sample_ID)) +
  theme_bw() + 
  geom_line(size=0.8) +
  geom_vline(xintercept=abs.min, color= "red", linetype='dashed') + 
  theme(axis.text= element_text(color="black",size=12),
        axis.title=element_text(size=12))+
  xlab("Sample Size") + ylab('Species')+
  scale_color_manual(values=mycolors2)+ggtitle("Rarefaction Curve")+
  theme(plot.title = element_text(hjust = 0.5,size=14))+xlim(0,6e7)
dev.off()

rm(xy,out,protox,seq.no1,seq.no2,seq.no3,seq.no4,mycolors2,abs.min)

  ##alpha diversity all samples together####
alpha <- microbiome::alpha(abs.spe.rare) %>% 
  rownames_to_column (var="ID") %>% merge(metadata,by="ID") %>% column_to_rownames(var="ID")
write.csv(alpha,"alpha_diversity.csv")
rm(abs.spe.rare)

  ####weaning age rarefy and alpha diversity analysis####
    ##W4 subset the dataset####
abs.spe.W4 <- abs.spe %>% dplyr::select(contains("W4"))

    ##W4 rarefy to min counts####
abs.min <- min(colSums(abs.spe.W4))
abs.spe.W4.rare <-vegan::rrarefy(t(abs.spe.W4), abs.min) %>% t()

    ##W4 rare plot by ggplot2####

out <- rarecurve(t(abs.spe.W4), step = 6000, sample = abs.min, col="blue",cex=0.6)
names(out) <-  colnames(abs.spe.W4)

protox <- mapply(FUN = function(x, y) {
  mydf <- as.data.frame(x)
  colnames(mydf) <- "value"
  mydf$SampleID <- y
  mydf$subsample <- attr(x, "Subsample")
  mydf
}, x = out, y = as.list(names(out)), SIMPLIFY = F)

xy <- do.call(rbind, protox)
row.names(xy) <- NULL

xy <- data.frame(xy, metadata[match(xy$SampleID, metadata$ID),])

seq.no1 <-seq(1:5)
seq.no2 <-seq(1:6) 

xy <- xy %>% mutate(Sample_ID=factor(Sample_ID, 
                                     levels=c(paste("W4.WT",seq.no1,sep="."), 
                                              paste("W4.CF",seq.no2,sep=".")))) 

mycolors3 <- colorRampPalette(c("blue4","mediumpurple", "olivedrab",
                                "mediumspringgreen" ,"firebrick","orange","gold"))(11)

tiff("rarecurve2.tiff", unit="in", width=6, height =4.5, res=600, compression = 'lzw' )
ggplot(xy, aes(x = subsample, y = value, color = Sample_ID)) +
  theme_bw() + 
  geom_line(size=0.8) +
  geom_vline(xintercept=abs.min, color= "red", linetype='dashed') + 
  theme(axis.text= element_text(color="black",size=14),
        axis.title=element_text(size=14))+
  xlab("Sample Size") + ylab('Species')+
  scale_color_manual(values=mycolors3)+ggtitle("4-week-old Rarefaction Curve")+
  theme(plot.title = element_text(hjust = 0.5,size=16))
dev.off()

rm(xy,out,protox,seq.no1,seq.no2,mycolors3,abs.min)

    ##W4 alpha diversity####
alpha.W4 <- microbiome::alpha(abs.spe.W4.rare) %>% 
  rownames_to_column (var="ID") %>% merge(metadata,by="ID") %>% column_to_rownames(var="ID")

write.csv(alpha.W4,"alpha_diversity_W4.csv")

rm(abs.spe.W4.rare)

  ####Adult rarefy and alpha diversity analysis####
    ##W14 subset the dataset####
abs.spe.W14 <- abs.spe %>% dplyr::select(contains("W14"))

    ##W14 rarefy to min counts####
abs.min <- min(colSums(abs.spe.W14))
abs.spe.W14.rare <-vegan::rrarefy(t(abs.spe.W14), abs.min) %>% t()

    ##W14 rare plot by ggplot2####

out <- rarecurve(t(abs.spe.W14), step = 6000, sample = abs.min, col="blue",cex=0.6)
names(out) <-  colnames(abs.spe.W14)

protox <- mapply(FUN = function(x, y) {
  mydf <- as.data.frame(x)
  colnames(mydf) <- "value"
  mydf$SampleID <- y
  mydf$subsample <- attr(x, "Subsample")
  mydf
}, x = out, y = as.list(names(out)), SIMPLIFY = F)

xy <- do.call(rbind, protox)
row.names(xy) <- NULL

xy <- data.frame(xy, metadata[match(xy$SampleID, metadata$ID),])

seq.no1 <-seq(1:9)
seq.no2 <-seq(1:7) 

xy <- xy %>% mutate(Sample_ID=factor(Sample_ID, 
                                     levels=c(paste("W14.WT",seq.no1,sep="."), 
                                              paste("W14.CF",seq.no2,sep=".")))) 

mycolors3 <- colorRampPalette(c("blue4","mediumpurple", "olivedrab",
                                "mediumspringgreen" ,"firebrick","orange","gold"))(16)

tiff("rarecurve3.tiff", unit="in", width=6, height =4.5, res=600, compression = 'lzw' )
ggplot(xy, aes(x = subsample, y = value, color = Sample_ID)) +
  theme_bw() + 
  geom_line(size=0.8) +
  geom_vline(xintercept=abs.min, color= "red", linetype='dashed') + 
  theme(axis.text= element_text(color="black",size=14),
        axis.title=element_text(size=14))+
  xlab("Sample Size") + ylab('Species')+
  scale_color_manual(values=mycolors3)+ggtitle("14-week-old Rarefaction Curve")+
  theme(plot.title = element_text(hjust = 0.5,size=16))+
  xlim(0,6e7)
dev.off()

rm(xy,out,protox,seq.no1,seq.no2,mycolors3,abs.min)

    ##W14 alpha diversity####
alpha.W14 <- microbiome::alpha(abs.spe.W14.rare) %>% 
  rownames_to_column (var="ID") %>% merge(metadata,by="ID") %>% column_to_rownames(var="ID")

write.csv(alpha.W14,"alpha_diversity_W14.csv")

rm(abs.spe.W14.rare)



####Import taxa relative abundance####
Abs.Filenames <- list.files("./Abs", full.names = T)

Rel <- lapply(Abs.Filenames, function(x) {
  df <- suppressMessages(read_tsv(x, col_names = T, skip=5) %>% dplyr::select(c(1,3)))
  col <- x %>% basename %>% gsub (".abs.txt","",.) %>% gsub ("337","W14",.) %>% gsub("352","W4",.)
  colnames(df) <- c("taxa",col)
  return(df) })

Rel.all <- purrr:: reduce (Rel, full_join, by="taxa")

#export all samples relative abundance
write.csv(Rel.all,"Rel.all.csv")
rm(Rel.all, Abs.Filenames, Rel)

####import SGB.rel table####
SGB <- read.table("SGB.rel.txt",header=T,row.names = 1) %>% t()
SGB[is.na(SGB)] <- 0

####PCoA all samples together####
  ##bray####

dist_bray <- vegan::vegdist(SGB, method = "bray")

pcoa_bray <- ape::pcoa(dist_bray)
pcoa_bray.df <- data.frame(pcoa1 = pcoa_bray$vectors[,1], 
                           pcoa2 = pcoa_bray$vectors[,2],
                           pcoa3 = pcoa_bray$vectors[,3])
pcoa_bray.df <- pcoa_bray.df %>% rownames_to_column(var="ID") %>% 
  merge(metadata,by="ID") %>% mutate(Genotype_Age=factor(Genotype_Age,level=c("WT_4wk","CF_4wk",
                                                                              "WT_14wk","CF_14wk")))

pcoa_bray.pc <- pcoa_bray$values %>% as.data.frame %>% 
  rownames_to_column("Axis") %>% mutate(PCvar=round (100*Relative_eig,2))

#bray plot#
#bray plot axis1 axis2#
pcoa_bray_plot <- ggplot(pcoa_bray.df, 
                         aes(x=pcoa1, y=pcoa2,color=Genotype_Age)) +
  geom_point(size=4) + theme_minimal() +
  xlab(paste("PCoA1"," (", paste(pcoa_bray.pc$PCvar[1]), "%)"))+
  ylab(paste("PCoA2"," (", paste(pcoa_bray.pc$PCvar[2]), "%)"))+
  scale_color_manual(values=c("#056943","#585858","#3951A2","#A80326"))+
  theme(legend.key = element_rect(fill = NA, colour = NA, size = 2),
        legend.title = element_text(color="black",size=12.5),
        legend.text = element_text(color="black",size=11.5),
        panel.background = element_blank(),
        plot.title = element_text(size=14.5, hjust = 0.5, face="bold"), 
        axis.title= element_text(color="black",size=12.5, face="bold"), 
        axis.text= element_text(color="black",size=11.5, face="bold"), 
        panel.border = element_rect(colour = "black", fill=NA))+
  ggtitle("Bray Curtis")

tiff("Figures/pcoa.bray.tiff", units="in", width = 5, height=4, res=600, compression="lzw")
pcoa_bray_plot
dev.off()

tiff("Figures/pcoa.bray2.tiff", units="in", width = 5, height=4, res=600, compression="lzw")
pcoa_bray_plot+stat_ellipse(level=0.95)
dev.off()

#bray plot axis2 axis3#
pcoa_bray_plot23 <- ggplot(pcoa_bray.df, 
                         aes(x=pcoa2, y=pcoa3,color=Genotype_Age)) +
  geom_point(size=4) + theme_minimal() +
  xlab(paste("PCoA2"," (", paste(pcoa_bray.pc$PCvar[2]), "%)"))+
  ylab(paste("PCoA3"," (", paste(pcoa_bray.pc$PCvar[3]), "%)"))+
  scale_color_manual(values=c("#056943","#585858","#3951A2","#A80326"))+
  theme(legend.key = element_rect(fill = NA, colour = NA, size = 2),
        legend.title = element_text(color="black",size=12.5),
        legend.text = element_text(color="black",size=11.5),
        panel.background = element_blank(),
        plot.title = element_text(size=14.5, hjust = 0.5, face="bold"), 
        axis.title= element_text(color="black",size=12.5, face="bold"), 
        axis.text= element_text(color="black",size=11.5, face="bold"), 
        panel.border = element_rect(colour = "black", fill=NA))+
  ggtitle("Bray Curtis")

tiff("Figures/pcoa.brayA2A3.tiff", units="in", width = 5, height=4, res=600, compression="lzw")
pcoa_bray_plot23
dev.off()

tiff("Figures/pcoa.brayA2A32.tiff", units="in", width = 5, height=4, res=600, compression="lzw")
pcoa_bray_plot23+stat_ellipse(level=0.95)
dev.off()




  ##jaccard ####
dist_jac <- vegan::vegdist(SGB, method = "jaccard")

pcoa_jac <- ape::pcoa(dist_jac)
pcoa_jac.df <- data.frame(pcoa1 = pcoa_jac$vectors[,1], 
                          pcoa2 = pcoa_jac$vectors[,2])
pcoa_jac.df <- pcoa_jac.df %>% rownames_to_column(var="ID") %>% 
  merge(metadata,by="ID") %>% mutate(Genotype_Age=factor(Genotype_Age,level=c("WT_4wk","CF_4wk",
                                                                              "WT_14wk","CF_14wk")))

pcoa_jac.pc <- pcoa_jac$values %>% as.data.frame %>% 
  rownames_to_column("Axis") %>% mutate(PCvar=round (100*Relative_eig,2))

#jac plot#
pcoa_jac_plot <- ggplot(pcoa_jac.df, 
                        aes(x=pcoa1, y=pcoa2,color=Genotype_Age)) +
  geom_point(size=4) + theme_minimal() +
  xlab(paste("PCoA1"," (", paste(pcoa_jac.pc$PCvar[1]), "%)"))+
  ylab(paste("PCoA2"," (", paste(pcoa_jac.pc$PCvar[2]), "%)"))+
  scale_color_manual(values=c("#056943","#585858","#3951A2","#A80326"))+
  theme(legend.key = element_rect(fill = NA, colour = NA, size = 2),
        legend.title = element_text(color="black",size=12.5),
        legend.text = element_text(color="black",size=11.5),
        panel.background = element_blank(),
        plot.title = element_text(size=14.5, hjust = 0.5, face="bold"), 
        axis.title= element_text(color="black",size=12.5, face="bold"), 
        axis.text= element_text(color="black",size=11.5, face="bold"), 
        panel.border = element_rect(colour = "black", fill=NA))+
  ggtitle("Jaccard")

tiff("Figures/pcoa.jac.tiff", units="in", width = 5, height=4, res=600, compression="lzw")
pcoa_jac_plot
dev.off()

tiff("Figures/pcoa.jac2.tiff", units="in", width = 5, height=4, res=600, compression="lzw")
pcoa_jac_plot+stat_ellipse(level=0.95)
dev.off()

  ##construct tree and SGB table for uniFrac analysis ####

#import phylo tree
tree <- ape::read.tree("mpa_vOct22_CHOCOPhlAnSGB_202212.nwk")

#construct tree and SGB table for uniFrac analysis
SGB_UniFrac <- t(SGB) %>% as.data.frame() %>% rownames_to_column(var="SGB")
SGB_UniFrac[,1] <- gsub("t__SGB", "", SGB_UniFrac[,1])
SGB_UniFrac <- SGB_UniFrac %>% column_to_rownames(var="SGB")
SGB_removed <- setdiff(rownames(SGB_UniFrac), tree$tip.label)
tree_filt <- ape::keep.tip(tree, setdiff(rownames(SGB_UniFrac), SGB_removed))
SGB_UniFrac <- SGB_UniFrac[tree_filt$tip.label,]
rm(tree,SGB_removed)

  ##weighted unifrac####
dist_w.UniFrac <- rbiom::beta.div(as.matrix(SGB_UniFrac), 
                                  tree=tree_filt, method="unifrac", weighted=T)

pcoa_w.UniFrac <- ape::pcoa(dist_w.UniFrac)

pcoa_w.UniFrac.df <- data.frame(pcoa1 = pcoa_w.UniFrac$vectors[,1], 
                                pcoa2 = pcoa_w.UniFrac$vectors[,2],
                                pcoa3 = pcoa_w.UniFrac$vectors[,3])

pcoa_w.UniFrac.df <- pcoa_w.UniFrac.df %>% rownames_to_column(var="ID") %>% 
  merge(metadata,by="ID") %>% mutate(Genotype_Age=factor(Genotype_Age,level=c("WT_4wk","CF_4wk",
                                                                              "WT_14wk","CF_14wk")))


pcoa_w.UniFrac.pc <- pcoa_w.UniFrac$values %>% as.data.frame %>% 
  rownames_to_column("Axis") %>% mutate(PCvar=round (100*Relative_eig,2))

#plot for weighted unifrac
pcoa_w.UniFrac_plot <- ggplot(pcoa_w.UniFrac.df, 
                              aes(x=pcoa1, y=pcoa2,color=Genotype_Age)) +
  geom_point(size=4) + theme_minimal() +
  xlab(paste("PCoA1"," (", paste(pcoa_w.UniFrac.pc$PCvar[1]), "%)"))+
  ylab(paste("PCoA2"," (", paste(pcoa_w.UniFrac.pc$PCvar[2]), "%)"))+
  scale_color_manual(values=c("#056943","#585858","#3951A2","#A80326"))+
  theme(legend.key = element_rect(fill = NA, colour = NA, size = 2),
        legend.title = element_text(color="black",size=12.5),
        legend.text = element_text(color="black",size=11.5),
        panel.background = element_blank(),
        plot.title = element_text(size=14.5, hjust = 0.5, face="bold"), 
        axis.title= element_text(color="black",size=12.5, face="bold"), 
        axis.text= element_text(color="black",size=11.5, face="bold"), 
        panel.border = element_rect(colour = "black", fill=NA))+
  ggtitle("Weighted UniFrac")

tiff("Figures/pcoa.w.UniFrac.tiff", units="in", width = 5, height=4, res=600, compression="lzw")
pcoa_w.UniFrac_plot
dev.off()

tiff("Figures/pcoa.w.UniFrac.2.tiff", units="in", width = 5, height=4, res=600, compression="lzw")
pcoa_w.UniFrac_plot+stat_ellipse(level=0.95)
dev.off()

  ##unweighted unifrac####
dist_un.UniFrac <- rbiom::beta.div(as.matrix(SGB_UniFrac), 
                                   tree=tree_filt, method="unifrac", weighted=F)

pcoa_un.UniFrac <- ape::pcoa(dist_un.UniFrac)

pcoa_un.UniFrac.df <- data.frame(pcoa1 = pcoa_un.UniFrac$vectors[,1], 
                                 pcoa2 = pcoa_un.UniFrac$vectors[,2],
                                 pcoa3 = pcoa_un.UniFrac$vectors[,3])

pcoa_un.UniFrac.df <- pcoa_un.UniFrac.df %>% rownames_to_column(var="ID") %>% 
  merge(metadata,by="ID") %>% mutate(Genotype_Age=factor(Genotype_Age,level=c("WT_4wk","CF_4wk",
                                                                              "WT_14wk","CF_14wk")))


pcoa_un.UniFrac.pc <- pcoa_un.UniFrac$values %>% as.data.frame %>% 
  rownames_to_column("Axis") %>% mutate(PCvar=round (100*Relative_eig,2))

#plot for weighted unifrac
pcoa_un.UniFrac_plot <- ggplot(pcoa_un.UniFrac.df, 
                               aes(x=pcoa1, y=pcoa2,color=Genotype_Age)) +
  geom_point(size=4) + theme_minimal() +
  xlab(paste("PCoA1"," (", paste(pcoa_un.UniFrac.pc$PCvar[1]), "%)"))+
  ylab(paste("PCoA2"," (", paste(pcoa_un.UniFrac.pc$PCvar[2]), "%)"))+
  scale_color_manual(values=c("#056943","#585858","#3951A2","#A80326"))+
  theme(legend.key = element_rect(fill = NA, colour = NA, size = 2),
        legend.title = element_text(color="black",size=12.5),
        legend.text = element_text(color="black",size=11.5),
        panel.background = element_blank(),
        plot.title = element_text(size=14.5, hjust = 0.5, face="bold"), 
        axis.title= element_text(color="black",size=12.5, face="bold"), 
        axis.text= element_text(color="black",size=11.5, face="bold"), 
        panel.border = element_rect(colour = "black", fill=NA))+
  ggtitle("Unweighted UniFrac")

tiff("Figures/pcoa.un.UniFrac.tiff", units="in", width = 5, height=4, res=600, compression="lzw")
pcoa_un.UniFrac_plot
dev.off()

tiff("Figures/pcoa.un.UniFrac.2.tiff", units="in", width = 5, height=4, res=600, compression="lzw")
pcoa_un.UniFrac_plot+stat_ellipse(level=0.95)
dev.off()

tiff("Figures/pcoa.un.UniFrac.3.tiff", units="in", width = 5, height=4, res=600, compression="lzw")
pcoa_un.UniFrac_plot+stat_ellipse(level=0.95)+ theme(legend.justification =  "top" )
dev.off()

tiff("Figures/pcoa.un.UniFrac.4.tiff", units="in", width = 5, height=4, res=600, compression="lzw")
pcoa_un.UniFrac_plot+stat_ellipse(level=0.95)+ 
  theme(legend.position =  "bottom")+
  guides(color=guide_legend(nrow=2, byrow=TRUE)) 
dev.off()



####4-week-old Adonis####
  ##4-week-old distance calculation####
SGB.W4 <- read.table("SGB.rel.txt",header=T,row.names = 1) %>% 
  dplyr::select(contains("W4")) %>% t()
SGB.W4[is.na(SGB.W4)] <- 0

dist_bray.W4 <- vegan::vegdist(SGB.W4, method = "bray")
dist_jac.W4 <- vegan::vegdist(SGB.W4, method = "jaccard")

tree <- ape::read.tree("mpa_vOct22_CHOCOPhlAnSGB_202212.nwk")

#construct tree and SGB table for uniFrac analysis
SGB_UniFrac.W4 <- t(SGB.W4) %>% as.data.frame() %>% rownames_to_column(var="SGB")
SGB_UniFrac.W4[,1] <- gsub("t__SGB", "", SGB_UniFrac.W4[,1])
SGB_UniFrac.W4 <- SGB_UniFrac.W4 %>% column_to_rownames(var="SGB")
SGB_removed <- setdiff(rownames(SGB_UniFrac.W4), tree$tip.label)
tree_filt <- ape::keep.tip(tree, setdiff(rownames(SGB_UniFrac.W4), SGB_removed))
SGB_UniFrac.W4 <- SGB_UniFrac.W4[tree_filt$tip.label,]
rm(tree,SGB_removed)

dist_w.UniFrac.W4 <- rbiom::beta.div(as.matrix(SGB_UniFrac.W4), 
                                     tree=tree_filt, method="unifrac", weighted=T)

dist_un.UniFrac.W4 <- rbiom::beta.div(as.matrix(SGB_UniFrac.W4), 
                                      tree=tree_filt, method="unifrac", weighted=F)

  ##4-week-old ADONIS for beta-diversity####
metadata.W4 <- metadata %>% filter(Age=="Weaning")

ADONIS<-tibble()
for (x in c("dist_bray.W4", "dist_jac.W4","dist_w.UniFrac.W4","dist_un.UniFrac.W4")) {
  td <- get(x)
  
  ad <- adonis2(td~Genotype,data=metadata.W4, permutations=999,parallel = 1) %>% 
    as.data.frame() %>% 
    rename(Pvalue=`Pr(>F)`) %>%  
    mutate(dist=paste0(x)) %>% 
    dplyr::select(dist, everything())
  
  ADONIS<-bind_rows(ADONIS, ad)
  
  rm(td,ad,x)
}

####14-week-old Adonis####
  ##14-week-old distance calculation####
SGB.W14 <- read.table("SGB.rel.txt",header=T,row.names = 1) %>% dplyr::select(contains("W14")) %>% t()
SGB.W14[is.na(SGB.W14)] <- 0

dist_bray.W14 <- vegan::vegdist(SGB.W14, method = "bray")
dist_jac.W14 <- vegan::vegdist(SGB.W14, method = "jaccard")

tree <- ape::read.tree("mpa_vOct22_CHOCOPhlAnSGB_202212.nwk")

#construct tree and SGB table for uniFrac analysis
SGB_UniFrac.W14 <- t(SGB.W14) %>% as.data.frame() %>% rownames_to_column(var="SGB")
SGB_UniFrac.W14[,1] <- gsub("t__SGB", "", SGB_UniFrac.W14[,1])
SGB_UniFrac.W14 <- SGB_UniFrac.W14 %>% column_to_rownames(var="SGB")
SGB_removed <- setdiff(rownames(SGB_UniFrac.W14), tree$tip.label)
tree_filt <- ape::keep.tip(tree, setdiff(rownames(SGB_UniFrac.W14), SGB_removed))
SGB_UniFrac.W14 <- SGB_UniFrac.W14[tree_filt$tip.label,]
rm(tree,SGB_removed)

dist_w.UniFrac.W14 <- rbiom::beta.div(as.matrix(SGB_UniFrac.W14), 
                                      tree=tree_filt, method="unifrac", weighted=T)

dist_un.UniFrac.W14 <- rbiom::beta.div(as.matrix(SGB_UniFrac.W14), 
                                       tree=tree_filt, method="unifrac", weighted=F)

  ##14-week-old ADONIS####
metadata.W14 <- metadata %>% filter(Age=="Adult")

ADONIS.W14<-tibble()
for (x in c("dist_bray.W14", "dist_jac.W14","dist_w.UniFrac.W14","dist_un.UniFrac.W14")) {
  td <- get(x)
  
  ad <- adonis2(td~Genotype,data=metadata.W14, permutations=999,parallel = 1) %>% 
    as.data.frame() %>% 
    rename(Pvalue=`Pr(>F)`) %>%  
    mutate(dist=paste0(x)) %>% 
    dplyr::select(dist, everything())
  
  ADONIS.W14<-bind_rows(ADONIS.W14, ad)
  
  rm(td,ad,x)
}

####weighted UniFrac NMDS####
nmds.w.unifrac <- vegan::metaMDS(dist_w.UniFrac, autotransform = F, k=3)

nmds.w.unifrac.df <- data.frame(nmds1 = nmds.w.unifrac$points[,1], 
                                nmds2 = nmds.w.unifrac$points[,2],
                                nmds3 = nmds.w.unifrac$points[,3])

nmds.w.unifrac.df <- nmds.w.unifrac.df %>% 
  rownames_to_column(var="ID") %>% 
  merge(metadata,by="ID") %>% 
  mutate(Genotype_Age=factor(Genotype_Age,
                             level=c("WT_4wk","CF_4wk","WT_14wk","CF_14wk")))

w.unifrac.centroid <- nmds.w.unifrac.df  %>%
group_by(Genotype_Age) %>%
summarize(centroid_nmds1 = mean(nmds1), 
          centroid_nmds2 = mean(nmds2),
          centroid_nmds3 = mean(nmds3))
                                                                              
nmds.w.unifrac.df  <- merge(nmds.w.unifrac.df, w.unifrac.centroid, by = "Genotype_Age", all.x = TRUE)                                                        
                                                                                                                                                            
nmds_w.UniFrac_plot_P1P2 <- ggplot() +
  geom_point(data=nmds.w.unifrac.df, aes(x=nmds1, y=nmds2, color=Genotype_Age), size=4) + theme_minimal() +
  xlab(paste("NMDS1"))+ 
  ylab(paste("NMDS2"))+ 
  scale_color_manual(values=c("#056943","#585858","#3951A2","#A80326"))+
  theme(legend.key = element_rect(fill = NA, colour = NA, size = 2),
        legend.title = element_text(color="black",size=12.5),
        legend.text = element_text(color="black",size=11.5),
        legend.position = "none",
        panel.background = element_blank(),
        plot.title = element_text(size=14.5, hjust = 0.5, face="bold"), 
        axis.title= element_text(color="black",size=12.5, face="bold"), 
        axis.text= element_text(color="black",size=11.5, face="bold"), 
        panel.border = element_rect(colour = "black", fill=NA))+
  ggtitle("Weighted UniFrac")+
  geom_segment(data = nmds.w.unifrac.df, 
               aes(x = nmds1, y = nmds2, 
                   xend = centroid_nmds1, yend = centroid_nmds2, 
                   group = Genotype_Age, color=Genotype_Age), 
               linetype = "dashed")  +
  geom_label (data = nmds.w.unifrac.df, aes(x = centroid_nmds1, 
                                            y = centroid_nmds2, 
                                            label=Genotype_Age,
                                            color=Genotype_Age))


nmds_w.UniFrac_plot_P1P3 <- ggplot() +
  geom_point(data=nmds.w.unifrac.df, aes(x=nmds1, y=nmds3, color=Genotype_Age), size=4) + theme_minimal() +
  xlab(paste("NMDS1"))+ 
  ylab(paste("NMDS3"))+ 
  scale_color_manual(values=c("#056943","#585858","#3951A2","#A80326"))+
  theme(legend.key = element_rect(fill = NA, colour = NA, size = 2),
        legend.title = element_text(color="black",size=12.5),
        legend.text = element_text(color="black",size=11.5),
        legend.position = "none",
        panel.background = element_blank(),
        plot.title = element_text(size=14.5, hjust = 0.5, face="bold"), 
        axis.title= element_text(color="black",size=12.5, face="bold"), 
        axis.text= element_text(color="black",size=11.5, face="bold"), 
        panel.border = element_rect(colour = "black", fill=NA))+
  ggtitle("Weighted UniFrac")+
  geom_segment(data = nmds.w.unifrac.df, 
               aes(x = nmds1, y = nmds3, 
                   xend = centroid_nmds1, yend = centroid_nmds3, 
                   group = Genotype_Age, color=Genotype_Age), 
               linetype = "dashed")  +
  geom_label (data = nmds.w.unifrac.df, 
              aes(x = centroid_nmds1, y = centroid_nmds3, label=Genotype_Age, color=Genotype_Age))


tiff("Figures/nmds_w.UniFrac_plot_P1P3_2.tiff", units="in", width = 4, height=4, res=600, compression="lzw")
nmds_w.UniFrac_plot_P1P3
dev.off()

tiff("Figures/nmds_w.UniFrac_plot_P1P2_2.tiff", units="in", width = 4, height=4, res=600, compression="lzw")
nmds_w.UniFrac_plot_P1P2
dev.off()

plot_ly(x=nmds.w.unifrac.df$nmds1, 
        y=nmds.w.unifrac.df$nmds2, 
        z=nmds.w.unifrac.df$nmds3, 
        color=nmds.w.unifrac.df$Genotype_Age,
        size=1)



plot_ly(x=nmds.un.unifrac.df$nmds1, 
        y=nmds.un.unifrac.df$nmds2, 
        z=nmds.un.unifrac.df$nmds3, 
        color=nmds.un.unifrac.df$Genotype_Age,
        size=1)


####un-weighted NDMS####
nmds.un.unifrac <- vegan::metaMDS(dist_un.UniFrac, autotransform = F, k=3)

nmds.un.unifrac.df <- data.frame(nmds1 = nmds.un.unifrac$points[,1], 
                                nmds2 = nmds.un.unifrac$points[,2],
                                nmds3 = nmds.un.unifrac$points[,3])

nmds.un.unifrac.df <- nmds.un.unifrac.df %>% rownames_to_column(var="ID") %>% 
  merge(metadata,by="ID") %>% mutate(Genotype_Age=factor(Genotype_Age,level=c("WT_4wk","CF_4wk",
                                                                              "WT_14wk","CF_14wk")))

un.unifrac.centroid <- nmds.un.unifrac.df  %>%
  group_by(Genotype_Age) %>%
  summarize(centroid_nmds1 = mean(nmds1), 
            centroid_nmds2 = mean(nmds2),
            centroid_nmds3 = mean(nmds3))

nmds.un.unifrac.df  <- merge(nmds.un.unifrac.df, un.unifrac.centroid, by = "Genotype_Age", all.x = TRUE)  


nmds_un.UniFrac_plot_P1P3 <- ggplot(nmds.un.unifrac.df, 
                              aes(x=nmds1, y=nmds3,color=Genotype_Age)) +
  geom_point(size=4) + theme_minimal() +
  xlab(paste("NMDS1"))+
  ylab(paste("NMDS3"))+
  scale_color_manual(values=c("#056943","#585858","#3951A2","#A80326"))+
  theme(legend.key = element_rect(fill = NA, colour = NA, size = 2),
        legend.title = element_text(color="black",size=12.5),
        legend.text = element_text(color="black",size=11.5),
        legend.position = "none",
        panel.background = element_blank(),
        plot.title = element_text(size=14.5, hjust = 0.5, face="bold"), 
        axis.title= element_text(color="black",size=12.5, face="bold"), 
        axis.text= element_text(color="black",size=11.5, face="bold"), 
        panel.border = element_rect(colour = "black", fill=NA))+
  ggtitle("Unweighted UniFrac")+
  geom_segment(data = nmds.un.unifrac.df, 
               aes(x = nmds1, y = nmds3, 
                   xend = centroid_nmds1, yend = centroid_nmds3, 
                   group = Genotype_Age, color=Genotype_Age), 
               linetype = "dashed")  +
  stat_ellipse(data = nmds.un.unifrac.df, aes(x = nmds1, y = nmds3, color = Genotype_Age), 
               level = 0.95)+
  geom_label (data = nmds.un.unifrac.df, 
              aes(x = centroid_nmds1, y = centroid_nmds3, label=Genotype_Age, color=Genotype_Age))



nmds_un.UniFrac_plot_P1P2 <- ggplot(nmds.un.unifrac.df, 
                                    aes(x=nmds1, y=nmds2,color=Genotype_Age)) +
  geom_point(size=4) + theme_minimal() +
  xlab(paste("NMDS1"))+
  ylab(paste("NMDS2"))+
  scale_color_manual(values=c("#056943","#585858","#3951A2","#A80326"))+
  theme(legend.key = element_rect(fill = NA, colour = NA, size = 2),
        legend.title = element_text(color="black",size=12.5),
        legend.text = element_text(color="black",size=11.5),
        legend.position = "none",
        panel.background = element_blank(),
        plot.title = element_text(size=14.5, hjust = 0.5, face="bold"), 
        axis.title= element_text(color="black",size=12.5, face="bold"), 
        axis.text= element_text(color="black",size=11.5, face="bold"), 
        panel.border = element_rect(colour = "black", fill=NA))+
  ggtitle("Unweighted UniFrac")+
  geom_segment(data = nmds.un.unifrac.df, 
               aes(x = nmds1, y = nmds2, 
                   xend = centroid_nmds1, yend = centroid_nmds2, 
                   group = Genotype_Age, color=Genotype_Age), 
               linetype = "dashed")  + 
  stat_ellipse(data = nmds.un.unifrac.df, aes(x = nmds1, y = nmds2, color = Genotype_Age), 
               level = 0.95)+
  geom_label (data = nmds.un.unifrac.df, 
              aes(x = centroid_nmds1, 
                  y = centroid_nmds2, 
                  label=Genotype_Age, 
                  color=Genotype_Age))



tiff("Figures/nmds_un.UniFrac_plot_P1P2_2.tiff", units="in", width = 4, height=4, res=600, compression="lzw")
nmds_un.UniFrac_plot_P1P2 
dev.off()

tiff("Figures/nmds_un.UniFrac_plot_P1P3_2.tiff", units="in", width = 4, height=4, res=600, compression="lzw")
nmds_un.UniFrac_plot_P1P3
dev.off()

####dispersion analysis#####
disp.w.uni.4wk <- betadisper(dist_w.UniFrac.W4, group =metadata.W4$Genotype)
permutest(disp.w.uni.4wk, pairwise=TRUE, permutations=1000)

disp.un.uni.4wk <- betadisper(dist_un.UniFrac.W4, group =metadata.W4$Genotype)
permutest(disp.un.uni.4wk, pairwise=TRUE, permutations=1000)

disp.w.uni.14wk <- betadisper(dist_w.UniFrac.W14, group =metadata.W14$Genotype)
permutest(disp.w.uni.14wk, pairwise=TRUE, permutations=1000)

disp.un.uni.14wk <- betadisper(dist_un.UniFrac.W14, group =metadata.W14$Genotype)
permutest(disp.un.uni.14wk, pairwise=TRUE, permutations=1000)

####PCoA plot weighted####
pcoa_w.UniFrac_plot_P1P2 <- ggplot(pcoa_w.UniFrac.df, 
                              aes(x=pcoa1, y=pcoa2,color=Genotype_Age)) +
  geom_point(size=4) + theme_minimal() +
  xlab(paste("PCoA1"," (", paste(pcoa_w.UniFrac.pc$PCvar[1]), "%)"))+
  ylab(paste("PCoA2"," (", paste(pcoa_w.UniFrac.pc$PCvar[2]), "%)"))+
  scale_color_manual(values=c("#056943","#585858","#3951A2","#A80326"))+
  theme(legend.key = element_rect(fill = NA, colour = NA, size = 2),
        legend.title = element_text(color="black",size=12.5),
        legend.text = element_text(color="black",size=11.5),
        panel.background = element_blank(),
        plot.title = element_text(size=14.5, hjust = 0.5, face="bold"), 
        axis.title= element_text(color="black",size=12.5, face="bold"), 
        axis.text= element_text(color="black",size=11.5, face="bold"), 
        panel.border = element_rect(colour = "black", fill=NA))+
  ggtitle("Weighted UniFrac")

pcoa_w.UniFrac_plot_P1P3 <- ggplot(pcoa_w.UniFrac.df, 
                                   aes(x=pcoa1, y=pcoa3,color=Genotype_Age)) +
  geom_point(size=4) + theme_minimal() +
  xlab(paste("PCoA1"," (", paste(pcoa_w.UniFrac.pc$PCvar[1]), "%)"))+
  ylab(paste("PCoA3"," (", paste(pcoa_w.UniFrac.pc$PCvar[3]), "%)"))+
  scale_color_manual(values=c("#056943","#585858","#3951A2","#A80326"))+
  theme(legend.key = element_rect(fill = NA, colour = NA, size = 2),
        legend.title = element_text(color="black",size=12.5),
        legend.text = element_text(color="black",size=11.5),
        panel.background = element_blank(),
        plot.title = element_text(size=14.5, hjust = 0.5, face="bold"), 
        axis.title= element_text(color="black",size=12.5, face="bold"), 
        axis.text= element_text(color="black",size=11.5, face="bold"), 
        panel.border = element_rect(colour = "black", fill=NA))+
  ggtitle("Weighted UniFrac")


tiff("Figures/pcoa.w.UniFrac.P1P2.tiff", units="in", width = 5, height=4, res=600, compression="lzw")
pcoa_w.UniFrac_plot_P1P2+stat_ellipse(level=0.95)
dev.off()

tiff("Figures/pcoa.w.UniFrac.P1P3.tiff", units="in", width = 5, height=4, res=600, compression="lzw")
pcoa_w.UniFrac_plot_P1P3+stat_ellipse(level=0.95)
dev.off()

####PCoA plot unweighted#####
pcoa_un.UniFrac_plot_P1P2 <- ggplot(pcoa_un.UniFrac.df, 
                               aes(x=pcoa1, y=pcoa2, color=Genotype_Age)) +
  geom_point(size=4) + theme_minimal() +
  xlab(paste("PCoA1"," (", paste(pcoa_un.UniFrac.pc$PCvar[1]), "%)"))+
  ylab(paste("PCoA2"," (", paste(pcoa_un.UniFrac.pc$PCvar[2]), "%)"))+
  scale_color_manual(values=c("#056943","#585858","#3951A2","#A80326"))+
  theme(legend.key = element_rect(fill = NA, colour = NA, size = 2),
        legend.title = element_text(color="black",size=12.5),
        legend.text = element_text(color="black",size=11.5),
        panel.background = element_blank(),
        plot.title = element_text(size=14.5, hjust = 0.5, face="bold"), 
        axis.title= element_text(color="black",size=12.5, face="bold"), 
        axis.text= element_text(color="black",size=11.5, face="bold"), 
        panel.border = element_rect(colour = "black", fill=NA))+
  ggtitle("Unweighted UniFrac")

pcoa_un.UniFrac_plot_P1P3 <- ggplot(pcoa_un.UniFrac.df, 
                                    aes(x=pcoa1, y=pcoa3, color=Genotype_Age)) +
  geom_point(size=4) + theme_minimal() +
  xlab(paste("PCoA1"," (", paste(pcoa_un.UniFrac.pc$PCvar[1]), "%)"))+
  ylab(paste("PCoA3"," (", paste(pcoa_un.UniFrac.pc$PCvar[3]), "%)"))+
  scale_color_manual(values=c("#056943","#585858","#3951A2","#A80326"))+
  theme(legend.key = element_rect(fill = NA, colour = NA, size = 2),
        legend.title = element_text(color="black",size=12.5),
        legend.text = element_text(color="black",size=11.5),
        panel.background = element_blank(),
        plot.title = element_text(size=14.5, hjust = 0.5, face="bold"), 
        axis.title= element_text(color="black",size=12.5, face="bold"), 
        axis.text= element_text(color="black",size=11.5, face="bold"), 
        panel.border = element_rect(colour = "black", fill=NA))+
  ggtitle("Unweighted UniFrac")

tiff("Figures/pcoa.un.UniFrac_P1P2.tiff", units="in", width = 5, height=4, res=600, compression="lzw")
pcoa_un.UniFrac_plot_P1P2+ stat_ellipse(level=0.95)
dev.off()

tiff("Figures/pcoa.un.UniFrac_P1P3.tiff", units="in", width = 5, height=4, res=600, compression="lzw")
pcoa_un.UniFrac_plot_P1P3+ stat_ellipse(level=0.95)
dev.off()



