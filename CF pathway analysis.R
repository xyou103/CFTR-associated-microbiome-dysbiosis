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

####import functional pathway####
path <- read.table("path_cpm_unstratified.txt",sep="\t",quote = "", header=T)

unclassified <- path[1,2:length(path)]+path[2,2:length(path)]
unclassified$path_ID <- c("Unclassified")
path <- rbind(unclassified,path) %>% dplyr::select(path_ID,everything())
path <- path %>% filter(!path_ID=="UNMAPPED") %>% filter(!path_ID=="UNINTEGRATED")

path_des <- read.csv("path_cpm_unstrat_des.csv",header=T)

path_df <- path %>%  
  column_to_rownames(var="path_ID") %>% t() %>% as.data.frame() %>% 
  rownames_to_column("ID") %>% merge(metadata,by="ID") %>% 
  dplyr::select("Sample_ID","Genotype","Age","Genotype_Age",everything()) %>% 
  column_to_rownames("ID")

rm(unclassified)

####PCA plot####
path.pca <- prcomp(path_df[,6:length(path_df)],scale=T)
summary(path.pca) 

percentage <- round(path.pca$sdev^2/sum(path.pca$sdev^2)*100,2)

path.pca.df <- as.data.frame(path.pca$x)

percentage <- paste(colnames(path.pca.df),
                    "(",paste(as.character(percentage),
                              "%",")"))

path.pca.df <- path.pca.df %>% 
  rownames_to_column("ID")

path.pca.x <- merge(path.pca.df,metadata,by="ID") %>% 
  mutate(Genotype_Age=factor(Genotype_Age,level=c("WT_4wk","CF_4wk",
                                                  "WT_14wk","CF_14wk")))

####export PCA plot####
path.pca.plot <- ggplot(path.pca.x,aes(x=PC1,y=PC2,color=Genotype_Age))+ 
  geom_point(size=3)+ 
  xlab(percentage[1]) + 
  ylab(percentage[2])+
  theme_bw()+scale_color_manual(values=c("#056943",
                                         "#585858",
                                         "#3951A2",
                                         "#A80326"))+
  theme_minimal() +
  theme(legend.key = element_rect(fill = NA, colour = NA, size = 2),
        legend.title = element_text(color="black",size=12.5),
        legend.text = element_text(color="black",size=11.5),
        legend.position = "none",
        panel.background = element_blank(),
        plot.title = element_text(size=14.5, hjust = 0.5, face="bold"), 
        axis.title= element_text(color="black",size=12.5, face="bold"), 
        axis.text= element_text(color="black",size=11.5, face="bold"), 
        panel.border = element_rect(colour = "black", fill=NA))+
  ggtitle("PCA_Pathway")


tiff("Figures/path.pca.tiff", units="in", width = 4, height=4, res=600, compression="lzw")
path.pca.plot+ stat_ellipse(geom="polygon", 
                            alpha = 0, 
                            show.legend = F, 
                            level = 0.95)
dev.off()

####path adonis####
path.4wk.dist.df <- path %>% filter(!path_ID=="Unclassified") %>% dplyr::select(!contains("W14")) %>% column_to_rownames(var="path_ID") %>% t()

path.dist.4wk <- vegan::vegdist(path.4wk.dist.df, method = "euclidean")

path.4wk.adnois <- adonis2(path.dist.4wk~Genotype,data=metadata.W4, permutations=999,parallel = 1) %>% 
  as.data.frame()


path.14wk.dist.df <- path %>% filter(!path_ID=="Unclassified") %>% dplyr::select(!contains("W4")) %>% column_to_rownames(var="path_ID") %>% t()

path.dist.14wk <- vegan::vegdist(path.14wk.dist.df, method = "euclidean")

path.14wk.adnois <- adonis2(path.dist.14wk~Genotype,data=metadata.W14, permutations=999,parallel = 1) %>% 
  as.data.frame()

#WT vs CF 4wk U test #####
path.4wk <- path_df %>% 
  dplyr::filter(Age=="Weaning") %>% 
  mutate(Genotype=factor(Genotype, levels=c("CF","WT")))

#filter out less than 100 cpm per sample
path.4wk <- cbind(path.4wk[,1:4],
                  path.4wk[6:length(path.4wk)]
                  [,colSums(path.4wk[6:length(path.4wk)])>100])

path.4wk.res <- data.frame()

for (i in 5:length(path.4wk)){
  path_ID<-colnames(path.4wk)[i]
  
  mod <- broom::tidy(wilcox.test(path.4wk[,i] ~ path.4wk$Genotype, data=path.4wk, 
                                 conf.int=TRUE, conf.level=0.95, exact = FALSE))
  
  mod$path_ID <- path_ID
  
  path.4wk.res <- rbind(path.4wk.res, mod)
  rm(mod,path_ID)
}

rm(i)

path.4wk.res$adj.p <- p.adjust(path.4wk.res$p.value, method ="BH" )

path.4wk.res <- dplyr::select (path.4wk.res,path_ID,everything()) %>% 
  merge(path_des,by="path_ID")

path.4wk.res.sig <- path.4wk.res %>% 
  mutate(significant=case_when(
    adj.p<0.05 & estimate > 0 ~ "increased in CF",
    adj.p<0.05 & estimate < 0 ~ "increased in WT",
    TRUE~"not significant"
  ))%>% filter(adj.p<0.05)

#WT vs CF 14wk U test####
path.14wk <- path_df %>% 
  dplyr::filter(Age=="Adult") %>% 
  mutate(Genotype=factor(Genotype, levels=c("CF","WT")))

#filter out less than 100 cpm per sample
path.14wk <- cbind(path.14wk[,1:4],
                   path.14wk[6:length(path.14wk)]
                   [,colSums(path.14wk[6:length(path.14wk)])>100])

path.14wk.res <- data.frame()

for (i in 5:length(path.14wk)){
  path_ID<-colnames(path.14wk)[i]
  
  mod <- broom::tidy(wilcox.test(path.14wk[,i] ~ path.14wk$Genotype, data=path.14wk, 
                                 conf.int=TRUE, conf.level=0.95, exact = FALSE))
  
  mod$path_ID <- path_ID
  
  path.14wk.res <- rbind(path.14wk.res, mod)
  rm(mod,path_ID)
}

rm(i)

path.14wk.res$adj.p <- p.adjust(path.14wk.res$p.value, method ="BH" )

path.14wk.res <- dplyr::select (path.14wk.res,path_ID,everything()) %>% 
  merge(path_des,by="path_ID")

path.14wk.res.sig <- path.14wk.res %>% 
  mutate(significant=case_when(
    adj.p<0.05 & estimate > 0 ~ "increased in CF",
    adj.p<0.05 & estimate < 0 ~ "increased in WT",
    TRUE~"not significant"
  ))%>% filter(adj.p<0.05)



####path venn plot####
path.sig <- list()
path.sig$Weaning <- path.4wk.res.sig$path_ID
path.sig$Adult <- path.14wk.res.sig$path_ID

venn.diagram(
  x = list(A=path.sig$Weaning,B=path.sig$Adult),
  category.names = c("" , ""),
  filename = 'Path.venn_diagramm.tiff',
  output=TRUE,
  
  ext.text=F,
  fontfamily = "arial",
  cat.fontfamily = "arial",
  cex=1.5,
  cat.cex=1.5,
  
  #output feature
  imagetype="tiff", units = "in",
  height = 4 , 
  width = 4 , 
  resolution = 600,
  compression = "lzw",
  
  lwd = 2,
  lty = 'blank',
  fill = c("#8EC58B", "#F8C471"),
  rotation.degree = 180
  
)

####CHObio####
path.CHObio.df <- path.forest.df %>% filter(str_detect(superclass, "Carbohydrate Biosynthesis")) %>% 
  mutate(group2=factor(group2,levels=c("14wk only","4wk only","4wk & 14wk"))) %>% 
  arrange(group2) %>%
  arrange(desc(superclass)) %>% 
  mutate(path_ID = fct_inorder(factor(path_ID)))

path.CHObio.plot <-  ggplot(path.CHObio.df, aes(x=estimate, xmin=conf.low, 
                                                xmax=conf.high, y=path_ID, 
                                                col=group, fill=group)) + 
  geom_linerange(linewidth=0.3, position=position_dodge(width = 0.6)) +
  geom_vline(xintercept = 0, linetype="dashed", color="grey50") +
  geom_point(size=3, stroke = 0, position=position_dodge(width = 0.6)) +
  theme_minimal()  +
  theme(panel.border = element_blank(), axis.line = element_line()) +
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  xlab("mean of difference (95% CI)")+ 
  theme(axis.text = element_text(color="black",size=10),
        axis.title=element_text(size=11), 
        axis.title.y = element_blank(),
        legend.text = element_text(size=10),
        legend.title = element_blank(),
        plot.title = element_text(size=14, hjust = 0.5, face="bold"),
        legend.position = "none",
        legend.background = element_rect(colour = 'grey90', fill = 'white', linetype='solid'),legend.key = element_blank())+
  theme(panel.grid.minor = element_blank(),
        axis.line.y  = element_blank())+ 
  scale_color_manual(values = c("#558552", "#f7b54a"))+
  ggtitle("Carbohydrate Biosynthesis")+
  scale_x_continuous(breaks=seq(0,320,40))

tiff('path.CHObio.tiff', units="in", width=6.8, height=7, res=600, compression = 'lzw')
path.CHObio.plot
dev.off()

####CHOdeg####
path.CHOdeg.df <- path.forest.df %>% filter(str_detect(superclass, "Carbohydrate Degradation")) %>% 
  mutate(group2=factor(group2,levels=c("14wk only","4wk only","4wk & 14wk"))) %>% 
  arrange(group2) %>%
  arrange(desc(superclass)) %>% 
  mutate(path_ID = fct_inorder(factor(path_ID)))

path.CHOdeg.plot <-  ggplot(path.CHOdeg.df, aes(x=estimate, xmin=conf.low, 
                                                xmax=conf.high, y=path_ID, 
                                                col=group, fill=group)) + 
  geom_linerange(linewidth=0.3, position=position_dodge(width = 0.6)) +
  geom_vline(xintercept = 0, linetype="dashed", color="grey50") +
  geom_point(size=3, stroke = 0, position=position_dodge(width = 0.6)) +
  theme_minimal()  +
  theme(panel.border = element_blank(), axis.line = element_line()) +
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  xlab("mean of difference (95% CI)")+ 
  theme(axis.text = element_text(color="black",size=10),
        axis.title=element_text(size=11), 
        axis.title.y = element_blank(),
        legend.text = element_text(size=10),
        legend.title = element_blank(),
        plot.title = element_text(size=14, hjust = 0.5, face="bold"),
        legend.position = "none",
        legend.background = element_rect(colour = 'grey90', fill = 'white', linetype='solid'),legend.key = element_blank())+
  theme(panel.grid.minor = element_blank(),
        axis.line.y  = element_blank())+ 
  scale_color_manual(values = c("#558552", "#f7b54a"))+
  ggtitle("Carbohydrate Degradation")+
  scale_x_continuous(breaks=seq(0,320,40))

tiff('path.CHOdeg.tiff', units="in", width=6.8, height=7, res=600, compression = 'lzw')
path.CHOdeg.plot
dev.off()

####FAbio####
path.FAbio.df <- path.forest.df %>% filter(str_detect(superclass, "Fatty Acid and Lipid Biosynthesis")) %>% 
  mutate(group2=factor(group2,levels=c("14wk only","4wk only","4wk & 14wk"))) %>% 
  arrange(estimate) %>% 
  arrange(group2) %>%
  arrange(desc(superclass)) %>% 
  mutate(path_ID = fct_inorder(factor(path_ID)))

path.FAbio.plot <-  ggplot(path.FAbio.df, aes(x=estimate, xmin=conf.low, 
                                              xmax=conf.high, y=path_ID, 
                                              col=group, fill=group)) + 
  geom_linerange(linewidth=0.3, position=position_dodge(width = 0.6)) +
  geom_vline(xintercept = 0, linetype="dashed", color="grey50") +
  geom_point(size=3, stroke = 0, position=position_dodge(width = 0.6)) +
  theme_minimal()  +
  theme(panel.border = element_blank(), axis.line = element_line()) +
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  xlab("mean of difference (95% CI)")+ 
  theme(axis.text = element_text(color="black",size=10),
        axis.title=element_text(size=11), 
        axis.title.y = element_blank(),
        legend.text = element_text(size=10),
        legend.title = element_blank(),
        plot.title = element_text(size=14, hjust = 0.5, face="bold"),
        legend.position = "none",
        legend.background = element_rect(colour = 'grey90', fill = 'white', linetype='solid'),legend.key = element_blank())+
  theme(panel.grid.minor = element_blank(),
        axis.line.y  = element_blank())+ 
  scale_color_manual(values = c("#558552", "#f7b54a"))+
  ggtitle("Fatty Acid and Lipid Biosynthesis")+
  scale_x_continuous(breaks=seq(0,320,40))

tiff('path.FAbio.tiff', units="in", width=6.8, height=5.5, res=600, compression = 'lzw')
path.FAbio.plot
dev.off()

####NuBio####
path.Nubio.df <- path.forest.df %>% filter(str_detect(superclass, "Nucleoside and Nucleotide Biosynthesis")) %>% 
  mutate(group2=factor(group2,levels=c("14wk only","4wk only","4wk & 14wk"))) %>% 
  arrange(group2) %>%
  arrange(desc(superclass)) %>% 
  mutate(path_ID = fct_inorder(factor(path_ID)))

path.Nubio.plot <-  ggplot(path.Nubio.df, aes(x=estimate, xmin=conf.low, 
                                              xmax=conf.high, y=path_ID, 
                                              col=group, fill=group)) + 
  geom_linerange(linewidth=0.3, position=position_dodge(width = 0.6)) +
  geom_vline(xintercept = 0, linetype="dashed", color="grey50") +
  geom_point(size=3, stroke = 0, position=position_dodge(width = 0.6)) +
  theme_minimal()  +
  theme(panel.border = element_blank(), axis.line = element_line()) +
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  xlab("mean of difference (95% CI)")+ 
  theme(axis.text = element_text(color="black",size=10),
        axis.title=element_text(size=11), 
        axis.title.y = element_blank(),
        legend.text = element_text(size=10),
        legend.title = element_blank(),
        plot.title = element_text(size=14, hjust = 0.5, face="bold"),
        legend.position = "none",
        legend.background = element_rect(colour = 'grey90', fill = 'white', linetype='solid'),legend.key = element_blank())+
  theme(panel.grid.minor = element_blank(),
        axis.line.y  = element_blank())+ 
  scale_color_manual(values = c("#558552", "#f7b54a"))+
  ggtitle("Nucleoside and Nucleotide Biosynthesis")+
  scale_x_continuous(breaks=seq(0,320,40))

tiff('path.Nubio.tiff', units="in", width=6.8, height=7, res=600, compression = 'lzw')
path.Nubio.plot
dev.off()

####AAbio####
path.AAbio.df <- path.forest.df %>% filter(str_detect(superclass, "Amino Acid Biosynthesis")) %>% 
  mutate(group2=factor(group2,levels=c("14wk only","4wk only","4wk & 14wk"))) %>% 
  arrange(group2) %>%
  arrange(desc(superclass)) %>% 
  mutate(path_ID = fct_inorder(factor(path_ID)))

path.AAbio.plot <-  ggplot(path.AAbio.df, aes(x=estimate, xmin=conf.low, 
                                              xmax=conf.high, y=path_ID, 
                                              col=group, fill=group)) + 
  geom_linerange(linewidth=0.3, position=position_dodge(width = 0.6)) +
  geom_vline(xintercept = 0, linetype="dashed", color="grey50") +
  geom_point(size=3, stroke = 0, position=position_dodge(width = 0.6)) +
  theme_minimal()  +
  theme(panel.border = element_blank(), axis.line = element_line()) +
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  xlab("mean of difference (95% CI)")+ 
  theme(axis.text = element_text(color="black",size=10),
        axis.title=element_text(size=11), 
        axis.title.y = element_blank(),
        legend.text = element_text(size=10),
        legend.title = element_blank(),
        plot.title = element_text(size=14, hjust = 0.5, face="bold"),
        legend.position = "none",
        legend.background = element_rect(colour = 'grey90', fill = 'white', linetype='solid'),legend.key = element_blank())+
  theme(panel.grid.minor = element_blank(),
        axis.line.y  = element_blank())+ 
  scale_color_manual(values = c("#558552", "#f7b54a"))+
  ggtitle("Amino Acid Biosynthesis")+
  scale_x_continuous(breaks=seq(0,320,40))

tiff('path.AAbio.tiff', units="in", width=6.8, height=7, res=600, compression = 'lzw')
path.AAbio.plot
dev.off()

####inflammation####
path.inflammation.df <- path.forest.df %>% filter(superclass2=="Inflammation") %>% 
  mutate(group2=factor(group2,levels=c("14wk only","4wk only","4wk & 14wk"))) %>% 
  arrange(estimate) %>% arrange(group2) %>% mutate(path_ID = fct_inorder(factor(path_ID)))

path.inflammation.plot <-  ggplot(path.inflammation.df, aes(x=estimate, xmin=conf.low, 
                                                            xmax=conf.high, y=path_ID, 
                                                            col=group, fill=group)) + 
  geom_linerange(linewidth=0.3, position=position_dodge(width = 0.6)) +
  geom_vline(xintercept = 0, linetype="dashed", color="grey50") +
  geom_point(size=3, stroke = 0, position=position_dodge(width = 0.6)) +
  theme_minimal()  +
  theme(panel.border = element_blank(), axis.line = element_line()) +
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  xlab("mean of difference (95% CI)")+ 
  theme(axis.text = element_text(color="black",size=10),
        axis.title=element_text(size=11), 
        axis.title.y = element_blank(),
        legend.text = element_text(size=10),
        legend.title = element_blank(),
        plot.title = element_text(size=14, hjust = 0.5, face="bold"),
        legend.position = "none",
        legend.background = element_rect(colour = 'grey90', fill = 'white', linetype='solid'),legend.key = element_blank())+
  theme(panel.grid.minor = element_blank(),
        axis.line.y  = element_blank())+ 
  scale_color_manual(values = c("#558552", "#f7b54a"))+
  ggtitle("LPS Biosynthesis")+
  scale_x_continuous(breaks=seq(0,320,40)) 

tiff('path.inflammation.tiff', units="in", width=6.8, height=2, res=600, compression = 'lzw')
path.inflammation.plot
dev.off()

####fermentation####
path.fermentation.df <- path.forest.df %>% filter(superclass2=="Fermentation") %>% 
  mutate(group2=factor(group2,levels=c("14wk only","4wk only","4wk & 14wk"))) %>% 
  arrange(group2) %>%  arrange(desc(superclass)) %>% mutate(path_ID = fct_inorder(factor(path_ID)))

path.fermentation.plot <-  ggplot(path.fermentation.df, aes(x=estimate, xmin=conf.low, 
                                                            xmax=conf.high, y=path_ID, 
                                                            col=group, fill=group)) + 
  geom_linerange(linewidth=0.3, position=position_dodge(width = 0.6)) +
  geom_vline(xintercept = 0, linetype="dashed", color="grey50") +
  geom_point(size=3, stroke = 0, position=position_dodge(width = 0.6)) +
  theme_minimal()  +
  theme(panel.border = element_blank(), axis.line = element_line()) +
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  xlab("mean of difference (95% CI)")+ 
  theme(axis.text = element_text(color="black",size=10),
        axis.title=element_text(size=11), 
        axis.title.y = element_blank(),
        legend.text = element_text(size=10),
        legend.title = element_blank(),
        plot.title = element_text(size=14, hjust = 0.5, face="bold"),
        legend.position = "none",
        legend.background = element_rect(colour = 'grey90', fill = 'white', linetype='solid'),legend.key = element_blank())+
  theme(panel.grid.minor = element_blank(),
        axis.line.y  = element_blank())+ 
  scale_color_manual(values = c("#558552", "#f7b54a"))+
  ggtitle("SCFA Metabolism")+
  scale_x_continuous(breaks=seq(0,320,40))

tiff('path.SCFA Metabolism.tiff', units="in", width=6.8, height=5.2, res=600, compression = 'lzw')
path.fermentation.plot
dev.off()







