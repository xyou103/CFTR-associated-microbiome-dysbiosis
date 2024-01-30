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

## taxa analysis at phylum level####
## phylum relative abundance df####
phyl <- read.table("phylum.rel.txt",header=T)
phyl[is.na(phyl)] <- 0

phyl <- as.data.frame(phyl)%>% rownames_to_column("Phylum") %>%  
  gather(-Phylum, key="ID", value="Abundance") %>% left_join(metadata)

phyl.tp <- phyl %>% group_by(Phylum) %>% summarize(mean=mean(Abundance)) %>%
  arrange(desc(mean))

phyl.order <- phyl %>%
  filter(Phylum=="Firmicutes") %>%
  arrange(desc(Abundance)) %>%
  pull(ID)

phyl.plot.df <- phyl %>%
  mutate(Phylum=factor(Phylum, levels = rev(phyl.tp$Phylum))) %>%
  left_join(metadata) %>% 
  mutate(Phylum=fct_relevel(Phylum, "Bacteria_unclassified")) %>% 
  mutate(ID=factor(ID, levels=phyl.order)) %>% 
  mutate(Genotype_Age=factor(Genotype_Age,level=c("WT_4wk","CF_4wk",
                                                  "WT_14wk","CF_14wk")))

## phylum plot####
phyl.plot <- ggplot(phyl.plot.df, aes(x=ID, y=Abundance, fill=Phylum, width=1)) +
  geom_bar(stat="identity") +
  facet_grid(~ Genotype_Age, scales="free", space = "free")+ylab("Relative abundance (%)")+
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_blank()) +
  theme(axis.text.y = element_text(color="black",size=11),
        axis.text.x=element_blank(),axis.title.x = element_blank(),
        axis.title.y=element_text(size=11), 
        strip.text = element_text(size=11))+
  scale_fill_manual(values=rev(c(
    "firebrick",
    "gold",
    "darkorchid",
    "steelblue2",
    "chartreuse1",
    "darkorange",
    "cadetblue",
    "grey")))

tiff('phyl.plot2.tiff', units="in", width=6, height=4, res=600, compression = 'lzw')
phyl.plot
dev.off()


## phylum U test 4wk ####

phy.4wk <- spread(phyl,key = Phylum,value=Abundance) %>% as.data.frame() %>% 
  dplyr::filter(Age=="Weaning") %>% 
  mutate(Genotype=factor(Genotype, levels=c("CF","WT")))

phy.4wk.res <- data.frame()

for (i in 6:length(phy.4wk)){
  taxa<-colnames(phy.4wk)[i]
  
  mod <- broom::tidy(wilcox.test(phy.4wk[,i] ~ phy.4wk$Genotype, data=phy.4wk, conf.int=TRUE, conf.level=0.95, exact = FALSE))
  
  mod$taxa <- taxa
  
  phy.4wk.res <- rbind(phy.4wk.res, mod)
  rm(mod,taxa)
}

phy.4wk.res$adj.p <- p.adjust(phy.4wk.res$p.value, method ="BH" )

phy.4wk.res <- dplyr::select (phy.4wk.res,taxa,everything())

phy.4wk.res.sig <- phy.4wk.res %>% 
  mutate(significant=case_when(
    adj.p<0.05 & estimate > 0 ~ "increased in CF",
    adj.p<0.05 & estimate < 0 ~ "increased in WT",
    TRUE~"not significant"
  ))%>% filter(adj.p<0.05)

## phylum U test 14wk ####

phy.14wk <- spread(phyl,key = Phylum,value=Abundance) %>% as.data.frame() %>% 
  dplyr::filter(Age=="Adult") %>% 
  mutate(Genotype=factor(Genotype, levels=c("CF","WT")))

phy.14wk.res <- data.frame()

for (i in 6:length(phy.14wk)){
  taxa<-colnames(phy.14wk)[i]
  
  mod <- broom::tidy(wilcox.test(phy.14wk[,i] ~ phy.14wk$Genotype, data=phy.14wk, conf.int=TRUE, conf.level=0.95, exact = FALSE))
  
  mod$taxa <- taxa
  
  phy.14wk.res <- rbind(phy.14wk.res, mod)
  rm(mod,taxa)
}

phy.14wk.res$adj.p <- p.adjust(phy.14wk.res$p.value, method ="BH" )

phy.14wk.res <- dplyr::select (phy.14wk.res,taxa,everything())

phy.14wk.res.sig <- phy.14wk.res %>% 
  mutate(significant=case_when(
    adj.p<0.05 & estimate > 0 ~ "increased in CF",
    adj.p<0.05 & estimate < 0 ~ "increased in WT",
    TRUE~"not significant"
  ))%>% filter(adj.p<0.05)




## taxa analysis at phylum level####
## genus relative abundance df####
genus <- read.table("genus.rel.txt",header=T)
genus[is.na(genus)] <- 0
genus <- as.data.frame(genus)%>% rownames_to_column("Genus") %>%  
  gather(-Genus, key="ID", value="Abundance") %>% left_join(metadata)


#select mean abundance over 0.25% in ctrl and abx group respectively  
genus.tp.4wk.WT <- genus %>% filter(Genotype=="WT") %>% group_by(Genus) %>% summarize(mean=mean(Abundance)) %>%
  arrange(desc(mean)) %>% filter(mean>0.25) %>% dplyr::select(Genus)

genus.tp.4wk.CF <- genus %>% filter(Genotype=="CF") %>% group_by(Genus) %>% summarize(mean=mean(Abundance)) %>%
  arrange(desc(mean)) %>% filter(mean>0.25) %>% dplyr::select(Genus)

genus2phy <- read.table("phy2genus.txt",header=T)

genus2phy$Genus2Phy <- paste(genus2phy$Phylum, ";", genus2phy$Genus)

genus.tp <- rbind(genus.tp.4wk.WT, genus.tp.4wk.CF) %>% unique() %>% filter(!Genus=="Bacteria_unclassified") %>% 
  bind_rows(., tibble(Genus="Bacteria_unclassified")) %>% left_join(genus2phy,"Genus")

genus.tp$Genus2Phy[genus.tp$Genus2Phy=="Bacteria_unclassified ; Bacteria_unclassified"] <- "Bacteria_unclassified"

genus.order <- genus %>%
  filter(Genus=="Ligilactobacillus") %>%
  arrange(desc(Abundance)) %>%
  pull(ID)

genus.plot.df <- genus %>% mutate(Genus=if_else(Genus %in% genus.tp$Genus, Genus, "Bacteria_unclassified")) %>%
  left_join(genus.tp,by = join_by(Genus)) %>% 
  mutate(Genus=factor(Genus2Phy, levels = rev(genus.tp$Genus2Phy))) %>% 
  left_join(metadata) %>% 
  mutate(ID=factor(ID, levels=genus.order)) %>% 
  mutate(Genotype_Age=factor(Genotype_Age,level=c("WT_4wk","CF_4wk","WT_14wk","CF_14wk")))

## genus plot####
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(58)

genus.plot <- ggplot(genus.plot.df, aes(x=ID, y=Abundance, fill=Genus, width=1)) +
  geom_bar(stat="identity") +
  facet_grid(~ Genotype_Age, scales="free", space = "free")+
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_blank()) +ylab("Relative abundance (%)")+
  theme(axis.text.y = element_text(color="black",size=11),
        axis.text.x=element_blank(),axis.title.x = element_blank(),
        axis.title.y=element_text(size=11), 
        strip.text = element_text(size=11),
        legend.text = element_text(size=6.5),
        legend.key.size = unit(0.5, 'cm'),legend.position="right")+
  scale_fill_manual(values=mycolors)

tiff('genus.plot.tiff', units="in", width=12, height=4, res=600, compression = 'lzw')
genus.plot
dev.off()

## genus U test 4wk ####

#import genus table#
genus <- read.table("genus.rel.txt",header=T)
genus[is.na(genus)] <- 0
genus <- as.data.frame(genus)%>% rownames_to_column("Genus") %>%  
  gather(-Genus, key="ID", value="Abundance") %>% left_join(metadata)

#filter 4wk genus data#
genus.4wk <- spread(genus,key = Genus,value=Abundance) %>% as.data.frame() %>% 
  dplyr::filter(Age=="Weaning") %>% 
  mutate(Genotype=factor(Genotype, levels=c("CF","WT"))) %>% dplyr::select(-Clostridioides)

genus.4wk.res <- data.frame()

for (i in 6:length(genus.4wk)){
  taxa<-colnames(genus.4wk)[i]
  
  mod <- broom::tidy(wilcox.test(genus.4wk[,i] ~ genus.4wk$Genotype, data=genus.4wk, conf.int=T, conf.level=0.95, exact = FALSE))
  
  mod$taxa <- taxa
  
  genus.4wk.res <- rbind(genus.4wk.res, mod)
  rm(mod,taxa)
}

genus.4wk.res$adj.p <- p.adjust(genus.4wk.res$p.value, method ="BH" )

genus.4wk.res <- dplyr::select (genus.4wk.res,taxa,everything())

genus.4wk.res.sig <- genus.4wk.res %>% 
  mutate(significant=case_when(
    adj.p<0.05 & estimate > 0 ~ "increased in CF",
    adj.p<0.05 & estimate < 0 ~ "increased in WT",
    TRUE~"not significant"
  ))%>% filter(adj.p<0.05)

## genus U test 14wk ####

genus.14wk <- spread(genus,key = Genus,value=Abundance) %>% as.data.frame() %>% 
  dplyr::filter(Age=="Adult") %>% 
  mutate(Genotype=factor(Genotype, levels=c("CF","WT"))) 

col_sums <- colSums(genus.14wk[6:131])

genus.14wk <- genus.14wk %>% select(-names(col_sums[col_sums == 0]))

genus.14wk.res <- data.frame()

for (i in 6:length(genus.14wk)){
  taxa<-colnames(genus.14wk)[i]
  
  mod <- broom::tidy(wilcox.test(genus.14wk[,i] ~ genus.14wk$Genotype, data=genus.14wk, conf.int=T, conf.level=0.95, exact = FALSE))
  
  mod$taxa <- taxa
  
  genus.14wk.res <- rbind(genus.14wk.res, mod)
  rm(mod,taxa)
}

genus.14wk.res$adj.p <- p.adjust(genus.14wk.res$p.value, method ="BH" )

genus.14wk.res <- dplyr::select (genus.14wk.res,taxa,everything())

genus.14wk.res.sig <- genus.14wk.res %>% 
  mutate(significant=case_when(
    adj.p<0.05 & estimate > 0 ~ "increased in CF",
    adj.p<0.05 & estimate < 0 ~ "increased in WT",
    TRUE~"not significant"
  ))%>% filter(adj.p<0.05)
## select sig genus taxa ####

genus.tp.4wk <- genus.4wk.res.sig %>% dplyr::select(taxa) %>% rename(Genus=taxa)

genus.tp.14wk <- genus.14wk.res.sig %>% dplyr::select(taxa)%>% rename(Genus=taxa)


genus2phy <- read.table("phy2genus.txt",header=T)

genus2phy$Genus2Phy <- paste(genus2phy$Phylum, ";", genus2phy$Genus)

genus.tp <- rbind(genus.tp.4wk, genus.tp.14wk) %>% unique() %>% filter(!Genus=="Bacteria_unclassified") %>% 
  bind_rows(., tibble(Genus="Bacteria_unclassified")) %>% left_join(genus2phy,"Genus")

genus.tp$Genus2Phy[genus.tp$Genus2Phy=="Bacteria_unclassified ; Bacteria_unclassified"] <- "Bacteria_unclassified"

genus.order <- genus %>%
  filter(Genus=="Bifidobacterium") %>%
  arrange(desc(Abundance)) %>%
  pull(ID)

genus.plot.df <- genus %>% mutate(Genus=if_else(Genus %in% genus.tp$Genus, Genus, "Bacteria_unclassified")) %>%
  left_join(genus.tp,by = join_by(Genus)) %>% 
  mutate(Genus=factor(Genus2Phy, levels = rev(genus.tp$Genus2Phy))) %>% 
  left_join(metadata) %>% 
  mutate(ID=factor(ID, levels=genus.order)) %>% 
  mutate(Genotype_Age=factor(Genotype_Age,level=c("WT_4wk","CF_4wk","WT_14wk","CF_14wk")))

mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(69)

genus.plot <- ggplot(genus.plot.df, aes(x=ID, y=Abundance, fill=Genus, width=1)) +
  geom_bar(stat="identity") +
  facet_grid(~ Genotype_Age, scales="free", space = "free")+
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_blank()) +ylab("Relative abundance (%)")+
  theme(axis.text.y = element_text(color="black",size=11),
        axis.text.x=element_blank(),axis.title.x = element_blank(),
        axis.title.y=element_text(size=11), 
        strip.text = element_text(size=11),
        legend.text = element_text(size=6.5),
        legend.key.size = unit(0.5, 'cm'),legend.position="right")+
  scale_fill_manual(values=mycolors)

tiff('genus.plot.sigtaxa.tiff', units="in", width=12, height=4, res=600, compression = 'lzw')
genus.plot
dev.off()

## genus venn plot####
genus.sig <- list()
genus.4wk.res.sig2 <- genus.4wk.res.sig %>% filter(taxa!="Bacteria_unclassified")
genus.14wk.res.sig2 <- genus.14wk.res.sig %>% filter(taxa!="Bacteria_unclassified")

genus.sig$Weaning <- genus.4wk.res.sig2$taxa

genus.sig$Adult <- genus.14wk.res.sig2$taxa

venn.diagram(
  x = genus.sig,
  category.names = c("" , ""),
  filename = 'Genus.venn_diagramm.tiff',
  output=TRUE,
  
  ext.text=F,
  fontfamily = "arial",
  cat.fontfamily = "arial",
  cex=1.5,
  cat.cex=1.5,
  
  imagetype="tiff", units = "in",
  height = 4 , 
  width = 4 , 
  resolution = 600,
  compression = "lzw",
  
  lwd = 2,
  lty = 'blank',
  fill = c("#bdb5e1", "#f9d580")
  
)






## species U test 4wk ####
#import species data#
sps <- read.table("species.rel.txt",header=T)
sps[is.na(sps)] <- 0
sps <- as.data.frame(sps)%>% rownames_to_column("Species") %>%  
  gather(-Species, key="ID", value="Abundance") %>% left_join(metadata)

#filter 4wk sps data#
sps.4wk <- spread(sps,key = Species,value=Abundance) %>% as.data.frame() %>% 
  dplyr::filter(Age=="Weaning") %>% 
  mutate(Genotype=factor(Genotype, levels=c("CF","WT"))) 

col_sums <- colSums(sps.4wk[6:162])

sps.4wk <- sps.4wk %>% select(-names(col_sums[col_sums == 0]))

#U test#
sps.4wk.res <- data.frame()

for (i in 6:length(sps.4wk)){
  taxa<-colnames(sps.4wk)[i]
  
  mod <- broom::tidy(wilcox.test(sps.4wk[,i] ~ sps.4wk$Genotype, data=sps.4wk, conf.int=T, conf.level=0.95, exact = FALSE))
  
  mod$taxa <- taxa
  
  sps.4wk.res <- rbind(sps.4wk.res, mod)
  rm(mod,taxa)
}

sps.4wk.res$adj.p <- p.adjust(sps.4wk.res$p.value, method ="BH" )

sps.4wk.res <- dplyr::select (sps.4wk.res,taxa,everything())

sps.4wk.res.sig <- sps.4wk.res %>% 
  mutate(significant=case_when(
    adj.p<0.05 & estimate > 0 ~ "increased in CF",
    adj.p<0.05 & estimate < 0 ~ "increased in WT",
    TRUE~"not significant"
  ))%>% filter(adj.p<0.05)

## species U test 14wk ####

#filter 14wk sps data#
sps.14wk <- spread(sps,key = Species,value=Abundance) %>% as.data.frame() %>% 
  dplyr::filter(Age=="Adult") %>% 
  mutate(Genotype=factor(Genotype, levels=c("CF","WT"))) 

col_sums <- colSums(sps.14wk[6:162])

sps.14wk <- sps.14wk %>% select(-names(col_sums[col_sums == 0]))

#U test#
sps.14wk.res <- data.frame()

for (i in 6:length(sps.14wk)){
  taxa<-colnames(sps.14wk)[i]
  
  mod <- broom::tidy(wilcox.test(sps.14wk[,i] ~ sps.14wk$Genotype, data=sps.14wk, conf.int=T, conf.level=0.95, exact = FALSE))
  
  mod$taxa <- taxa
  
  sps.14wk.res <- rbind(sps.14wk.res, mod)
  rm(mod,taxa)
}

sps.14wk.res$adj.p <- p.adjust(sps.14wk.res$p.value, method ="BH" )

sps.14wk.res <- dplyr::select (sps.14wk.res,taxa,everything())

sps.14wk.res.sig <- sps.14wk.res %>% 
  mutate(significant=case_when(
    adj.p<0.05 & estimate > 0 ~ "increased in CF",
    adj.p<0.05 & estimate < 0 ~ "increased in WT",
    TRUE~"not significant"
  ))%>% filter(adj.p<0.05)


## species plot for sig taxa####
species.tp.4wk <- sps.4wk.res.sig %>% dplyr::select(taxa) %>% rename(Species=taxa)

species.tp.14wk <- sps.14wk.res.sig %>% dplyr::select(taxa)%>% rename(Species=taxa)

species.tp <- rbind(species.tp.4wk, species.tp.14wk) %>% unique() %>% 
  bind_rows(., tibble(Species="Bacteria_unclassified")) 

species.order <- sps %>%
  filter(Species=="Bifidobacterium_animalis") %>%
  arrange(desc(Abundance)) %>%
  pull(ID)

species.plot.df <- sps %>% 
  mutate(Species=if_else(Species %in% species.tp$Species, Species, "Bacteria_unclassified")) %>%
  left_join(species.tp, by = join_by(Species)) %>% 
  left_join(metadata) %>% 
  mutate(Species=fct_relevel(Species, "Bacteria_unclassified")) %>%
  mutate(ID=factor(ID, levels=species.order)) %>% 
  mutate(Genotype_Age=factor(Genotype_Age,level=c("WT_4wk","CF_4wk","WT_14wk","CF_14wk")))

mycolors <- colorRampPalette(brewer.pal(12, "Set3"))(88)

species.plot <- ggplot(species.plot.df, aes(x=ID, y=Abundance, fill=Species, width=1)) +
  geom_bar(stat="identity") +
  facet_grid(~ Genotype_Age, scales="free", space = "free")+
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_blank()) +ylab("Relative abundance (%)")+
  theme(axis.text.y = element_text(color="black",size=11),
        axis.text.x=element_blank(),axis.title.x = element_blank(),
        axis.title.y=element_text(size=11), 
        strip.text = element_text(size=11),
        legend.text = element_text(size=6.5),
        legend.key.size = unit(0.5, 'cm'),legend.position="right")+
  scale_fill_manual(values=mycolors)

tiff('species.plot.sigtaxa2.tiff', units="in", width=12.5, height=4, res=600, compression = 'lzw')
species.plot
dev.off()

## species venn plot####
sps.sig <- list()
sps.sig$Weaning <- sps.4wk.res.sig$taxa
sps.sig$Adult <- sps.14wk.res.sig$taxa

venn.diagram(
  x = sps.sig,
  category.names = c("" , ""),
  filename = 'Species.venn_diagramm.tiff',
  output=TRUE,
  
  ext.text=F,
  fontfamily = "arial",
  cat.fontfamily = "arial",
  cex=1.5,
  cat.cex=1.5,
  
  imagetype="tiff", units = "in",
  height = 4 , 
  width = 4 , 
  resolution = 600,
  compression = "lzw",
  
  lwd = 2,
  lty = 'blank',
  fill = c("#99b9e9", "#eca680")
)

## common species taxa####
names(sps.4wk.res.sig)[-1] <- paste("4wk.",names(sps.4wk.res.sig)[-1])
names(sps.14wk.res.sig)[-1] <- paste("14wk.",names(sps.14wk.res.sig)[-1])

sps.res.full <- full_join(sps.4wk.res.sig,sps.14wk.res.sig,by="taxa") %>% 
  dplyr::select("taxa","4wk. significant","14wk. significant",everything())

sps.res.inner <- inner_join(sps.4wk.res.sig,sps.14wk.res.sig,by="taxa") %>% 
  dplyr::select("taxa","4wk. significant","14wk. significant",everything())

sps.4wk.uni <- anti_join(sps.4wk.res.sig,sps.14wk.res.sig,by="taxa")
sps.14wk.uni <- anti_join(sps.14wk.res.sig,sps.4wk.res.sig,by="taxa")

## species forest plot#####
sps.4wk.inner <- sps.4wk.res %>% 
  mutate(significant=case_when(
    adj.p<0.05 & estimate > 0 ~ "increased in CF",
    adj.p<0.05 & estimate < 0 ~ "increased in WT",
    TRUE~"not significant"
  ))%>% filter(adj.p<0.05) %>% filter(!is.na(match(taxa, sps.res.inner$taxa)))%>% 
  mutate(group="4wk")

sps.14wk.inner <- sps.14wk.res %>% 
  mutate(significant=case_when(
    adj.p<0.05 & estimate > 0 ~ "increased in CF",
    adj.p<0.05 & estimate < 0 ~ "increased in WT",
    TRUE~"not significant"
  ))%>% filter(adj.p<0.05) %>% 
  filter(!is.na(match(taxa, sps.res.inner$taxa))) %>% 
  mutate(group="14wk")

sps.forest2 <- rbind(sps.4wk.inner,sps.14wk.inner) %>% 
  mutate(group=factor(group,level=c("4wk","14wk")))

####species forest plot ggplot2####

sps.forest.plot2 <-  ggplot(sps.forest2, aes(x=estimate, xmin=conf.low, 
                                             xmax=conf.high, y=taxa, 
                                             col=group, fill=group)) +
  geom_linerange(linewidth=0.3, position=position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, linetype="dashed", color="grey50") +
  geom_point(size=3, stroke = 0, position=position_dodge(width = 0.5)) +
  theme_minimal()  +
  theme(panel.border = element_blank(), axis.line = element_line()) +
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  xlab("mean of difference (95% CI)")+ 
  theme(axis.text = element_text(color="black",size=10),
        axis.title=element_text(size=11), 
        axis.title.y = element_blank(),
        legend.text = element_text(size=10),
        legend.title = element_blank(),
        plot.title = element_text(size=14, hjust = 3, face="bold"),
        legend.position = c(0.85,0.5),
        legend.background = element_rect(colour = 'grey90', fill = 'white', linetype='solid'))+
  theme(panel.grid.minor = element_blank(),
        axis.line.y  = element_blank())+ 
  scale_color_manual(values = c("#4d8ae8", "#e37f64"))+
  ggtitle("Differentially abundant species (4wk & 14wk)")+
  scale_x_continuous(breaks=seq(-10,15,5))

tiff('species.forestplot2.tiff', units="in", width=6.8, height=5, res=600, compression = 'lzw')
sps.forest.plot2
dev.off()

####forest plot filter by 1%####
sps.forest3 <- sps.forest2 %>% filter(abs(estimate)>1) %>% filter(!taxa=="Escherichia_coli")
sps.forest3 <- bind_rows(sps.forest3, sps.forest2[sps.forest2$taxa=="Escherichia_coli",]) %>% 
  arrange(estimate) %>% mutate(taxa=fct_inorder(factor(taxa)))

sps.forest.plot3 <-  ggplot(sps.forest3, aes(x=estimate, xmin=conf.low, 
                                             xmax=conf.high, y=taxa, 
                                             col=group, fill=group)) +
  geom_linerange(linewidth=0.3, position=position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, linetype="dashed", color="grey50") +
  geom_point(size=3, stroke = 0, position=position_dodge(width = 0.5)) +
  theme_minimal()  +
  theme(panel.border = element_blank(), axis.line = element_line()) +
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  xlab("mean of difference (95% CI)")+ 
  theme(axis.text = element_text(color="black",size=10),
        axis.title=element_text(size=11), 
        axis.title.y = element_blank(),
        legend.text = element_text(size=10),
        legend.title = element_blank(),
        plot.title = element_text(size=14, hjust = 0, face="bold"),
        legend.position = c(0.85,0.5),
        legend.background = element_rect(colour = 'grey90', fill = 'white', linetype='solid'))+
  theme(panel.grid.minor = element_blank(),
        axis.line.y  = element_blank())+ 
  scale_color_manual(values = c("#4d8ae8", "#e37f64"))+
  ggtitle("Differentially abundant species (4wk & 14wk)")+
  scale_x_continuous(breaks=seq(-10,15,5))

tiff('species.forestplot3.tiff', units="in", width=6.8, height=2.2, res=600, compression = 'lzw')
sps.forest.plot3
dev.off()

####forest plot CF age-specific species####
sps.4wk.uni2 <- sps.4wk.res %>% 
  mutate(significant=case_when(
    adj.p<0.05 & estimate > 0 ~ "increased in CF",
    adj.p<0.05 & estimate < 0 ~ "increased in WT",
    TRUE~"not significant"
  ))%>% filter(adj.p<0.05) %>% filter(is.na(match(taxa, sps.res.inner$taxa)))%>% 
  mutate(group="4wk")

sps.14wk.uni2 <- sps.14wk.res %>% 
  mutate(significant=case_when(
    adj.p<0.05 & estimate > 0 ~ "increased in CF",
    adj.p<0.05 & estimate < 0 ~ "increased in WT",
    TRUE~"not significant"
  ))%>% filter(adj.p<0.05) %>% 
  filter(is.na(match(taxa, sps.res.inner$taxa))) %>% 
  mutate(group="14wk")

sps.age.uni <- rbind(sps.4wk.uni2,sps.14wk.uni2) %>% 
  mutate(group=factor(group,level=c("4wk","14wk")))


sps.age.forest.plot2 <-  ggplot(sps.age.uni, aes(x=estimate, xmin=conf.low, 
                                                 xmax=conf.high, y=taxa, 
                                                 col=group, fill=group)) +
  geom_linerange(linewidth=0.3, position=position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, linetype="dashed", color="grey50") +
  geom_point(size=3, stroke = 0, position=position_dodge(width = 0.5)) +
  theme_minimal()  +
  theme(panel.border = element_blank(), axis.line = element_line()) +
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  xlab("mean of difference (95% CI)")+ 
  theme(axis.text = element_text(color="black",size=10),
        axis.title=element_text(size=11), 
        axis.title.y = element_blank(),
        legend.text = element_text(size=10),
        legend.title = element_blank(),
        plot.title = element_text(size=14, hjust = 0, face="bold"),
        legend.position = c(0.85,0.5),
        legend.background = element_rect(colour = 'grey90', fill = 'white', linetype='solid'))+
  theme(panel.grid.minor = element_blank(),
        axis.line.y  = element_blank())+ 
  scale_color_manual(values = c("#4d8ae8", "#e37f64"))+
  ggtitle("Age specific differentially abundant species")+
  scale_x_continuous(breaks=seq(-10,15,5))

tiff('species.age.forest.tiff', units="in", width=7.5, height=13, res=600, compression = 'lzw')
sps.age.forest.plot2
dev.off()


####4wk species forest plot filter by 1%####
sps.4wk.uni2 <- sps.4wk.res %>% 
  mutate(significant=case_when(
    adj.p<0.05 & estimate > 0 ~ "increased in CF",
    adj.p<0.05 & estimate < 0 ~ "increased in WT",
    TRUE~"not significant"
  ))%>% filter(adj.p<0.05) %>% filter(is.na(match(taxa, sps.res.inner$taxa)))%>% 
  mutate(group="4wk") %>% filter(abs(estimate) > 1) %>% 
  arrange(estimate) %>% mutate(taxa = fct_inorder(factor(taxa)))

sps.4wk.forest.plot <-  ggplot(sps.4wk.uni2, aes(x=estimate, xmin=conf.low, 
                                                 xmax=conf.high, y=taxa, 
                                                 col=group, fill=group)) +
  geom_linerange(linewidth=0.3, position=position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, linetype="dashed", color="grey50") +
  geom_point(size=3, stroke = 0, position=position_dodge(width = 0.5)) +
  theme_minimal()  +
  theme(panel.border = element_blank(), axis.line = element_line()) +
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  xlab("mean of difference (95% CI)")+ 
  theme(axis.text = element_text(color="black",size=10),
        axis.title=element_text(size=11), 
        axis.title.y = element_blank(),
        legend.text = element_text(size=10),
        legend.title = element_blank(),
        plot.title = element_text(size=14, hjust = 20 , face="bold"),
        legend.position = "none",
        legend.background = element_rect(colour = 'grey90', fill = 'white', linetype='solid'))+
  theme(panel.grid.minor = element_blank(),
        axis.line.y  = element_blank())+ 
  scale_color_manual(values = c("#4d8ae8"))+
  ggtitle("Differentially abundant species (4wk only)")+
  scale_x_continuous(breaks=seq(-10,15,5))

tiff('species.4wk.forest.tiff', units="in", width=6.8, height=5, res=600, compression = 'lzw')
sps.4wk.forest.plot
dev.off()

####14wk species forest plot filter by 1%####
sps.14wk.uni2 <- sps.14wk.res %>% 
  mutate(significant=case_when(
    adj.p<0.05 & estimate > 0 ~ "increased in CF",
    adj.p<0.05 & estimate < 0 ~ "increased in WT",
    TRUE~"not significant"
  ))%>% filter(adj.p<0.05) %>% filter(is.na(match(taxa, sps.res.inner$taxa)))%>% 
  mutate(group="14wk") %>% filter(abs(estimate) > 1) %>% 
  arrange(estimate) %>% mutate(taxa=fct_inorder(factor(taxa)))

sps.14wk.forest.plot <-  ggplot(sps.14wk.uni2, aes(x=estimate, xmin=conf.low, 
                                                   xmax=conf.high, y=taxa, 
                                                   col=group, fill=group)) +
  geom_linerange(linewidth=0.3, position=position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, linetype="dashed", color="grey50") +
  geom_point(size=3, stroke = 0, position=position_dodge(width = 0.5)) +
  theme_minimal()  +
  theme(panel.border = element_blank(), axis.line = element_line()) +
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  xlab("mean of difference (95% CI)")+ 
  theme(axis.text = element_text(color="black",size=10),
        axis.title=element_text(size=11), 
        axis.title.y = element_blank(),
        legend.text = element_text(size=10),
        legend.title = element_blank(),
        plot.title = element_text(size=14, hjust = 0 , face="bold"),
        legend.position = "none",
        legend.background = element_rect(colour = 'grey90', fill = 'white', linetype='solid'))+
  theme(panel.grid.minor = element_blank(),
        axis.line.y  = element_blank())+ 
  scale_color_manual(values = c("#e37f64"))+
  ggtitle("Differentially abundant species (14wk only)")+
  scale_x_continuous(breaks=seq(-10,15,5))

tiff('species.14wk.forest.tiff', units="in", width=6.8, height=2, res=600, compression = 'lzw')
sps.14wk.forest.plot
dev.off()



