#environment 
library(tidyverse)
library(ggplot2)
library(lmerTest)
library(broom.mixed)
library(emmeans)

#data organization####
#import metadata
metadata <- read.csv("metadata.csv",header=T)

#spf weight data import
spf_weight_df <- read.csv("SPF_weight.csv",header=T)
spf_weight_df <- spf_weight_df %>% gather(mouse,weight,-Age)
spf_weight_df <- merge(spf_weight_df, metadata, by="mouse") %>% 
  mutate(genotype=relevel(as.factor(genotype),ref = "WT")) 

#gf weight data import
gf_weight_df <- read.csv("gf_weight.csv",header=T)
gf_weight_df <- gf_weight_df %>% gather(mouse,weight,-Age)
gf_weight_df <- merge(gf_weight_df, metadata, by="mouse")%>% 
  mutate(genotype=relevel(as.factor(genotype),ref = "WT")) 

#spf length data import
spf_length_df <- read.csv("SPF_length.csv",header=T)
spf_length_df <- spf_length_df %>% gather(mouse,length,-Age)
spf_length_df <- merge(spf_length_df, metadata, by="mouse")%>% 
  mutate(genotype=relevel(as.factor(genotype),ref = "WT")) 

#gf length data import
gf_length_df <- read.csv("GF_length.csv",header=T)
gf_length_df <- gf_length_df %>% gather(mouse,length,-Age)
gf_length_df <- merge(gf_length_df, metadata, by="mouse")%>% 
  mutate(genotype=relevel(as.factor(genotype),ref = "WT")) 

#for loop to get the slope and 95% CI for male length ####
df <- rbind(gf_length_df,spf_length_df)
df <- df %>% rename(age=Age) %>% filter(sex=="M")
breakpoints <- c(28, 56)

results <- data.frame(genotype = character(), 
                      condition = character(), 
                      segment = integer(), 
                      slope = numeric(),
                      conf.low = numeric(), 
                      conf.high = numeric(),
                      stringsAsFactors = FALSE)

for (g in unique(df$genotype)) {
  for (c in unique(df$condition[df$genotype == g])) {
    subset_data <- df[df$genotype == g & df$condition == c, ]
    for (i in 1:(length(breakpoints)+1)) {
      if (i == 1) {
        min_age <- min(subset_data$age)
        max_age <- breakpoints[i]
      } else if (i == length(breakpoints)+1) {
        min_age <- breakpoints[i-1]
        max_age <- max(subset_data$age)
      } else {
        min_age <- breakpoints[i-1]
        max_age <- breakpoints[i]
      }
      segment_data <- subset_data[subset_data$age >= min_age & subset_data$age <= max_age, ]
      model <- lm(length ~ age, data = segment_data)
      results <- rbind(results, data.frame(genotype = g, condition = c, age_segment = i, 
                                           slope = coef(model)[2],
                                           conf.low = confint(model)[2,1], 
                                           conf.high = confint(model)[2,2], 
                                           stringsAsFactors = FALSE))
    }
  }
}

male.res.length.slope <- results

rm(c,g,i,max_age,min_age,subset_data,segment_data,model,results)

#DID between SPF vs GF using lmer for each age_segment for male length #####
#data organization
df <- rbind(gf_length_df,spf_length_df)
df <- df %>% rename(age=Age) %>% filter(sex=="M")
df$age_segment <- cut(df$age, breaks = c(-1, 28, 56, 105),
                      labels = c("age1", "age2", "age3"))
#model construction
model_list <- list()
for (i in 1:3) {
  df_i <- df[df$age_segment == paste0("age", i),]
  model_i <- lmer(length ~ age + genotype * condition + genotype * age + (1|mouse), data = df_i)
  
  model_list[[i]] <- model_i
}

#result summary
DID_combined <- map_dfr(model_list, ~{
  tidy_model <- tidy(.x, conf.int = TRUE)
  
  i <- which(sapply(model_list, identical, .x))
  
  tidy_model %>% mutate(age_segment = paste0("age", i))
})

male_DID_length <- DID_combined[, c("age_segment", "term", "estimate", "std.error", "p.value", "conf.low", "conf.high")] %>% 
  mutate(significance = case_when(
    0.01 < p.value & p.value < 0.05 ~ "*",
    0.001 < p.value & p.value <= 0.01 ~ "**",
    p.value <= 0.001 ~ "***",
    TRUE ~ ""
  ))

rm(i,model_i, model_list,df_i,DID_combined)

#for loop to get the slope and 95% CI for male weight ####
df <- rbind(gf_weight_df,spf_weight_df)
df <- df %>% rename(age=Age) %>% filter(sex=="M")
breakpoints <- c(28, 56)

results <- data.frame(genotype = character(), 
                      condition = character(), 
                      segment = integer(), 
                      slope = numeric(),
                      conf.low = numeric(), 
                      conf.high = numeric(),
                      stringsAsFactors = FALSE)

for (g in unique(df$genotype)) {
  for (c in unique(df$condition[df$genotype == g])) {
    subset_data <- df[df$genotype == g & df$condition == c, ]
    for (i in 1:(length(breakpoints)+1)) {
      if (i == 1) {
        min_age <- min(subset_data$age)
        max_age <- breakpoints[i]
      } else if (i == length(breakpoints)+1) {
        min_age <- breakpoints[i-1]
        max_age <- max(subset_data$age)
      } else {
        min_age <- breakpoints[i-1]
        max_age <- breakpoints[i]
      }
      segment_data <- subset_data[subset_data$age >= min_age & subset_data$age <= max_age, ]
      model <- lm(weight ~ age, data = segment_data)
      results <- rbind(results, data.frame(genotype = g, condition = c, age_segment = i, 
                                           slope = coef(model)[2],
                                           conf.low = confint(model)[2,1], 
                                           conf.high = confint(model)[2,2], 
                                           stringsAsFactors = FALSE))
    }
  }
}

male.res.weight.slope <- results

rm(c,g,i,max_age,min_age,subset_data,segment_data,model,results)

#DID between SPF vs GF using lmer for each age_segment for male weight #####

df <- rbind(gf_weight_df,spf_weight_df)
df <- df %>% rename(age=Age) %>% filter(sex=="M")

df$age_segment <- cut(df$age, breaks = c(-1, 28, 56, 105), 
                      labels = c("age1", "age2", "age3"))

model_list <- list()

for (i in 1:3) {
  df_i <- df[df$age_segment == paste0("age", i),]
  model_i <- lmer(weight ~ age + genotype * condition + genotype * age + (1|mouse), data = df_i)
  
  model_list[[i]] <- model_i
}

DID_combined <- map_dfr(model_list, ~{
  tidy_model <- tidy(.x, conf.int = TRUE)
  
  i <- which(sapply(model_list, identical, .x))
  
  tidy_model %>% mutate(age_segment = paste0("age", i))
})

male_DID_weight <- DID_combined[, c("age_segment", "term", "estimate", "std.error", "p.value", "conf.low", "conf.high")] %>% 
  mutate(significance = case_when(
    0.01 < p.value & p.value < 0.05 ~ "*",
    0.001 < p.value & p.value <= 0.01 ~ "**",
    p.value <= 0.001 ~ "***",
    TRUE ~ ""
  ))

rm(i,model_i, model_list,df_i,DID_combined)

#export table of slope and DID####
#DID
male_DID_length <- male_DID_length %>% mutate(sex="M") %>% mutate (para="length")
male_DID_weight <- male_DID_weight %>% mutate(sex="M") %>% mutate (para="weight")

DID_res <- bind_rows(male_DID_length, male_DID_weight)
write.csv(DID_res,"DID_res.csv")

#slope
male.res.length.slope <- male.res.length.slope %>% mutate(sex="M") %>% mutate (para="length")
male.res.weight.slope <- male.res.weight.slope %>% mutate(sex="M") %>% mutate (para="weight")

Slope_res <- bind_rows(male.res.length.slope, male.res.weight.slope)
write.csv(Slope_res,"Slope_res.csv")

#male length plot####
#data organization
df <- rbind(gf_length_df,spf_length_df) %>% rename(age=Age) %>% 
  filter(sex=="M") %>% 
  mutate(condition=relevel(as.factor(condition),ref = "SPF"))

df1 <- df %>% filter(age<=28)
df2 <- df %>% filter(age>=28 & age <=56)
df3 <- df %>% filter (age>=56)


#plot
male_length_plot <- ggplot(df, aes(x = age, y = length, color=genotype,shape=genotype)) +
  geom_point(size=1.2,shape=1) +
  geom_smooth(data=df1,aes(x = age, y = length), method="lm", se=F, linetype = "dashed") + 
  geom_smooth(data=df2,aes(x = age, y = length), method="lm", se=F, linetype = "dashed") + 
  geom_smooth(data=df3,aes(x = age, y = length), method="lm", se=F,linetype = "dashed") +
  labs(x = "Age (day)", y = "Body length (mm)") + 
  scale_colour_manual(values=c("royalblue4", "firebrick"))+
  theme(text=element_text(size=14),
        axis.title= element_text(color="black",size=15),
        axis.text= element_text(color="black",size=14), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        legend.position="right", 
        plot.title = element_text(size=16, hjust = 0.5),
        panel.border = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        legend.key = element_rect(colour = NA, fill = NA),
        strip.background = element_rect(fill="white", color = "white"))+
  scale_x_continuous(limits = c(0, 105), breaks = c(0, 14, 28, 42, 56, 70,84,98))+
  scale_y_continuous(limits = c(25, 105), breaks = c(25, 50, 75, 100))+
  facet_wrap(~ condition, scales = "free")+
  ggtitle("Male body length")+
  guides(shape = F)

tiff('length_plot_male.tiff', units="in", width=8, height=6, res=300, compression = 'lzw')
male_length_plot
dev.off()

rm(df,df1,df2,df3)

#male weight plot####
df <- rbind(gf_weight_df,spf_weight_df) %>% rename(age=Age) %>% 
  filter(sex=="M") %>% 
  mutate(condition=relevel(as.factor(condition),ref = "SPF"))

df1 <- df %>% filter(age<=28)
df2 <- df %>% filter(age>=28 & age <=56)
df3 <- df %>% filter (age>=56)


#plot
male_weight_plot <- ggplot(df, aes(x = age, y = weight, color=genotype,shape=genotype)) +
  geom_point(size=1.2,shape=1) +
  geom_smooth(data=df1,aes(x = age, y = weight), method="lm", se=F, linetype = "dashed") + 
  geom_smooth(data=df2,aes(x = age, y = weight), method="lm", se=F, linetype = "dashed") + 
  geom_smooth(data=df3,aes(x = age, y = weight), method="lm", se=F,linetype = "dashed") +
  labs(x = "Age (day)", y = "Body weight (mm)") + 
  scale_colour_manual(values=c("royalblue4", "firebrick"))+
  theme(text=element_text(size=14),
        axis.title= element_text(color="black",size=15),
        axis.text= element_text(color="black",size=14), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        legend.position="right", 
        plot.title = element_text(size=16, hjust = 0.5),
        panel.border = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        legend.key = element_rect(colour = NA, fill = NA),
        strip.background = element_rect(fill="white", color = "white"))+
  scale_x_continuous(limits = c(0, 105), breaks = c(0, 14, 28, 42, 56, 70,84,98))+
  scale_y_continuous(limits = c(0, 34), breaks = c(0, 10, 20, 30))+
  facet_wrap(~ condition, scales = "free")+
  ggtitle("Male body weight")+
  guides(shape = F)

tiff('weight_plot_male.tiff', units="in", width=8, height=6, res=300, compression = 'lzw')
male_weight_plot
dev.off()

rm(df,df1,df2,df3)


####SPF length####
#data organization
df <- spf_length_df %>% rename(age=Age) %>% 
  filter(sex=="M") 

df1 <- df %>% filter(age<=28) 
df2 <- df %>% filter(age>=28 & age <=56) 
df3 <- df %>% filter (age>=56) 


#plot
male_length_plot.spf <- ggplot(df, aes(x = age, y = length, color=genotype,shape=genotype)) +
  geom_point(size=1.2,shape=1) +
  geom_smooth(data=df1,aes(x = age, y = length), method="lm", se=F, linetype = "dashed") + 
  geom_smooth(data=df2,aes(x = age, y = length), method="lm", se=F, linetype = "dashed") + 
  geom_smooth(data=df3,aes(x = age, y = length), method="lm", se=F,linetype = "dashed") +
  labs(x = "Age (day)", y = "Body length (mm)") + 
  scale_colour_manual(values=c("#4B57A2", "#A80326"))+
  theme(text=element_text(size=14),
        axis.title= element_text(color="black",size=15),
        axis.text= element_text(color="black",size=14), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        legend.position="right", 
        plot.title = element_text(size=16, hjust = 0.5),
        panel.border = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        legend.key = element_rect(colour = NA, fill = NA),
        strip.background = element_rect(fill="white", color = "white"))+
  scale_x_continuous(limits = c(0, 105), breaks = c(0, 14, 28, 42, 56, 70,84,98))+
  scale_y_continuous(limits = c(25, 105), breaks = c(25, 50, 75, 100))+
  ggtitle("SPF")+
  guides(shape = F)

tiff('length_plot_male.spf.tiff', units="in", width=4.5, height=6, res=300, compression = 'lzw')
male_length_plot.spf
dev.off()

rm(df,df1,df2,df3)

####GF length####
#data organization
df <- gf_length_df %>% rename(age=Age) %>% filter(sex=="M") 
df1 <- df %>% filter(age<=28) 
df2 <- df %>% filter(age>=28 & age <=56) 
df3 <- df %>% filter (age>=56) 

#plot
male_length_plot.gf <- ggplot(df, aes(x = age, y = length, color=genotype,shape=genotype)) +
  geom_point(size=1.2,shape=1) +
  geom_smooth(data=df1,aes(x = age, y = length), method="lm", se=F, linetype = "dashed") + 
  geom_smooth(data=df2,aes(x = age, y = length), method="lm", se=F, linetype = "dashed") + 
  geom_smooth(data=df3,aes(x = age, y = length), method="lm", se=F,linetype = "dashed") +
  labs(x = "Age (day)", y = "Body length (mm)") + 
  scale_colour_manual(values=c("#2D6A64", "#D58337"))+
  theme(text=element_text(size=14),
        axis.title= element_text(color="black",size=15),
        axis.text= element_text(color="black",size=14), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        legend.position="right", 
        plot.title = element_text(size=16, hjust = 0.5),
        panel.border = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        legend.key = element_rect(colour = NA, fill = NA),
        strip.background = element_rect(fill="white", color = "white"))+
  scale_x_continuous(limits = c(0, 105), breaks = c(0, 14, 28, 42, 56, 70,84,98))+
  scale_y_continuous(limits = c(25, 105), breaks = c(25, 50, 75, 100))+
  ggtitle("GF")+
  guides(shape = F)

tiff('length_plot_male.gf.tiff', units="in", width=4.5, height=6, res=300, compression = 'lzw')
male_length_plot.gf
dev.off()

rm(df,df1,df2,df3)



####SPF weight####
df <- spf_weight_df %>% rename(age=Age) %>% filter(sex=="M")

df1 <- df %>% filter(age<=28)
df2 <- df %>% filter(age>=28 & age <=56) 
df3 <- df %>% filter (age>=56) 


#plot
male_weight_plot.spf <- ggplot(df, aes(x = age, y = weight, color=genotype,shape=genotype)) +
  geom_point(size=1.2,shape=1) +
  geom_smooth(data=df1,aes(x = age, y = weight), method="lm", se=F, linetype = "dashed") + 
  geom_smooth(data=df2,aes(x = age, y = weight), method="lm", se=F, linetype = "dashed") + 
  geom_smooth(data=df3,aes(x = age, y = weight), method="lm", se=F,linetype = "dashed") +
  labs(x = "Age (day)", y = "Body weight (mm)") + 
  scale_colour_manual(values=c("#4B57A2", "#A80326"))+
  theme(text=element_text(size=14),
        axis.title= element_text(color="black",size=15),
        axis.text= element_text(color="black",size=14), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        legend.position="right", 
        plot.title = element_text(size=16, hjust = 0.5),
        panel.border = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        legend.key = element_rect(colour = NA, fill = NA),
        strip.background = element_rect(fill="white", color = "white"))+
  scale_x_continuous(limits = c(0, 105), breaks = c(0, 14, 28, 42, 56, 70,84,98))+
  scale_y_continuous(limits = c(0, 34), breaks = c(0, 10, 20, 30))+
  ggtitle("SPF")+
  guides(shape = F)

tiff('weight_plot_male.spf.tiff', units="in", width=4.5, height=6, res=300, compression = 'lzw')
male_weight_plot.spf
dev.off()

rm(df,df1,df2,df3)

####GF weight####
df <- gf_weight_df %>% rename(age=Age) %>% filter(sex=="M")

df1 <- df %>% filter(age<=28) 
df2 <- df %>% filter(age>=28 & age <=56) 
df3 <- df %>% filter (age>=56)

#plot
male_weight_plot.gf <- ggplot(df, aes(x = age, y = weight, color=genotype,shape=genotype)) +
  geom_point(size=1.2,shape=1) +
  geom_smooth(data=df1,aes(x = age, y = weight), method="lm", se=F, linetype = "dashed") + 
  geom_smooth(data=df2,aes(x = age, y = weight), method="lm", se=F, linetype = "dashed") + 
  geom_smooth(data=df3,aes(x = age, y = weight), method="lm", se=F,linetype = "dashed") +
  labs(x = "Age (day)", y = "Body weight (mm)") + 
  scale_colour_manual(values=c("#2D6A64", "#D58337"))+
  theme(text=element_text(size=14),
        axis.title= element_text(color="black",size=15),
        axis.text= element_text(color="black",size=14), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        legend.position="right", 
        plot.title = element_text(size=16, hjust = 0.5),
        panel.border = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        legend.key = element_rect(colour = NA, fill = NA),
        strip.background = element_rect(fill="white", color = "white"))+
  scale_x_continuous(limits = c(0, 105), breaks = c(0, 14, 28, 42, 56, 70,84,98))+
  scale_y_continuous(limits = c(0, 34), breaks = c(0, 10, 20, 30))+
  ggtitle("GF")+
  guides(shape = F)

tiff('weight_plot_male.gf.tiff', units="in", width=4.5, height=6, res=300, compression = 'lzw')
male_weight_plot.gf
dev.off()

rm(df,df1,df2,df3)


























