### This script file will read in the SAM model output across sites
### and make summary graphs

# Load packages
library(tidyverse)
library(gridExtra)
library(ggforce) # for facet grids
library(cowplot) # for ggsave
library(ggh4x) # for facet_nested
# Load self-made functions
source("./scripts/functions.R")

# Create necessary folders if they do not already exist
if(!file.exists("plots")) { dir.create("plots")}

path_out = "./plots/" # set save path

# mod1
sitename = c("lae", "crk", "geb", "lnf", "yat")
varname = c("ET", "GPP", "SIF_O2A", "WUE_GPP", "WUE_SIF", "T_TEA", "WUE_GPP_TEA", "WUE_SIF_TEA")
df_list <- list()
df_list2 <- list()
for(i in 1:length(sitename)){
  for(j in 1:length(varname)){
    
    #dffilename <- paste("./models/model1/", sitename[i],"/coda/df_mod1_", sitename[i], "_", varname[j], ".csv", sep = "")
    dffilename <- paste("./output_model1/df/df_", sitename[i], "_", varname[j],"_hour", ".csv", sep = "")
    df_temp <- tryCatch(read.csv(dffilename), error=function(err) NA)
    
    if(is.na(df_temp)){
      next
    }
    
    df_temp <- df_temp %>%
      mutate(response = varname[j], site = sitename[i])
    
    df_list[[j]] = df_temp
  }
  df_list2[[i]] = bind_rows(df_list)
}
df = bind_rows(df_list2)

# mod2
sitename = c("lae", "crk", "geb", "lnf", "yat")
varname = c("WUE_GPP", "WUE_SIF")
df_list <- list()
df_list2 <- list()
for(i in 1:length(sitename)){
  for(j in 1:length(varname)){
    
    #dffilename <- paste("./models/model1/", sitename[i],"/coda/df_mod1_", sitename[i], "_", varname[j], ".csv", sep = "")
    dffilename <- paste("./output_model2/df/df_", sitename[i], "_", varname[j],"_hour", ".csv", sep = "")
    df_temp <- tryCatch(read.csv(dffilename), error=function(err) NA)
    
    if(is.na(df_temp)){
      next
    }
    
    df_temp <- df_temp %>%
      mutate(response = varname[j], site = sitename[i])
    
    df_list[[j]] = df_temp
  }
  df_list2[[i]] = bind_rows(df_list)
}
df2 = bind_rows(df_list2)



p <- df %>%
  filter(var %in% c("wSs","wT","wV", "wPAR")) %>%
  ggplot(aes(x=ID1, y=mean)) + 
  geom_pointrange(aes(ymax = pc97.5, ymin = pc2.5), position = position_dodge(width = 1), fatten = .5) +
  #geom_point(pch=21, size = .5) +
  scale_color_brewer(palette = "Dark2") +
  #facet_grid(param~var) +
  facet_nested(site ~ response + var, scales="free_x") + 
  #geom_text(data = df_p, aes(label=w_letters), position = position_dodge(width = .7), vjust = -3, size=3) +
  labs(title = NULL, y = "weight", x = NULL)+
  theme_bw() +
  theme(legend.position = "top",
        legend.text=element_text(size=12),
        text = element_text(size=12),
        axis.text.x = element_text(size = 11, angle = 90, vjust = 0.6, hjust = 1),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.spacing=unit(.05, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        strip.background = element_rect(color = "black", fill = "white", size = 1))
p

ggsave2("p_weights.png", plot = p, path = path_out, width = 9, height = 5)

df_p <- df %>%
  filter(var == "beta1")
df_p$ID <-factor(df_p$ID1 , levels = c("1","2","3","4", "5"), labels = c("VPD", "Tair", "Sshall", "Sdeep", "PAR"))
p <- df_p %>%
  filter(ID %in% c("VPD", "Tair","PAR", "Sshall", "Sdeep")) %>%
  ggplot(aes(x=mean, y=ID)) + 
  geom_pointrange(aes(xmax = pc97.5, xmin = pc2.5, alpha = ifelse(pc2.5 <= 0 & pc97.5 >= 0, .9, 1), color = ifelse(pc2.5 < 0, "blue", "red")), position = position_dodge(width = .7), fatten = .5) +
  geom_vline(xintercept = 0, linetype = 2) + 
  #scale_color_manual(values=c("red", "tomato2", "sienna2", "orange","light blue","blue", "dark blue")) +
  #facet_grid(var ~ site, scales = "free_x") +
  #facet_row("var", scales = "free") + 
  facet_nested(site ~ response + var, scales="free_x") +
  #geom_text(data = df_p, aes(label=m_season_letters), position = position_dodge(width = .7), vjust = -.5, size=3) +
  labs(title = NULL, y = "main effects", x = NULL)+
  theme_bw() +
  theme(legend.position = "none",
        legend.text=element_text(size=10),
        text = element_text(size=10),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))
p

ggsave2("p_main_effects.png", plot = p, path = path_out, width = 8, height = 6)

df_p <- df %>%
  filter(var == "beta2")
# key: 1 VPD # 2 Tair # 3 PAR # 4 Sshall # 5 Sdeep
# X1a = cbind(X1[,1]^2, X1[,2]^2)
# X2 = X1[,1]*X1[,2], X1[,1]*X1[,4], X1[,1]*X1[,5], X1[,2]*X1[,4], X1[,2]*X1[,5], X1[,4]*X1[,5]
df_p$ID <-factor(df_p$ID1 , levels = c("1","2","3","4","5", "6"), labels = c("VPD*Tair", "Tair*Sdeep", "VPD*Sdeep", "PAR*VPD", "PAR*Tair", "PAR*Sdeep"))
p <- df_p %>%
  ggplot(aes(x=mean, y=ID)) + 
  geom_pointrange(aes(xmax = pc97.5, xmin = pc2.5, alpha = ifelse(pc2.5 <= 0 & pc97.5 >= 0, .9, 1)), position = position_dodge(width = .7), fatten = .5) +
  geom_vline(xintercept = 0, linetype = 2) + 
  #scale_color_manual(values=c("red", "tomato2", "sienna2", "orange","light blue","blue", "dark blue")) +
  #facet_grid(var ~ site, scales = "free_x") +
  #facet_row("var", scales = "free") + 
  facet_nested(site ~ response + var, scales="free_x") +
  #geom_text(data = df_p, aes(label=m_season_letters), position = position_dodge(width = .7), vjust = -.5, size=3) +
  labs(title = NULL, y = "main effects", x = NULL)+
  theme_bw() +
  theme(legend.position = "none",
        legend.text=element_text(size=10),
        text = element_text(size=10),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))
p

ggsave2("p_int_effects.png", plot = p, path = path_out, width = 8, height = 6)



df_p <- df %>%
  filter(var == "dYdX")
df_p$ID2 <-factor(df_p$ID2 , levels = c("1","2","3","4", "5"), labels = c("VPD", "Tair", "Sshall", "Sdeep", "PAR"))
df_p$TIMESTAMP <- as.POSIXct(df_p$TIMESTAMP)

for(i in 1:length(sitename)){
  p <- df_p %>%
    filter(site==sitename[i]) %>%
    ggplot(aes(x=ID1, y=mean)) +
    geom_pointrange(aes(ymax = pc97.5, ymin = pc2.5, color = ifelse(pc2.5 <= 0 & pc97.5 >= 0, "nonsignificant", "significant")), position = position_dodge(width = 1), fatten = .5) +
    scale_color_manual(values = c("gray", "black")) +
    geom_hline(yintercept = 0, linetype = 2) +
    #geom_rect(aes(xmin = date, xmax = dplyr::lead(date), ymin = -Inf, ymax = Inf, fill = Season), alpha = 0.1) +
    #scale_fill_brewer(palette = "Dark2") +
    facet_grid(response~ID2, scales = "free") +
    #facet_nested_wrap(ID2 ~ response + site, scales="free", nrow=4) +
    labs(title = sitename[i], y = "dY/dX", x = "Half-hour Timestep")+
    theme_bw() +
    theme(legend.position = "top",
          legend.text=element_text(size=10),
          axis.text.x = element_blank(),
          #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          text = element_text(size=10),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5))

  ggname <- paste("p_netsens_", i,".png", sep="")
  ggsave2(ggname, plot = p, path = path_out, width = 14, height = 8)
  
}



# Bayesian R2

p <- df %>%
  filter(var == "R2") %>%
  ggplot(aes(x=site, y=mean)) +
  geom_pointrange(aes(ymax = pc97.5, ymin = pc2.5), position = position_dodge(width = 1), fatten = .5) +
  facet_row("response") +
  #facet_nested_wrap(ID2 ~ response + site, scales="free", nrow=4) +
  labs(title = NULL, y = "R2", x = "Site")+
  theme_bw() +
  theme(legend.position = "top",
        legend.text=element_text(size=10),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size=10),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))
p

ggsave2("R2.png", plot = p, path = path_out, width = 8, height = 6)


# Decoupling timeseries graph

for(i in 1:length(sitename)){
  dfname <- paste("./data_clean/",sitename[i],"dat.RData", sep="")
  load(dfname) #dat
  
  dat$TIMESTAMP <- as.POSIXct(dat$TIMESTAMP)
  
  p1 <- dat %>%
    ggplot(aes(x=TIMESTAMP)) +
    geom_line(aes(y = T_TEA), position = position_dodge(width = 1)) +
    #geom_line(aes(y = P, color = "P"), position = position_dodge(width = 1)) +
    labs(title = NULL, y = "T", x = NULL)+
    theme_bw() +
    theme(legend.position = "top",
          legend.text=element_text(size=10),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          text = element_text(size=10),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5))
  p1.5 <- dat %>%
    ggplot(aes(x=TIMESTAMP)) +
    geom_line(aes(y = ET), position = position_dodge(width = 1)) +
    #geom_line(aes(y = P, color = "P"), position = position_dodge(width = 1)) +
    labs(title = NULL, y = "ET", x = NULL)+
    ylim(c(0,30))+
    theme_bw() +
    theme(legend.position = "top",
          legend.text=element_text(size=10),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          text = element_text(size=10),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5))
  p2 <- dat %>%
    ggplot(aes(x=TIMESTAMP)) +
    geom_line(aes(y = GPP), position = position_dodge(width = 1)) +
    labs(title = NULL, y = "GPP", x = NULL)+
    theme_bw() +
    theme(legend.position = "top",
          legend.text=element_text(size=10),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          text = element_text(size=10),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5))
  p3 <- dat %>%
    ggplot(aes(x=TIMESTAMP)) +
    geom_line(aes(y = SIF_O2A), position = position_dodge(width = 1)) +
    labs(title = NULL, y = "far-red SIF", x = NULL)+
    theme_bw() +
    theme(legend.position = "top",
          legend.text=element_text(size=10),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          text = element_text(size=10),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5))
  p4 <- dat %>%
    ggplot(aes(x=TIMESTAMP)) +
    #geom_line(aes(y = T_TEA), position = position_dodge(width = 1)) +
    geom_line(aes(y = P), position = position_dodge(width = 1)) +
    labs(title = NULL, y = "P", x = NULL)+
    theme_bw() +
    theme(legend.position = "top",
          legend.text=element_text(size=10),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          text = element_text(size=10),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5))
  p <- grid.arrange(p1,p1.5,p2,p3,p4,nrow=5)

  savename <- paste(sitename[i],"view.png",sep="")
  ggsave2(savename, plot = p, path = path_out, width = 9, height = 5)
  
}



################## model 2



p <- df2 %>%
  filter(var %in% c("wSs","wT","wV", "wPAR")) %>%
  ggplot(aes(x=ID1, y=mean)) + 
  geom_pointrange(aes(ymax = pc97.5, ymin = pc2.5), position = position_dodge(width = 1), fatten = .5) +
  #geom_point(pch=21, size = .5) +
  scale_color_brewer(palette = "Dark2") +
  #facet_grid(param~var) +
  facet_nested(site ~ response + var, scales="free_x") + 
  #geom_text(data = df_p, aes(label=w_letters), position = position_dodge(width = .7), vjust = -3, size=3) +
  labs(title = NULL, y = "weight", x = NULL)+
  theme_bw() +
  theme(legend.position = "top",
        legend.text=element_text(size=12),
        text = element_text(size=12),
        axis.text.x = element_text(size = 11, angle = 90, vjust = 0.6, hjust = 1),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.spacing=unit(.05, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        strip.background = element_rect(color = "black", fill = "white", size = 1))
p

ggsave2("p_weights_mod2.png", plot = p, path = path_out, width = 8, height = 6)

df_p <- df2 %>%
  filter(var == "beta1")
df_p$ID <-factor(df_p$ID1 , levels = c("1","2","3","4"), labels = c("VPD", "Tair", "Sshall", "Sdeep"))
p <- df_p %>%
  filter(ID %in% c("VPD", "Tair","PAR", "Sshall", "Sdeep")) %>%
  ggplot(aes(x=mean, y=ID)) + 
  geom_pointrange(aes(xmax = pc97.5, xmin = pc2.5, alpha = ifelse(pc2.5 <= 0 & pc97.5 >= 0, .9, 1), color = ifelse(pc2.5 < 0, "blue", "red")), position = position_dodge(width = .7), fatten = .5) +
  geom_vline(xintercept = 0, linetype = 2) + 
  #scale_color_manual(values=c("red", "tomato2", "sienna2", "orange","light blue","blue", "dark blue")) +
  #facet_grid(var ~ site, scales = "free_x") +
  #facet_row("var", scales = "free") + 
  facet_nested(site ~ response + var, scales="free_x") +
  #geom_text(data = df_p, aes(label=m_season_letters), position = position_dodge(width = .7), vjust = -.5, size=3) +
  labs(title = NULL, y = "main effects", x = NULL)+
  theme_bw() +
  theme(legend.position = "none",
        legend.text=element_text(size=10),
        text = element_text(size=10),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))
p

ggsave2("p_main_effects_mod2.png", plot = p, path = path_out, width = 8, height = 6)

df_p <- df2 %>%
  filter(var == "beta2")
# key: 1 VPD # 2 Tair # 3 PAR # 4 Sshall # 5 Sdeep
# X1a = cbind(X1[,1]^2, X1[,2]^2)
# X2 = X1[,1]*X1[,2], X1[,1]*X1[,4], X1[,1]*X1[,5], X1[,2]*X1[,4], X1[,2]*X1[,5], X1[,4]*X1[,5]
df_p$ID <-factor(df_p$ID1 , levels = c("1","2","3","4","5", "6"), labels = c("VPD*Tair", "Tair*Sdeep", "VPD*Sdeep", "PAR*VPD", "PAR*Tair", "PAR*Sdeep"))
p <- df_p %>%
  ggplot(aes(x=mean, y=ID)) + 
  geom_pointrange(aes(xmax = pc97.5, xmin = pc2.5, alpha = ifelse(pc2.5 <= 0 & pc97.5 >= 0, .9, 1)), position = position_dodge(width = .7), fatten = .5) +
  geom_vline(xintercept = 0, linetype = 2) + 
  #scale_color_manual(values=c("red", "tomato2", "sienna2", "orange","light blue","blue", "dark blue")) +
  #facet_grid(var ~ site, scales = "free_x") +
  #facet_row("var", scales = "free") + 
  facet_nested(site ~ response + var, scales="free_x") +
  #geom_text(data = df_p, aes(label=m_season_letters), position = position_dodge(width = .7), vjust = -.5, size=3) +
  labs(title = NULL, y = "main effects", x = NULL)+
  theme_bw() +
  theme(legend.position = "none",
        legend.text=element_text(size=10),
        text = element_text(size=10),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))
p

ggsave2("p_int_effects_mod2.png", plot = p, path = path_out, width = 8, height = 6)


df_p <- df2 %>%
  filter(var == "dYdX")
df_p$ID2 <-factor(df_p$ID2 , levels = c("1","2","3","4", "5"), labels = c("VPD", "Tair", "Sshall", "Sdeep", "PAR"))

sitename <- c("yat", "crk", "lae", "lnf")
for(i in 2:length(sitename)){
  p <- df_p %>%
    filter(site==sitename[i]) %>%
    ggplot(aes(x=ID1, y=mean)) +
    geom_pointrange(aes(ymax = pc97.5, ymin = pc2.5, color = ifelse(pc2.5 <= 0 & pc97.5 >= 0, "nonsignificant", "significant")), position = position_dodge(width = 1), fatten = .5) +
    scale_color_manual(values = c("gray", "black")) +
    geom_hline(yintercept = 0, linetype = 2) +
    #geom_rect(aes(xmin = date, xmax = dplyr::lead(date), ymin = -Inf, ymax = Inf, fill = Season), alpha = 0.1) +
    #scale_fill_brewer(palette = "Dark2") +
    facet_grid(response~ID2, scales = "free") +
    #facet_nested_wrap(ID2 ~ response + site, scales="free", nrow=4) +
    labs(title = sitename[i], y = "dY/dX", x = "Hour Timestep")+
    theme_bw() +
    theme(legend.position = "top",
          legend.text=element_text(size=10),
          axis.text.x = element_blank(),
          #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          text = element_text(size=10),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5))
  #p
  
  ggname <- paste("p_netsens_mod2_", i,".png", sep="")
  ggsave2(ggname, plot = p, path = path_out, width = 14, height = 8)
  
}

#################################################### Compare Ts

tempsite <- "geb"

df_p <- df2 %>%
  filter(site==tempsite) %>%
  filter(var=="T.pred") %>%
  filter(response=="WUE_GPP")
df_p$TIMESTAMP <-as.POSIXct(df_p$TIMESTAMP, format = "%Y-%m-%d %H:%M:%S")

df_p$year <- format(df_p$TIMESTAMP, format = "%Y")
df_p$month <- format(df_p$TIMESTAMP, format = "%m")
df_p$day <- format(df_p$TIMESTAMP, format = "%d")
df_p$hour <- format(df_p$TIMESTAMP, format = "%H")

load(paste("./data_clean/", tempsite,"dat.RData", sep="")) #dat
dat2 <- dat

df_p2 <- dat2 %>%
  select(T_TEA, TIMESTAMP,ET)
df_p2$TIMESTAMP <- as.POSIXct(df_p2$TIMESTAMP, format = "%Y-%m-%d %H:%M:%S")
df_p2$year <- format(df_p2$TIMESTAMP, format = "%Y")
df_p2$month <- format(df_p2$TIMESTAMP, format = "%m")
df_p2$day <- format(df_p2$TIMESTAMP, format = "%d")
df_p2$hour <- format(df_p2$TIMESTAMP, format = "%H")

df_p3 <- left_join(df_p, df_p2, by = c("year","month","day","hour"))
df_p3$T_TEA_ratio <- df_p3$T_TEA/df_p3$ET
df_p3$T_ratio <- df_p3$mean/df_p3$ET

p <- df_p3 %>%
  filter(DOY > 160) %>%
  ggplot(aes(x=TIMESTAMP.x)) +
  geom_line(aes(y=mean, color = "DEPART T")) +
  geom_line(aes(y=T_TEA, color = "TEA T")) +
  scale_color_manual(values = c("red", "blue")) +
  labs(title = NULL, y = "T", x = "Hour Timestep")+
  theme_bw() +
  theme(legend.position = "top",
        legend.text=element_text(size=14),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size=14),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))
p

ggsave2(paste("p_T_comp_",tempsite,".png", sep=""), plot = p, path = path_out, width = 8, height = 5)

summary(lm(T_TEA ~ mean, data = df_p3))

library(ggpubr)
p <- df_p3 %>%
  ggplot(aes(x = T_TEA, y= mean)) +
  #geom_point() +
  geom_pointrange(aes(ymin=pc2.5, ymax=pc97.5), alpha=0.5)+
  geom_smooth(method="lm", se = F, color = "red") +
  geom_abline(slope=1, intercept=0, lty=2, col="blue", size=1.25)+
  stat_cor(aes(label = ..rr.label..), color = "red", geom = "label") +
  #stat_cor(aes(label = ..rr.label..), color = "red", label.x = 0.5, size = 3) +
  ylim(0,.6) + xlim(0,.6) +
  facet_row("site", strip.position = "top") +
  labs(title = NULL, x="TEA T", y="DEPART T") +
  #theme_classic(base_size = 12)+
  theme(legend.position = "right",
        legend.text=element_text(size=14),
        text = element_text(size=14),
        legend.title = element_blank(),
        panel.background = element_rect(fill="white"),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(colour="black"),
        aspect.ratio=1,
        plot.title = element_text(hjust = 0.5))
p
ggsave2(paste("p_T_R2_",tempsite,".png", sep=""), plot = p, path = path_out, width = 8, height = 5)


p <- df_p3 %>%
  filter(DOY > 160) %>%
  ggplot(aes(x=TIMESTAMP.x)) +
  geom_point(aes(y=T_ratio), color = "red") +
  geom_point(aes(y=T_TEA_ratio), color = "blue") +
  labs(title = NULL, y = "T/ET", x = "Hour Timestep")+
  theme_bw() +
  theme(legend.position = "top",
        legend.text=element_text(size=14),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size=14),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))
p

########################## Selected plots for presentations/posters


df_plot <- df %>%
  filter(var == "dYdX") %>%
  filter(response %nin% c("WUE_GPP", "WUE_SIF","WUE_GPP_TEA", "WUE_SIF_TEA", "ET"))%>%
  mutate(mod = 1)
df_plot2 <- df2 %>%
  filter(var == "dYdX") %>%
  mutate(mod = 2)

# df_p$TIMESTAMP <- as.character(df_p$TIMESTAMP) # timestamps should be in chatacter format for joining
# df_p <- left_join(df_p, df_wue, by = c("DOY", "site", "TIMESTAMP", "response"))

df_plot <- rbind(df_plot, df_plot2)
df_plot$ID2 <-factor(df_plot$ID2 , levels = c("1","2","3","4", "5"), labels = c("VPD", "Tair", "Sshall", "Sdeep", "PAR"))
df_plot$response <-factor(df_plot$response , levels = c("WUE_GPP", "WUE_SIF","GPP","SIF_O2A", "T_TEA"), labels = c("GPP/T (WUE)", "SIF/T (WUE)", "GPP", "SIF", "T_TEA"))
df_plot$cat <- ifelse(df_plot$ID2 %in% c("VPD", "Tair", "PAR"), "atmospheric driver", "soil moisture driver")
df_p <- df_plot
df_p$TIMESTAMP <- as.POSIXct(df_p$TIMESTAMP)


sitename = c("lae", "crk", "geb", "lnf", "yat")
df_list <- list()
for(i in c(1:length(sitename))){
  print(paste("i: ", i, sep=""))
  # load org data for each site
  load(paste("./data_clean/",sitename[i],"dat.RData", sep="")) # dat
  
  # aggregate by hour
  dat <- dat %>%
    mutate(hour = format(dat$TIMESTAMP, format ="%H"), DOY = floor(DOY)) %>%
    group_by(DOY,hour) %>%
    summarise(across(everything(), mean, na.rm = T)) %>%
    ungroup()
  
  dat$site <- sitename[i]
  
  df_list[[i]] <- dat
}
df_clean <- bind_rows(df_list)

# combine clean and output, but make sure TIMESTAMPS are in correct format
#df_clean$TIMESTAMP <- as.character(df_clean$TIMESTAMP)
df_clean$year <- format(df_clean$TIMESTAMP, format = "%Y")
df_clean$month <- format(df_clean$TIMESTAMP, format = "%m")
df_clean$day <- format(df_clean$TIMESTAMP, format = "%d")
df_clean$hour <- format(df_clean$TIMESTAMP, format = "%H")

#df_p$TIMESTAMP <- as.character(df_p$TIMESTAMP) # timestamps should be in chatacter format for joining
df_p$year <- format(df_p$TIMESTAMP, format = "%Y")
df_p$month <- format(df_p$TIMESTAMP, format = "%m")
df_p$day <- format(df_p$TIMESTAMP, format = "%d")
df_p$hour <- format(df_p$TIMESTAMP, format = "%H")
df_p <- left_join(df_p, df_clean, by = c("site", "year", "month", "day", "hour"))
df_p$time <- df_p$hour
df_p$TIMESTAMP <- as.POSIXct(df_p$TIMESTAMP.x)


p <- df_p %>%
  filter(site=="crk") %>%
  ggplot(aes(x=TIMESTAMP, y=mean)) +
  geom_pointrange(aes(ymax = pc97.5, ymin = pc2.5, color = ifelse(pc2.5 <= 0 & pc97.5 >= 0, "nonsignificant", "significant")), position = position_dodge(width = 1), fatten = .5) +
  scale_color_manual(values = c("gray", "black")) +
  geom_hline(yintercept = 0, linetype = 2) +
  #geom_rect(aes(xmin = date, xmax = dplyr::lead(date), ymin = -Inf, ymax = Inf, fill = Season), alpha = 0.1) +
  #scale_fill_brewer(palette = "Dark2") +
  #facet_grid(response~ID2, scales = "free") +
  facet_nested(response ~ cat + ID2, scales="free") +
  xlim(c(as.POSIXct("2016-05-28 07:30:00", format="%Y-%m-%d %H:%M:%S"),as.POSIXct("2016-06-06 16:00:00", format="%Y-%m-%d %H:%M:%S"))) +
  labs(title = "KR-Crk Rice Paddy Cropland", y = "dY/dX", x = "Half-hour Timestep in 2016")+
  theme_bw() +
  theme(legend.position = "top",
        legend.text=element_text(size=14),
        #axis.text.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        text = element_text(size=14),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.spacing=unit(.05, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        strip.background = element_rect(color = "black", fill = "white", size = 1))
p

ggsave2("p_poster_net_crk.png", plot = p, path = path_out, width = 8, height = 6)

# Take average of net sensitivities by time
df_p_grouped <- df_p %>%
  group_by(site,cat,response,time,ID2) %>%
  summarise_at(vars(mean, median, pc2.5, pc97.5), list(avg = mean))
df_p_grouped$time <- as.numeric(df_p_grouped$time)

p <- df_p_grouped %>% #CRK
  filter(site=="crk")%>%
  filter(response %in% c("GPP/T (WUE)", "SIF/T (WUE)")) %>%
  ggplot(aes(x=time, y=mean_avg)) +
  geom_pointrange(aes(ymax = pc97.5_avg, ymin = pc2.5_avg), color= "black", fatten = .7) +
  geom_pointrange(aes(ymax = pc97.5_avg, ymin = pc2.5_avg, color= time), fatten = .2) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_color_viridis_c(option = "D") +
  #geom_rect(aes(xmin = date, xmax = dplyr::lead(date), ymin = -Inf, ymax = Inf, fill = Season), alpha = 0.1) +
  #scale_fill_brewer(palette = "Dark2") +
  #facet_grid(response~ID2, scales = "free") +
  facet_nested(response ~ cat + ID2, scales = "free") +
  #xlim(c(as.POSIXct("2016-05-28 07:30:00", format="%Y-%m-%d %H:%M:%S"),as.POSIXct("2016-06-06 16:00:00", format="%Y-%m-%d %H:%M:%S"))) +
  labs(title = NULL, y = "Average Diurnal dY/dX", x = "Hour of Day")+
  theme_bw() +
  theme(legend.position = "none",
        legend.text=element_text(size=24),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        text = element_text(size=24),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.spacing=unit(.05, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        strip.background = element_rect(color = "black", fill = "white", size = 1))
p
ggsave2("p_poster_net_crk3.png", plot = p, path = path_out, width = 9, height = 5.5)


p <- df_p_grouped %>% #YAT
  filter(site=="geb")%>%
  filter(response %in% c("GPP/T (WUE)", "SIF/T (WUE)")) %>%
  ggplot(aes(x=time, y=median_avg, color= time)) +
  geom_point() +
  geom_pointrange(aes(ymax = pc97.5_avg, ymin = pc2.5_avg), color= "black", fatten = .7) +
  geom_pointrange(aes(ymax = pc97.5_avg, ymin = pc2.5_avg, color= time), fatten = .2) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_color_viridis_c(option = "D") +
  #geom_rect(aes(xmin = date, xmax = dplyr::lead(date), ymin = -Inf, ymax = Inf, fill = Season), alpha = 0.1) +
  #scale_fill_brewer(palette = "Dark2") +
  #facet_grid(response~ID2, scales = "free") +
  facet_nested(response ~ cat + ID2, scales = "free") +
  labs(title = NULL, y = "Average Diurnal dY/dX", x = "Hour of Day")+
  theme_bw() +
  theme(legend.position = "none",
        legend.text=element_text(size=24),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        text = element_text(size=24),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.spacing=unit(.05, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        strip.background = element_rect(color = "black", fill = "white", size = 1))
p
ggsave2("p_poster_net_geb.png", plot = p, path = path_out, width = 9, height = 5.5)

p <- df_p_grouped %>% #YAT
  filter(site=="lae")%>%
  filter(response %in% c("GPP/T (WUE)", "SIF/T (WUE)")) %>%
  ggplot(aes(x=time, y=median_avg, color= time)) +
  geom_point() +
  geom_pointrange(aes(ymax = pc97.5_avg, ymin = pc2.5_avg), color= "black", fatten = .7) +
  geom_pointrange(aes(ymax = pc97.5_avg, ymin = pc2.5_avg, color= time), fatten = .2) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_color_viridis_c(option = "D") +
  #geom_rect(aes(xmin = date, xmax = dplyr::lead(date), ymin = -Inf, ymax = Inf, fill = Season), alpha = 0.1) +
  #scale_fill_brewer(palette = "Dark2") +
  #facet_grid(response~ID2, scales = "free") +
  facet_nested(response ~ cat + ID2, scales = "free") +
  labs(title = NULL, y = "Average Diurnal dY/dX", x = "Hour of Day")+
  theme_bw() +
  theme(legend.position = "none",
        legend.text=element_text(size=24),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        text = element_text(size=24),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.spacing=unit(.05, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        strip.background = element_rect(color = "black", fill = "white", size = 1))
p
ggsave2("p_poster_net_lae.png", plot = p, path = path_out, width = 9, height = 5.5)


p <- df_p_grouped %>% # temp
  filter(site=="crk")%>%
  filter(response %in% c("GPP", "SIF")) %>%
  ggplot(aes(x=time, y=mean_avg)) +
  geom_pointrange(aes(ymax = pc97.5_avg, ymin = pc2.5_avg), color= "black", fatten = .7) +
  geom_pointrange(aes(ymax = pc97.5_avg, ymin = pc2.5_avg, color= time), fatten = .2) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_color_viridis_c(option = "D") +
  #geom_rect(aes(xmin = date, xmax = dplyr::lead(date), ymin = -Inf, ymax = Inf, fill = Season), alpha = 0.1) +
  #scale_fill_brewer(palette = "Dark2") +
  #facet_grid(response~ID2, scales = "free") +
  facet_nested(response ~ cat + ID2, scales = "free") +
  #xlim(c(as.POSIXct("2016-05-28 07:30:00", format="%Y-%m-%d %H:%M:%S"),as.POSIXct("2016-06-06 16:00:00", format="%Y-%m-%d %H:%M:%S"))) +
  labs(title = NULL, y = "Average Diurnal dY/dX", x = "Hour of Day")+
  theme_bw() +
  theme(legend.position = "none",
        legend.text=element_text(size=24),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        text = element_text(size=24),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.spacing=unit(.05, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        strip.background = element_rect(color = "black", fill = "white", size = 1))
p

ggsave2("p_poster_net_crk_phot.png", plot = p, path = path_out, width = 9, height = 5.5)

p <- df_p_grouped %>% # temp
  filter(site=="geb")%>%
  filter(response %in% c("GPP", "SIF")) %>%
  ggplot(aes(x=time, y=mean_avg, color = time)) +
  #geom_point() +
  geom_pointrange(aes(ymax = pc97.5_avg, ymin = pc2.5_avg), color= "black", fatten = .7) +
  geom_pointrange(aes(ymax = pc97.5_avg, ymin = pc2.5_avg, color= time), fatten = .2) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_color_viridis_c(option = "D") +
  #scale_color_distiller(palette = "YlOrBr") +
  #geom_rect(aes(xmin = date, xmax = dplyr::lead(date), ymin = -Inf, ymax = Inf, fill = Season), alpha = 0.1) +
  #scale_fill_brewer(palette = "Dark2") +
  #facet_grid(response~ID2, scales = "free") +
  facet_nested(response ~ cat + ID2, scales = "free") +
  #xlim(c(as.POSIXct("2016-05-28 07:30:00", format="%Y-%m-%d %H:%M:%S"),as.POSIXct("2016-06-06 16:00:00", format="%Y-%m-%d %H:%M:%S"))) +
  labs(title = NULL, y = "Average Diurnal dY/dX", x = "Hour of Day")+
  theme_bw() +
  theme(legend.position = "none",
        legend.text=element_text(size=24),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        text = element_text(size=24),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.spacing=unit(.05, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        strip.background = element_rect(color = "black", fill = "white", size = 1))
p

ggsave2("p_poster_net_geb_phot.png", plot = p, path = path_out, width = 9, height = 5.5)

p <- df_p %>%
  filter(site=="geb") %>%
  filter(response %in% c("GPP", "SIF")) %>%
  #filter(ID2 %in% c("Tair")) %>%
  ggplot(aes(x=TIMESTAMP, y=mean)) +
  geom_pointrange(aes(ymax = pc97.5, ymin = pc2.5, color = ifelse(pc2.5 <= 0 & pc97.5 >= 0, "nonsignificant", "significant")), position = position_dodge(width = 1), fatten = .5) +
  scale_color_manual(values = c("gray", "black")) +
  geom_hline(yintercept = 0, linetype = 2) +
  #geom_rect(aes(xmin = date, xmax = dplyr::lead(date), ymin = -Inf, ymax = Inf, fill = Season), alpha = 0.1) +
  facet_grid(response ~ ID2, scales="free") +
  #xlim(c(as.POSIXct("2016-05-28 07:30:00", format="%Y-%m-%d %H:%M:%S"),as.POSIXct("2016-06-06 16:00:00", format="%Y-%m-%d %H:%M:%S"))) +
  labs(title = "IL-Yat Mediterranean Evergreen Coniferous Forest", y = "dY/dX", x = "Half-hour Timestep in 2017")+
  theme_bw() +
  theme(legend.position = "top",
        legend.text=element_text(size=12),
        #axis.text.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        text = element_text(size=14),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.spacing=unit(.06, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        strip.background = element_rect(color = "black", fill = "white", size = 1))
p

ggsave2("p_poster_net_geb.png", plot = p, path = path_out, width = 7, height = 4)


p <- df_p %>%
  filter(site=="lae") %>%
  filter(response %in% c("GPP", "SIF")) %>%
  #filter(ID2 %in% c("Tair")) %>%
  ggplot(aes(x=TIMESTAMP, y=mean)) +
  geom_pointrange(aes(ymax = pc97.5, ymin = pc2.5, color = ifelse(pc2.5 <= 0 & pc97.5 >= 0, "nonsignificant", "significant")), position = position_dodge(width = 1), fatten = .5) +
  scale_color_manual(values = c("gray", "black")) +
  geom_hline(yintercept = 0, linetype = 2) +
  #geom_rect(aes(xmin = date, xmax = dplyr::lead(date), ymin = -Inf, ymax = Inf, fill = Season), alpha = 0.1) +
  facet_grid(response ~ ID2, scales="free") +
  #xlim(c(as.POSIXct("2018-07-01 07:15:00", format="%Y-%m-%d %H:%M:%S"),as.POSIXct("2018-08-30 16:00:00", format="%Y-%m-%d %H:%M:%S"))) +
  labs(title = "CH-Lae Mixed Deciduous Mountain Forest", y = "dY/dX", x = "Half-hour Timestep in 2018")+
  theme_bw() +
  theme(legend.position = "top",
        legend.text=element_text(size=12),
        #axis.text.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        text = element_text(size=14),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.spacing=unit(.06, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        strip.background = element_rect(color = "black", fill = "white", size = 1))
p

ggsave2("p_poster_net_lae.png", plot = p, path = path_out, width = 7, height = 4)


p <- df_p %>%
  filter(site=="geb") %>%
  filter(response %in% c("GPP", "SIF")) %>%
  #filter(ID2 %in% c("Tair")) %>%
  ggplot(aes(x=TIMESTAMP, y=mean)) +
  geom_pointrange(aes(ymax = pc97.5, ymin = pc2.5, color = ifelse(pc2.5 <= 0 & pc97.5 >= 0, "nonsignificant", "significant")), position = position_dodge(width = 1), fatten = .5) +
  scale_color_manual(values = c("gray", "black")) +
  geom_hline(yintercept = 0, linetype = 2) +
  #geom_rect(aes(xmin = date, xmax = dplyr::lead(date), ymin = -Inf, ymax = Inf, fill = Season), alpha = 0.1) +
  facet_grid(response ~ ID2, scales="free") +
  #xlim(c(as.POSIXct("2016-05-28 07:30:00", format="%Y-%m-%d %H:%M:%S"),as.POSIXct("2016-06-06 16:00:00", format="%Y-%m-%d %H:%M:%S"))) +
  labs(title = "DE-Geb", y = "dY/dX", x = "Half-hour Timestep in 2017")+
  theme_bw() +
  theme(legend.position = "top",
        legend.text=element_text(size=12),
        #axis.text.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        text = element_text(size=14),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.spacing=unit(.06, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        strip.background = element_rect(color = "black", fill = "white", size = 1))
p

ggsave2("p_poster_net_geb.png", plot = p, path = path_out, width = 7, height = 4)


p <- df %>%
  filter(site %in% c("lnf")) %>%
  filter(response %in% c("ET","GPP", "SIF_O2A")) %>%
  filter(var %in% c("wT","wV")) %>%
  ggplot(aes(x=ID1, y=mean)) + 
  geom_pointrange(aes(ymax = pc97.5, ymin = pc2.5), position = position_dodge(width = 1), fatten = .5) +
  #geom_point(pch=21, size = .5) +
  scale_color_brewer(palette = "Dark2") +
  facet_grid(var~response) +
  #facet_nested(site ~ response + var, scales="free_x") + 
  #geom_text(data = df_p, aes(label=w_letters), position = position_dodge(width = .7), vjust = -3, size=3) +
  labs(title = NULL, y = "weight", x = NULL)+
  theme_bw() +
  theme(legend.position = "top",
        legend.text=element_text(size=10),
        text = element_text(size=10),
        axis.text.x = element_text(size = 11, angle = 90, vjust = 0.6, hjust = 1),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))
p

ggsave2("p_weights_poster_lnf.png", plot = p, path = path_out, width = 6, height = 3.5)


p <- df %>%
  #filter(site %in% c("yat")) %>%
  filter(response %in% c("ET","GPP", "SIF_O2A")) %>%
  filter(var %in% c("wT")) %>%
  ggplot(aes(x=ID1, y=mean)) + 
  geom_pointrange(aes(ymax = pc97.5, ymin = pc2.5), position = position_dodge(width = 1), fatten = .5) +
  #geom_point(pch=21, size = .5) +
  scale_color_brewer(palette = "Dark2") +
  #facet_nested(var~response+site) +
  facet_nested(site~response+var) +
  #facet_nested(site ~ response + var, scales="free_x") + 
  #geom_text(data = df_p, aes(label=w_letters), position = position_dodge(width = .7), vjust = -3, size=3) +
  labs(title = NULL, y = "weight", x = NULL)+
  theme_bw() +
  theme(legend.position = "top",
        legend.text=element_text(size=10),
        text = element_text(size=10),
        axis.text.x = element_text(size = 11, angle = 90, vjust = 0.6, hjust = 1),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))
p

ggsave2("p_weights_poster_temp.png", plot = p, path = path_out, width = 6, height = 3.5)


# net sens tair

df_p1 <- df_p %>%
  filter(site=="crk") %>%
  filter(response %in% c("GPP", "SIF"))%>%
  filter(ID2=="Tair")
df_p1$time <- as.numeric(df_p1$time)
df_p2 <- df_p %>%
  filter(site=="geb") %>%
  filter(response %in% c("GPP", "SIF"))%>%
  filter(ID2=="Tair")
df_p2$time <- as.numeric(df_p2$time)

library(ggnewscale)
p <- ggplot() +
  #geom_pointrange(aes(ymax = pc97.5, ymin = pc2.5, color = site), fatten=.2,alpha=.3) +
  geom_pointrange(data = df_p1, aes(x=TA, y=mean, ymax = pc97.5, ymin = pc2.5, color = time), fatten=.2,alpha=.3) +
  scale_color_viridis_c(option = "D") +
  #scale_color_distiller(palette = "RdPu") +
  new_scale_color() +
  geom_pointrange(data = df_p2, aes(x=TA, y=mean, ymax = pc97.5, ymin = pc2.5, color = time), shape = 2, fatten=.2,alpha=.3) +
  scale_color_distiller(palette = "YlOrBr") +
  geom_hline(yintercept = 0, linetype = 2) +
  #facet_grid(response~ID2, scales = "free") +
  facet_wrap("response", scales = "free_y") +
  #facet_wrap(site~response, scales = "free_y") +
  labs(title = NULL, y = "dY/dX", x = "Air Temperature (°C)")+
  theme_bw() +
  theme(legend.position = "top",
        legend.text=element_text(size=14),
        #legend.title = "time",
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle , vjust = 1, hjust=1),
        text = element_text(size=14),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(color = "black", fill = "white", size = 1))
p

ggsave2("p_poster_net_tair.png", plot = p, path = path_out, width = 8, height = 6)


df_p3 <- df_p %>%
  filter(response %in% c("GPP", "SIF"))%>%
  filter(ID2=="Tair")
df_p3$time <- as.numeric(df_p3$time)
p <- ggplot(data = df_p3, aes(x=TA, y=mean)) +
  #geom_pointrange(aes(ymax = pc97.5, ymin = pc2.5, color = site), fatten=.2,alpha=.3) +
  geom_pointrange(data = df_p3 %>% filter(site=="crk"), aes(ymax = pc97.5, ymin = pc2.5, color = time, shape = site), fatten=.2,alpha=.3) +
  scale_color_viridis_c(option = "D") +
  #scale_color_distiller(palette = "RdPu") +
  new_scale_color() +
  geom_pointrange(data = df_p3 %>% filter(site=="geb"), aes(ymax = pc97.5, ymin = pc2.5, color = time, shape = site), fatten=.2,alpha=.3) +
  scale_color_distiller(palette = "YlOrBr") +
  scale_shape_manual(values = c(1,2)) +
  geom_hline(yintercept = 0, linetype = 2) +
  #facet_grid(response~ID2, scales = "free") +
  facet_wrap("response", scales = "free_y") +
  #facet_wrap(site~response, scales = "free_y") +
  labs(title = NULL, y = "dY/dX", x = "Air Temperature (°C)")+
  theme_bw() +
  theme(legend.position = "top",
        legend.text=element_text(size=18),
        #legend.title = element_blank(),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle , vjust = 1, hjust=1),
        text = element_text(size=24),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(color = "black", fill = "white", size = 1))
p

ggsave2("p_poster_net_tair2.png", plot = p, path = path_out, width = 8, height = 5)


df_p3 <- df_p %>%
  filter(response %in% c("GPP", "SIF", "T_TEA"))%>%
  filter(ID2=="Tair")
df_p3$time <- as.numeric(df_p3$time)
p <- ggplot(data = df_p3, aes(x=TA, y=mean)) +
  #geom_pointrange(aes(ymax = pc97.5, ymin = pc2.5, color = site), fatten=.2,alpha=.3) +
  geom_pointrange(data = df_p3 %>% filter(site=="lae"), aes(ymax = pc97.5, ymin = pc2.5, color = time, shape = site), fatten=.2,alpha=.3) +
  scale_color_viridis_c(option = "D") +
  #scale_color_distiller(palette = "RdPu") +
  new_scale_color() +
  geom_pointrange(data = df_p3 %>% filter(site=="geb"), aes(ymax = pc97.5, ymin = pc2.5, color = time, shape = site), fatten=.2,alpha=.3) +
  scale_color_distiller(palette = "YlOrBr") +
  scale_shape_manual(values = c(1,2)) +
  geom_hline(yintercept = 0, linetype = 2) +
  #facet_grid(response~ID2, scales = "free") +
  facet_wrap("response", scales = "free_y") +
  #facet_wrap(site~response, scales = "free_y") +
  labs(title = NULL, y = "dY/dTA", x = "Air Temperature (°C)")+
  theme_bw() +
  theme(legend.position = "top",
        legend.text=element_text(size=18),
        #legend.title = element_blank(),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle , vjust = 1, hjust=1),
        text = element_text(size=24),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(color = "black", fill = "white", size = 1))
p

ggsave2("p_poster_net_tair3.png", plot = p, path = path_out, width = 8, height = 5)


df_p3 <- df_p %>%
  filter(site %in% c("lae", "geb", "yat")) %>%
  filter(response %in% c("GPP", "SIF", "T_TEA"))%>%
  filter(ID2=="Tair")
df_p3$time <- as.numeric(df_p3$time)
p <- ggplot(data = df_p3, aes(x=TA, y=mean)) +
  #geom_pointrange(aes(ymax = pc97.5, ymin = pc2.5, color = site), fatten=.2,alpha=.3) +
  geom_pointrange(aes(ymax = pc97.5, ymin = pc2.5, color = time), fatten=.2,alpha=.3) +
  scale_color_viridis_c(option = "D") +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_grid2(site~response, scales = "free", independent = "y") +
  #facet_wrap("response", scales = "free_y") +
  labs(title = NULL, y = "dY/dTA", x = "Air Temperature (°C)", color = "hour of day")+
  theme_bw() +
  theme(legend.position = "top",
        legend.text=element_text(size=14),
        #legend.title = element_blank(),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle , vjust = 1, hjust=1),
        text = element_text(size=14),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(color = "black", fill = "white", size = 1))
p

ggsave2("p_poster_net_tair3_2.png", plot = p, path = path_out, width = 8, height = 5)

df_p3 <- df_p %>%
  filter(site %in% c("lae", "geb", "yat")) %>%
  filter(response %in% c("GPP", "SIF", "T_TEA"))%>%
  filter(ID2=="Tair")
df_p3$time <- as.numeric(df_p3$time)
p <- ggplot(data = df_p3, aes(x=TA, y=mean)) +
  #geom_pointrange(aes(ymax = pc97.5, ymin = pc2.5, color = site), fatten=.2,alpha=.3) +
  geom_pointrange(aes(ymax = pc97.5, ymin = pc2.5, color = DOY.y), fatten=.2,alpha=.3) +
  scale_color_viridis_c(option = "D") +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_grid2(site~response, scales = "free", independent = "y") +
  #facet_wrap("response", scales = "free_y") +
  labs(title = NULL, y = "dY/dTA", x = "Air Temperature (°C)", color = "DOY")+
  theme_bw() +
  theme(legend.position = "top",
        legend.text=element_text(size=14),
        #legend.title = element_blank(),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle , vjust = 1, hjust=1),
        text = element_text(size=14),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(color = "black", fill = "white", size = 1))
p
ggsave2("p_poster_net_tair_doy.png", plot = p, path = path_out, width = 8, height = 5)

df_p3 <- df_p %>%
  filter(site %in% c("lae", "geb", "crk")) %>%
  filter(response %in% c("GPP/T (WUE)", "SIF/T (WUE)"))%>%
  filter(ID2=="Tair")
df_p3$time <- as.numeric(df_p3$time)
p <- ggplot(data = df_p3, aes(x=TA, y=mean)) +
  #geom_pointrange(aes(ymax = pc97.5, ymin = pc2.5, color = site), fatten=.2,alpha=.3) +
  geom_pointrange(aes(ymax = pc97.5, ymin = pc2.5, color = DOY.y), fatten=.2,alpha=.3) +
  scale_color_viridis_c(option = "D") +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_grid2(site~response, scales = "free", independent = "y") +
  #facet_wrap("response", scales = "free_y") +
  labs(title = NULL, y = "dY/dTA", x = "Air Temperature (°C)", color = "DOY")+
  theme_bw() +
  theme(legend.position = "top",
        legend.text=element_text(size=14),
        #legend.title = element_blank(),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle , vjust = 1, hjust=1),
        text = element_text(size=14),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(color = "black", fill = "white", size = 1))
p
ggsave2("p_poster_net_tair_doy_wue.png", plot = p, path = path_out, width = 8, height = 5)


df_p3 <- df_p %>%
  filter(site %in% c("lae", "geb", "yat")) %>%
  filter(response %in% c("GPP", "SIF", "T_TEA"))%>%
  filter(ID2=="Sdeep")
df_p3$time <- as.numeric(df_p3$time)
p <- ggplot(data = df_p3, aes(x=SWC_deep, y=mean)) +
  #geom_pointrange(aes(ymax = pc97.5, ymin = pc2.5, color = site), fatten=.2,alpha=.3) +
  geom_pointrange(aes(ymax = pc97.5, ymin = pc2.5, color = time), fatten=.2,alpha=.3) +
  scale_color_viridis_c(option = "D") +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_grid2(site~response, scales = "free", independent = "all") +
  #facet_wrap("response", scales = "free_y") +
  labs(title = NULL, y = "dY/dTA", x = "Air Temperature (°C)", color = "time")+
  theme_bw() +
  theme(legend.position = "top",
        legend.text=element_text(size=14),
        #legend.title = element_blank(),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle , vjust = 1, hjust=1),
        text = element_text(size=14),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(color = "black", fill = "white", size = 1))
p
ggsave2("p_poster_net_sdeep_time.png", plot = p, path = path_out, width = 9, height = 5)


df_p3 <- df_p %>%
  filter(site=="lnf") %>%
  filter(response %in% c("GPP", "SIF"))%>%
  filter(ID2=="VPD")
df_p3$time <- as.numeric(df_p3$time)
p <- ggplot(data = df_p3, aes(x=SWC_deep, y=mean)) +
  #geom_pointrange(aes(ymax = pc97.5, ymin = pc2.5, color = site), fatten=.2,alpha=.3) +
  geom_pointrange(data = df_p3, aes(ymax = pc97.5, ymin = pc2.5, color = time), fatten=.2,alpha=.3) +
  scale_color_viridis_c(option = "D") +
  geom_hline(yintercept = 0, linetype = 2) +
  #facet_grid(response~ID2, scales = "free") +
  facet_grid2("response", scales = "free", independent = "all") +
  #facet_wrap(site~response, scales = "free_y") +
  labs(title = NULL, y = "dY/dVPD", x = "deep SWC")+
  theme_bw() +
  theme(legend.position = "right",
        legend.text=element_text(size=18),
        #legend.title = element_blank(),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle , vjust = 1, hjust=1),
        text = element_text(size=24),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(color = "black", fill = "white", size = 1))
p

#ggsave2("p_poster_net_vpd_swc.png", plot = p, path = path_out, width = 10, height = 4)


# Take average of net sensitivities by hour
df_p_grouped <- df_p %>%
  group_by(site,cat,response,time,ID2) %>%
  summarise_at(vars(mean, median, pc2.5, pc97.5), list(avg = mean))

df_p3 <- df_p_grouped %>%
  filter(response %in% c("GPP", "SIF", "T_TEA"))#%>%
  #filter(ID2=="VPD")
df_p3$time <- as.numeric(df_p3$time)
p <- ggplot(data = df_p3, aes(x=time, y=mean_avg)) +
  #geom_pointrange(aes(ymax = pc97.5, ymin = pc2.5, color = site), fatten=.2,alpha=.3) +
  geom_pointrange(data = df_p3, aes(ymax = pc97.5_avg, ymin = pc2.5_avg, color = site), fatten=.2) +
  #scale_color_viridis_c(option = "D") +
  scale_color_viridis_d(option = "D") +
  #scale_shape_manual(values = c(1,2)) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_nested(response ~ site + ID2, scales = "free") +
  #facet_grid2(site ~ response, scales = "free", independent = "all") +
  labs(title = NULL, y = "dY/dX", x = "hour of day")+
  theme_bw() +
  theme(legend.position = "right",
        legend.text=element_text(size=18),
        #legend.title = element_blank(),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle , vjust = 1, hjust=1),
        text = element_text(size=24),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(color = "black", fill = "white", size = 1))
p


ggsave2("p_poster_net_vs_time.png", plot = p, path = path_out, width = 15, height = 6)



### yatir heatwave


