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

sitename = c("lae", "crk", "lnf", "yat")
varname = c("ET", "GPP", "SIF_O2A", "WUE_GPP", "WUE_SIF")
df_list <- list()
df_list2 <- list()
for(i in 1:length(sitename)){
  for(j in 1:length(varname)){
    
    dffilename <- paste("./models/model1/", sitename[i],"/coda/df_mod1_", sitename[i], "_", varname[j], ".csv", sep = "")
    df_temp <- tryCatch(read.csv(dffilename), error=function(err) NA)
    
    if(is.na(df_temp)){
      next
    }
    
    df_temp <- df_temp %>%
      mutate(var = varname[j], site = sitename[i])
    
    df_list[[j]] = df_temp
  }
  df_list2[[i]] = bind_rows(df_list)
}
df = bind_rows(df_list2)



p <- df %>%
  filter(param %in% c("wSs","wT","wV")) %>%
  ggplot(aes(x=ID, y=mean)) + 
  geom_pointrange(aes(ymax = pc97.5, ymin = pc2.5), position = position_dodge(width = 1), fatten = .5) +
  #geom_point(pch=21, size = .5) +
  scale_color_brewer(palette = "Dark2") +
  #facet_grid(param~var) +
  facet_nested(site ~ var + param, scales="free_x") + 
  #geom_text(data = df_p, aes(label=w_letters), position = position_dodge(width = .7), vjust = -3, size=3) +
  labs(title = NULL, y = "weight", x = NULL)+
  theme_bw() +
  theme(legend.position = "top",
        legend.text=element_text(size=12),
        text = element_text(size=12),
        axis.text.x = element_text(size = 11, angle = 90, vjust = 0.6, hjust = 1),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))
p

ggsave2("p_weights.png", plot = p, path = path_out, width = 8, height = 6)

df_p <- df %>%
  filter(param == "beta1")
df_p$ID <-factor(df_p$ID , levels = c("1","2","3","4","5"), labels = c("VPD", "Tair","PAR", "Sshall", "Sdeep"))
p <- df_p %>%
  filter(ID %in% c("VPD", "Tair","PAR", "Sshall", "Sdeep")) %>%
  ggplot(aes(x=mean, y=ID)) + 
  geom_pointrange(aes(xmax = pc97.5, xmin = pc2.5, alpha = ifelse(pc2.5 <= 0 & pc97.5 >= 0, .9, 1), color = ifelse(pc2.5 < 0, "blue", "red")), position = position_dodge(width = .7), fatten = .5) +
  geom_vline(xintercept = 0, linetype = 2) + 
  #scale_color_manual(values=c("red", "tomato2", "sienna2", "orange","light blue","blue", "dark blue")) +
  #facet_grid(var ~ site, scales = "free_x") +
  #facet_row("var", scales = "free") + 
  facet_nested(site ~ var + param, scales="free_x") +
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
  filter(param == "beta2")
# key: 1 VPD # 2 Tair # 3 PAR # 4 Sshall # 5 Sdeep
# X1a = cbind(X1[,1]^2, X1[,2]^2)
# X2 = X1[,1]*X1[,2], X1[,1]*X1[,4], X1[,1]*X1[,5], X1[,2]*X1[,4], X1[,2]*X1[,5], X1[,4]*X1[,5]
df_p$ID <-factor(df_p$ID , levels = c("1","2","3","4","5", "6"), labels = c("VPD*Tair", "Tair*Sshall", "Tair*Sdeep", "VPD*Sshall", "VPD*Sdeep", "Sshall*Sdeep"))
p <- df_p %>%
  ggplot(aes(x=mean, y=ID)) + 
  geom_pointrange(aes(xmax = pc97.5, xmin = pc2.5, alpha = ifelse(pc2.5 <= 0 & pc97.5 >= 0, .9, 1)), position = position_dodge(width = .7), fatten = .5) +
  geom_vline(xintercept = 0, linetype = 2) + 
  #scale_color_manual(values=c("red", "tomato2", "sienna2", "orange","light blue","blue", "dark blue")) +
  #facet_grid(var ~ site, scales = "free_x") +
  #facet_row("var", scales = "free") + 
  facet_nested(site ~ var + param, scales="free_x") +
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


### Read in YINs for order

# correctly name the net sensitivities (need YIN)
sitename = c("lae", "crk", "lnf", "yat")
varname = c("ET", "GPP", "SIF_O2A", "WUE_GPP", "WUE_SIF")
df_list <- list()
df_list2 <- list()
for(i in 1:length(sitename)){
  
  # load org data for each site
  load(paste("./data_clean/",sitename[i],"dat.RData", sep=""))
  
  for(j in 1:length(varname)){
    df_p <- df %>%
      filter(param %in% c("dYdX")) %>%
      filter(site == sitename[i]) %>%
      filter(var == varname[j])
    
    if(nrow(df_p)==0){next} # if we don't have the results, skip
    
    df_p$ID1 <- c(1:nrow(df_p))
    df_p$ID2 <- NA
    df_p$coupled <- NA
    
    # Get the YIN for the site-variable combo
    YIN = dat %>%
      rowid_to_column("ind") %>% # ind: index to link the growing season Y variables back to the appropriate row in the covariate data set
      mutate(time = as.numeric(format(TIMESTAMP, "%H%M"))) %>%
      filter(between(time, 730, 1600)) %>% # filter from 7:30 am to 4:00 pm everyday
      filter(!is.na(SIF_O2A)) %>% # filter any starting/ending gaps for SIF out that we don't want to fill
      filter(!is.na(ET)) %>%
      filter(!is.na(GPP)) %>%
      filter(!is.na(TA)) %>% # filter any weird start/end gaps for envir variables, this won't cause gaps in the middle of the day, since we already gap-filled
      filter(!is.na(VPD)) %>%
      filter(!is.na(PAR)) %>%
      filter(!is.na(SWC_shall)) %>%
      filter(!is.na(SWC_deep)) %>%
      select(ind, all_of(varname[j]), rsquared, coupled) # select variables we will be using as Y variables
    
    start <- 1
    end <- nrow(YIN)
    df_p$ID2[start:end] <- "VPD"
    df_p$coupled[start:end] <- YIN$coupled[start:end]
    df_p$rsquared[start:end] <- YIN$rsquared[start:end]
    start <- end+1
    end <- start + nrow(YIN) - 1
    df_p$ID2[start:end] <- "Tair"
    start <- end+1
    end <- start + nrow(YIN) - 1
    df_p$ID2[start:end] <- "PAR"
    start <- end+1
    end <- start + nrow(YIN) - 1
    df_p$ID2[start:end] <- "Sshall"
    start <- end+1
    end <- start + nrow(YIN) - 1
    df_p$ID2[start:end] <- "Sdeep"
    
    df_list[[j]] <- df_p
  }
  df_list2[[i]] <- bind_rows(df_list)
  
}
df_p <- bind_rows(df_list2)


paramnames <- c("VPD", "Tair", "PAR", "Sshall", "Sdeep")
df_list <- list()
df_list2 <- list()
for(i in 1:length(varname)){
  df_param1 <- df_p %>%
    filter(var == varname[i])
  for(j in 1:length(paramnames)){
    df_param2 <- df_param1 %>%
      filter(ID2 == paramnames[j])

    df_param2$ID1 <- c(1:nrow(df_param2))
    df_list[[j]] <- df_param2
  }
  df_list2[[i]] <- bind_rows(df_list)
}
df_p <- bind_rows(df_list2)

p <- df_p %>%
  filter(ID2 != "PAR") %>%
  #filter(var %nin% "WUE_SIF") %>%
  filter(param == "dYdX") %>%
  ggplot(aes(x=ID1, y=mean)) +
  #ggplot(aes(x=ID, y=mean)) +
  geom_pointrange(aes(ymax = pc97.5, ymin = pc2.5, color = ifelse(pc2.5 <= 0 & pc97.5 >= 0, "nonsignificant", "significant")), position = position_dodge(width = 1), fatten = .5) +
  scale_color_manual(values = c("gray", "black")) +
  geom_hline(yintercept = 0, linetype = 2) +
  #geom_rect(aes(xmin = date, xmax = dplyr::lead(date), ymin = -Inf, ymax = Inf, fill = Season), alpha = 0.1) +
  #scale_fill_brewer(palette = "Dark2") +
  #facet_grid(var~ID2, scales = "free") +
  facet_nested_wrap(ID2 ~ var + site, scales="free", nrow=4) +
  labs(title = NULL, y = "dY/dX", x = "Half-hour Timestep")+
  theme_bw() +
  theme(legend.position = "top",
        legend.text=element_text(size=10),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size=10),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))
p

ggsave2("p_netsens.png", plot = p, path = path_out, width = 14, height = 8)


p <- df_p %>%
  filter(ID2 != "PAR") %>%
  #filter(var %nin% "WUE_SIF") %>%
  filter(param == "dYdX") %>%
  ggplot(aes(x=ID1, y=mean)) +
  #ggplot(aes(x=ID, y=mean)) +
  geom_pointrange(aes(ymax = pc97.5, ymin = pc2.5, color = rsquared), position = position_dodge(width = 1), fatten = .5) +
  #scale_color_manual(values = c("gray", "black")) +
  geom_hline(yintercept = 0, linetype = 2) +
  #geom_rect(aes(xmin = ID1, xmax = max(ID1), ymin = -Inf, ymax = Inf, fill = rsquared), alpha = 0.1) +
  #scale_fill_brewer(palette = "Dark2") +
  #facet_grid(var~ID2, scales = "free") +
  facet_nested_wrap(ID2 ~ var + site, scales="free", nrow=4) +
  labs(title = NULL, y = "dY/dX", x = "Half-hour Timestep")+
  theme_bw() +
  theme(legend.position = "top",
        legend.text=element_text(size=10),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size=10),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))
p

ggsave2("p_netsens2.png", plot = p, path = path_out, width = 14, height = 8)





