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

# mod1 daily
sitename = c("lae", "crk", "geb", "lnf", "yat")
varname = c("ET", "GPP", "SIF_O2A", "WUE_GPP", "WUE_SIF", "T_TEA", "WUE_GPP_TEA", "WUE_SIF_TEA")
df_list <- list()
df_list2 <- list()
for(i in 1:length(sitename)){
  for(j in 1:length(varname)){
    
    #dffilename <- paste("./models/model1/", sitename[i],"/coda/df_mod1_", sitename[i], "_", varname[j], ".csv", sep = "")
    dffilename <- paste("./output_model1/df/df_", sitename[i], "_", varname[j],"_daily", ".csv", sep = "")
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
df_daily = bind_rows(df_list2)

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
df_p$ID <-factor(df_p$ID1 , levels = c("1","2","3","4","5", "6", "7", "8", "9", "10"), labels = c("VPD*Tair", "Tair*Sshall", "Tair*Sdeep", "VPD*Sshall", "VPD*Sdeep", "Sshall*Sdeep", "PAR*VPD", "PAR*Tair", "PAR*Sshall", "PAR*Sdeep"))
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
df_p$ID <-factor(df_p$ID1 , levels = c("1","2","3","4","5", "6"), labels = c("VPD*Tair", "Tair*Sshall", "Tair*Sdeep", "VPD*Sshall", "VPD*Sdeep", "Sshall*Sdeep"))
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

sitename <- c("yat", "crk")
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
    labs(title = sitename[i], y = "dY/dX", x = "Half-hour Timestep")+
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
### Read in YINs for dates
# For comparing T_TEA and T.pred
# sitename = c("lae", "geb", "yat")
# responsename = c("WUE_GPP", "WUE_SIF")
# df_list <- list()
# df_list2 <- list()
# for(i in c(1:length(sitename))){
#   
#   # load org data for each site
#   load(paste("./data_clean/",sitename[i],"dat.RData", sep="")) # dat
#   
#   # Load site metadata for variables needed for process-based evaporation model
#   df_soil <- read.csv("./data_misc/soildata.csv")
#   
#   df_siteinfo <- df_soil %>% # filter for the site we're running the model for
#     filter(sitecode == sitename[i])
#   
#   Z = df_siteinfo$ref_height_m # reference height (m)
#   h = 	df_siteinfo$elev_m # elevation
#   fc = df_siteinfo$FC # field capacity
#   fclay = df_siteinfo$clay/100 # clay fraction
#   fsand = df_siteinfo$sand/100 # sand fraction
#   dat$SWC_shall <- dat$SWC_shall/100 # convert to fraction
#   dat$SWC_deep <- dat$SWC_deep/100 # convert to fraction
#   
#   # Calculate intermediate parameters needed for evaporation model
#   # Need dataframe with columns for PA, TA, TS, WS, RH, SWC_shall
#   # Need separate data for Z, fsand, fclay, fc
#   dat <- get_evap(dat, Z = Z, h = h, fc = fc, fclay = fclay, fsand = fsand)
#   
#   for(j in c(1:length(responsename))){
#     df_p <- df2 %>% # df2 is for model 2
#       filter(var %in% c("T.pred")) %>%
#       filter(site == sitename[i]) %>%
#       filter(response == responsename[j])
#     
#     if(nrow(df_p)==0){next} # if we don't have the results, skip
#     
#     
#     # If we are missing env variables we don't want to gapfill 5 timesteps into the past,
#     # make ET NA to skip over those periods
#     # we will filter out the NAs in ET later in our YIN df to skip over those rows in our weight calculations
#     for(k in c(1:nrow(dat))){
#       if(k < 7){
#         next
#       }
#       for(m in c(1:6)){
#         if(is.na(dat$SWC_shall[k-m])){
#           dat[k,"ET"] = NA
#         }
#         if(is.na(dat$SWC_deep[k-m])){
#           dat[k,"ET"] = NA
#         }
#       }
#       
#     }
#     
#     # Get the YIN for the site-variable combo
#     YIN = dat %>%
#       rowid_to_column("ind") %>% # ind: index to link the growing season Y variables back to the appropriate row in the covariate data set
#       mutate(time = as.numeric(format(TIMESTAMP, "%H%M"))) %>%
#       filter(between(time, 730, 1600)) %>% # filter from 7:30 am to 4:00 pm everyday
#       filter(!is.na(SIF_O2A)) %>% # filter any gaps for SIF out that we don't want to fill
#       filter(!is.na(ET)) %>%
#       filter(!is.na(GPP)) %>%
#       filter(!is.na(TA)) %>% # filter any weird start/end gaps for envir variables, this won't cause gaps in the middle of the day, since we already gap-filled
#       filter(!is.na(VPD)) %>%
#       filter(!is.na(SWC_shall)) %>%
#       filter(!is.na(WS)) %>% # filter away any evap equation-dependent variables we are missing or do not want to gapfill
#       filter(!is.na(Ri)) %>%
#       filter(!is.na(ea)) %>%
#       filter(!is.na(es))
#     
#     TIMESTAMP = YIN$TIMESTAMP
#     df_TEA = data.frame(X = NA,
#                         var = "T.pred",
#                         ID1 = c(1:nrow(YIN)),
#                         ID2 = NA,
#                         mean = YIN$T_TEA,
#                         median = NA,
#                         pc2.5 = NA,
#                         pc97.5 = NA,
#                         overlap0 = NA,
#                         gel = NA,
#                         response = "T_TEA",
#                         site = sitename[i],
#                         TIMESTAMP = YIN$TIMESTAMP)
#     df_p <- cbind(df_p, TIMESTAMP)
#     df_p <- rbind(df_p, df_TEA)
#     
#     
#     
#     df_list[[j]] <- df_p
#   }
#   df_list2[[i]] <- bind_rows(df_list)
#   
# }
# df_p <- bind_rows(df_list2)

for(i in 1:nrow(df_p)){
  if(df_p$response[i]=="WUE_GPP"){
    df_p$mean[i] = df_p$mean[i]/48
  }else if(df_p$response[i]=="WUE_SIF"){
    df_p$mean[i] = df_p$mean[i]/48
  }else{
    df_p$mean[i] = df_p$mean[i]
  }
}

p <- df_p %>%
  filter(site=="lae") %>%
  ggplot(aes(x=TIMESTAMP)) +
  geom_point(aes(y=mean, color = response), position = position_dodge(width = 1), fatten = .5) +
  #geom_pointrange(aes(y=mean, ymax = pc97.5, ymin = pc2.5, color = response), position = position_dodge(width = 1), fatten = .5) +
  facet_grid(response ~ site) +
  #xlim(c(102,197)) +
  #ylim(c(0,.4)) +
  labs(title = NULL, y = "T", x = "Half-hour Timestep")+
  theme_bw() +
  theme(legend.position = "top",
        legend.text=element_text(size=14),
        #axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size=14),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.spacing=unit(.05, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        strip.background = element_rect(color = "black", fill = "white", size = 1))
p

ggsave2("p_T_comp_lae.png", plot = p, path = path_out, width = 6, height = 8)


tempsite = "yat"
T_TEA = df_p %>%
  filter(site == tempsite) %>%
  filter(response == "T_TEA") %>%
  select(mean, TIMESTAMP) %>%
  rename(T_TEA = mean)
T_GPP = df_p %>%
  filter(site == tempsite) %>%
  filter(response == "WUE_GPP") %>%
  select(mean, TIMESTAMP) %>%
  rename(T_GPP = mean)
T_SIF = df_p %>%
  filter(site == tempsite) %>%
  filter(response == "WUE_SIF") %>%
  select(mean, TIMESTAMP) %>%
  rename(T_SIF = mean)

df_T <- left_join(T_GPP, T_SIF, by = "TIMESTAMP")
df_T <- left_join(df_T , T_TEA, by = "TIMESTAMP")

test <- lm(T_TEA ~ T_GPP, data = df_T)
summary(test)
test <- lm(T_TEA ~ T_SIF, data = df_T)
summary(test)

# library(ggpubr) # for 1:1 graphs stat_cor
# p1 <- df_T %>%
#   #pivot_wider(id_cols = c(site,response,ID1,TIMESTAMP), names_from = var, values_from = c(mean,median,pc2.5,pc97.5)) %>%
#   ggplot(aes(x = T_TEA, y= T_GPP)) +
#   #geom_point() +
#   #geom_pointrange(aes(ymin=pc2.5_ET.pred, ymax=pc97.5_ET.pred), alpha=0.5)+
#   geom_smooth(method="lm", se = F, color = "red") +
#   geom_abline(slope=1, intercept=0, lty=2, col="blue", size=1.25)+
#   stat_cor(aes(label = ..rr.label..), color = "red", geom = "label") +
#   #stat_cor(aes(label = ..rr.label..), color = "red", label.x = 0.5, size = 3) +
#   #ylim(0,6) + xlim(0,6) +
#   #facet_row("site", strip.position = "top") +
#   labs(title = NULL, x="observed ET", y="predicted ET") +
#   #theme_classic(base_size = 12)+
#   theme(legend.position = "right",
#         legend.text=element_text(size=14),
#         text = element_text(size=14),
#         legend.title = element_blank(),
#         panel.background = element_rect(fill="white"),
#         axis.line = element_line(color = "black"),
#         axis.text.x = element_text(colour="black"),
#         aspect.ratio=1,
#         plot.title = element_text(hjust = 0.5))
# 
# p1

########################## Selected plots for presentations/posters


df_plot <- df %>%
  filter(var == "dYdX") %>%
  filter(response %nin% c("WUE_GPP", "WUE_SIF","WUE_GPP_TEA", "WUE_SIF_TEA", "ET"))%>%
  mutate(mod = 1)
df_plot2 <- df2 %>%
  filter(var == "dYdX") %>%
  mutate(mod = 2)

df_p$TIMESTAMP <- as.character(df_p$TIMESTAMP) # timestamps should be in chatacter format for joining
df_p <- left_join(df_p, df_wue, by = c("DOY", "site", "TIMESTAMP", "response"))
df_p$mean <- df_p$mean * df_p$wue_mean
df_p$pc2.5 <- df_p$pc2.5 * df_p$wue_pc2.5
df_p$pc97.5 <- df_p$pc97.5 * df_p$wue_pc97.5
df_plot <- rbind(df_plot, df_plot2)
df_plot$ID2 <-factor(df_plot$ID2 , levels = c("1","2","3","4", "5"), labels = c("VPD", "Tair", "Sshall", "Sdeep", "PAR"))
df_plot$response <-factor(df_plot$response , levels = c("WUE_GPP", "WUE_SIF","GPP","SIF_O2A", "T_TEA"), labels = c("GPP/T (WUE)", "SIF/T (WUE)", "GPP", "SIF", "T_TEA"))
df_plot$cat <- ifelse(df_plot$ID2 %in% c("VPD", "Tair", "PAR"), "atmospheric driver", "soil moisture driver")
df_p <- df_plot


sitename = c("lae", "crk", "geb", "lnf", "yat")
df_list <- list()
for(i in c(1:length(sitename))){
  print(paste("i: ", i, sep=""))
  # load org data for each site
  load(paste("./data_clean/",sitename[i],"dat.RData", sep="")) # dat
  
  # aggregate by hour
  dat <- dat %>%
    select(-X) %>%
    mutate(hour = format(dat$TIMESTAMP, format ="%H"), DOY = floor(DOY)) %>%
    group_by(DOY,hour) %>%
    summarise(across(everything(), mean, na.rm = T)) %>%
    ungroup()
  
  dat$site <- sitename[i]
  
  df_list[[i]] <- dat
}
df_clean <- bind_rows(df_list)

# combine clean and output, but make sure TIMESTAMPS are in correct format
df_clean$TIMESTAMP <- as.character(df_clean$TIMESTAMP)
df_p$TIMESTAMP <- as.character(df_p$TIMESTAMP) # timestamps should be in chatacter format for joining
df_p <- left_join(df_p, df_clean, by = c("DOY", "site", "TIMESTAMP"))
df_p$TIMESTAMP <- as.POSIXct(df_p$TIMESTAMP)

dtime <- read.table(text=format(df_p$TIMESTAMP, "%H:%M:%S"), sep=":", header=FALSE)
dtime$V2 <- ifelse(dtime$V2==30, .5, 0)# convert minutes to hours
dtime$time <- dtime$V1 + dtime$V2
df_p$time <- dtime$time

# temp
df_wue <- df2 %>%
  filter(var == "WUE.pred") %>%
  select(mean, pc2.5, pc97.5, DOY, TIMESTAMP, response, site)
df_wue$response <-factor(df_wue$response , levels = c("WUE_GPP", "WUE_SIF"), labels = c("GPP/T (WUE)", "SIF/T (WUE)"))
df_wue <-  df_wue %>%
  rename(wue_mean = mean, wue_pc2.5 = pc2.5, wue_pc97.5 = pc97.5)

df_p$TIMESTAMP <- as.character(df_p$TIMESTAMP) # timestamps should be in chatacter format for joining
df_p <- left_join(df_p, df_wue, by = c("DOY", "site", "TIMESTAMP", "response"))
df_p$mean <- ifelse(df_p$response %in% c("GPP/T (WUE)", "SIF/T (WUE)"), df_p$mean * df_p$wue_mean, df_p$mean)
df_p$pc2.5 <- ifelse(df_p$response %in% c("GPP/T (WUE)", "SIF/T (WUE)"), df_p$pc2.5 * df_p$wue_pc2.5, df_p$pc2.5)
df_p$pc97.5 <- ifelse(df_p$response %in% c("GPP/T (WUE)", "SIF/T (WUE)"), df_p$pc97.5 * df_p$wue_pc97.5, df_p$pc97.5)
df_p$TIMESTAMP <- as.POSIXct(df_p$TIMESTAMP)


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

test <- df_p %>%
  filter(site == "yat") %>%
  filter(response %in% c("SIF/T (WUE)")) %>%
  select(wue_mean)

test2 <- df_clean %>%
  filter(site == "yat") %>%
  select(GPP, T_TEA) %>%
  mutate(wue_TEA = GPP/T_TEA)


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
  filter(site=="yat")%>%
  filter(response %in% c("GPP/T (WUE)", "SIF/T (WUE)")) %>%
  ggplot(aes(x=time, y=median_avg, color= time)) +
  geom_point() +
  #geom_pointrange(aes(ymax = pc97.5_avg, ymin = pc2.5_avg), color= "black", fatten = .7) +
  #geom_pointrange(aes(ymax = pc97.5_avg, ymin = pc2.5_avg, color= time), fatten = .2) +
  geom_hline(yintercept = 0, linetype = 2) +
  #scale_color_viridis_c(option = "D") +
  scale_color_distiller(palette = "YlOrBr") +
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
ggsave2("p_poster_net_yat3.png", plot = p, path = path_out, width = 9, height = 5.5)


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
  filter(site=="yat")%>%
  filter(response %in% c("GPP", "SIF")) %>%
  ggplot(aes(x=time, y=mean_avg, color = time)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = 2) +
  #scale_color_viridis_c(option = "D") +
  scale_color_distiller(palette = "YlOrBr") +
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

ggsave2("p_poster_net_yat_phot.png", plot = p, path = path_out, width = 9, height = 5.5)

p <- df_p %>%
  filter(site=="yat") %>%
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

ggsave2("p_poster_net_yat.png", plot = p, path = path_out, width = 7, height = 4)


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


# net sens tair

df_p1 <- df_p %>%
  filter(site=="crk") %>%
  filter(response %in% c("GPP", "SIF"))%>%
  filter(ID2=="Tair")
df_p2 <- df_p %>%
  filter(site=="yat") %>%
  filter(response %in% c("GPP", "SIF"))%>%
  filter(ID2=="Tair")

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
p <- ggplot(data = df_p3, aes(x=TA, y=mean)) +
  #geom_pointrange(aes(ymax = pc97.5, ymin = pc2.5, color = site), fatten=.2,alpha=.3) +
  geom_pointrange(data = df_p3 %>% filter(site=="crk"), aes(ymax = pc97.5, ymin = pc2.5, color = time, shape = site), fatten=.2,alpha=.3) +
  scale_color_viridis_c(option = "D") +
  #scale_color_distiller(palette = "RdPu") +
  new_scale_color() +
  geom_pointrange(data = df_p3 %>% filter(site=="yat"), aes(ymax = pc97.5, ymin = pc2.5, color = time, shape = site), fatten=.2,alpha=.3) +
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
p <- ggplot(data = df_p3, aes(x=TA, y=mean)) +
  #geom_pointrange(aes(ymax = pc97.5, ymin = pc2.5, color = site), fatten=.2,alpha=.3) +
  geom_pointrange(data = df_p3 %>% filter(site=="lae"), aes(ymax = pc97.5, ymin = pc2.5, color = time, shape = site), fatten=.2,alpha=.3) +
  scale_color_viridis_c(option = "D") +
  #scale_color_distiller(palette = "RdPu") +
  new_scale_color() +
  geom_pointrange(data = df_p3 %>% filter(site=="yat"), aes(ymax = pc97.5, ymin = pc2.5, color = time, shape = site), fatten=.2,alpha=.3) +
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
  filter(site %in% c("crk", "yat")) %>%
  filter(response %in% c("GPP", "SIF"))%>%
  filter(ID2=="Tair")
p <- ggplot(data = df_p3, aes(x=PAR, y=mean)) +
  #geom_pointrange(aes(ymax = pc97.5, ymin = pc2.5, color = site), fatten=.2,alpha=.3) +
  geom_pointrange(data = df_p3 %>% filter(site=="crk"), aes(ymax = pc97.5, ymin = pc2.5, color = TA, shape = site), fatten=.2,alpha=.3) +
  scale_color_viridis_c(option = "D") +
  #scale_color_distiller(palette = "RdPu") +
  new_scale_color() +
  geom_pointrange(data = df_p3 %>% filter(site=="yat"), aes(ymax = pc97.5, ymin = pc2.5, color = TA, shape = site), fatten=.2,alpha=.3) +
  scale_color_distiller(palette = "YlOrBr") +
  scale_shape_manual(values = c(1,2)) +
  geom_hline(yintercept = 0, linetype = 2) +
  #facet_grid(response~ID2, scales = "free") +
  facet_grid2(site ~ response, scales = "free", independent = "all") +
  #facet_wrap(site~response, scales = "free_y") +
  labs(title = NULL, y = "dY/dTA", x = "PAR")+
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

ggsave2("p_poster_net_tair_par2.png", plot = p, path = path_out, width = 10, height = 4)


df_p3 <- df_p %>%
  filter(site %in% c("crk", "yat")) %>%
  filter(response %in% c("GPP", "SIF"))%>%
  filter(ID2=="VPD")
p <- ggplot(data = df_p3, aes(x=SWC_shall, y=mean)) +
  #geom_pointrange(aes(ymax = pc97.5, ymin = pc2.5, color = site), fatten=.2,alpha=.3) +
  geom_pointrange(data = df_p3 %>% filter(site=="crk"), aes(ymax = pc97.5, ymin = pc2.5, color = time, shape = site), fatten=.2,alpha=.3) +
  scale_color_viridis_c(option = "D") +
  #scale_color_distiller(palette = "RdPu") +
  new_scale_color() +
  geom_pointrange(data = df_p3 %>% filter(site=="yat"), aes(ymax = pc97.5, ymin = pc2.5, color = time, shape = site), fatten=.2,alpha=.3) +
  scale_color_distiller(palette = "YlOrBr") +
  scale_shape_manual(values = c(1,2)) +
  geom_hline(yintercept = 0, linetype = 2) +
  #facet_grid(response~ID2, scales = "free") +
  facet_grid2(site ~ response, scales = "free", independent = "all") +
  #facet_wrap(site~response, scales = "free_y") +
  labs(title = NULL, y = "dY/dVPD", x = "shallow SWC")+
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

ggsave2("p_poster_net_vpd.png", plot = p, path = path_out, width = 10, height = 4)


# Take average of net sensitivities by temp
df_p_grouped <- df_p %>%
  group_by(site,cat,response,time,ID2) %>%
  summarise_at(vars(mean, median, pc2.5, pc97.5), list(avg = mean))

df_p3 <- df_p_grouped %>%
  filter(response %in% c("GPP", "SIF", "T_TEA"))#%>%
  #filter(ID2=="VPD")
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


