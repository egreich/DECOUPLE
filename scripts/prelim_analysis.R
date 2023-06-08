
library(tidyverse)
library(cowplot) # for ggsave2
library(ggpubr) # for stat_cor
library(ggh4x) # for facet_wrap2
library(biwavelet)

# Create necessary folders if they do not already exist
if(!file.exists("plots")) { dir.create("plots")}

path_out = "./plots/" # set save path

key <- c("seg", "ses", "wjs", "mpj", "vcp", "vcm1","vcm2", "vcs")
list_df_daily <- c()

for(i in c(1:8)){
  list_df_daily[[i]]<- read.csv(paste("../ETpart/output_dfs/df_daily_",key[i],".csv",sep=""))
}
df_daily <- bind_rows(list_df_daily)
df_daily$site <- factor(df_daily$site, levels =  c("seg", "ses", "wjs", "mpj", "vcp", "vcm1","vcm2", "vcs"),
                        labels = c("US-Seg", "US-Ses", "US-Wjs", "US-Mpj", "US-Vcp", "US-Vcm1","US-Vcm2", "US-Vcs"))

### some loose sanity cleaning
df_daily$E <- ifelse(df_daily$E > df_daily$ET, df_daily$ET, df_daily$E)
df_daily$E <- ifelse(df_daily$E < 0, 0, df_daily$E)

df_daily$T_check <- df_daily$ET - df_daily$E

reg <- lm(T_check ~ mean_T.pred, data = df_daily)
summary(reg)


p <- df_daily %>%
  ggplot(aes(x = T_check, y= mean_T.pred)) +
  geom_point(alpha=0.5,size=.05)+
  geom_smooth(method="lm", se = F, aes(color = "linear regression line")) +
  geom_abline(slope=1, intercept=0, lty=2, col="blue", size=.8)+
  stat_cor(aes(label = ..rr.label..), color = "red", geom = "label") +
  scale_color_manual(values = c("red", "blue")) + 
  ylim(0,6) + xlim(0,6) +
  facet_wrap2(vars(site), ncol = 2, strip = strip_nested(bleed = FALSE)) +
  labs(title = NULL, x="non-bayesian T", y="DEPART T") +
  theme(legend.position = "top",
        legend.text=element_text(size=12),
        text = element_text(size=14),
        legend.title = element_blank(),
        panel.background = element_rect(fill="white"),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(colour="black"),
        aspect.ratio=1,
        plot.title = element_text(hjust = 0.5))
p

ggsave2("test_fit_facet.png", plot = p, path = path_out, width = 6, height = 8)


p <- df_daily %>%
  ggplot(aes(x = T_check, y= GPP)) +
  geom_point(alpha=0.5,size=.05)+
  geom_smooth(method="lm", se = F, aes(color = "linear regression line")) +
  geom_abline(slope=1, intercept=0, lty=2, col="blue", size=.8)+
  stat_cor(aes(label = ..rr.label..), color = "red", geom = "label") +
  scale_color_manual(values = c("red", "blue")) + 
  ylim(0,6) + xlim(0,6) +
  facet_wrap2(vars(site), ncol = 2, strip = strip_nested(bleed = FALSE)) +
  labs(title = NULL, x="non-bayesian T", y="GPP") +
  theme(legend.position = "top",
        legend.text=element_text(size=12),
        text = element_text(size=14),
        legend.title = element_blank(),
        panel.background = element_rect(fill="white"),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(colour="black"),
        aspect.ratio=1,
        plot.title = element_text(hjust = 0.5))
p

ggsave2("test2_fit_facet.png", plot = p, path = path_out, width = 6, height = 8)





df <- df_daily %>%
  filter(site == "US-Mpj") %>%
  select(T_check,GPP)


t1 = cbind(c(1:nrow(df)), df$T_check)
t2 = cbind(c(1:nrow(df)), df$GPP)

nrands = 10

wtc.AB = wtc(t1, t2, nrands = nrands)

temp.period <- wtc.AB$period # take the longest period
temp.rsq <- as.data.frame(wtc.AB$rsq)
temp.rsq.mat <- as.matrix(wtc.AB$rsq)
temp.phase <- as.data.frame(wtc.AB$phase)
temp.signif <- as.data.frame(wtc.AB$signif)
#test <- as.data.frame(wtc.AB$coi)
for(m in c(1:nrow(temp.rsq))){
  for(n in c(1:ncol(temp.rsq))){
    if(temp.signif[m,n] < 0.7){
      temp.rsq[m,n] = NA
      temp.phase[m,n] = NA
    }
  }
}
temp.avg.rsq <- rowMeans(temp.rsq, na.rm = T)
temp.sd.rsq <- rowSds(wtc.AB$rsq, na.rm = T)
temp.avg.phase <- rowAngleMeans(temp.phase)
for(m in c(1:nrow(temp.rsq))){ # if it's antiphase, make the rsq negative to reflect that
  temp.avg.rsq[m] <- ifelse(temp.avg.phase[m] < 0, temp.avg.rsq[m] * -1, temp.avg.rsq[m])
}
#temp.avg.phase <- calclag(wtc.AB)
list_var_rsq[[k]] <- temp.avg.rsq # will save list of avg rsq
list_sd_rsq[[k]] <- temp.sd.rsq # will save list of avg rsq
list_var_phase[[k]] <- temp.avg.phase # will take avg of phases

# Plotting a graph

  png(paste("./plots/p_", "US-Mpj", "_TvsGPP", ".png",sep=""), width = 900, height = 700)
  
  par(oma = c(0, 0, 0, 1), mar = c(5, 4, 5, 5) + 0.1)
  plot(wtc.AB, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2,
       lwd.sig = 2, arrow.lwd = 0.06, arrow.len = 0.11,
       ylim=range(2:64),
       ylab = "Scale", xlab = "Week",
       plot.cb = TRUE, main = "US-Mpj Wavelet Coherence: T vs GPP")
  
  # Adding grid lines
  n = length(t1[, 1])
  abline(v = seq(52, n, 52), h = 1:16, col = "brown", lty = 1, lwd = 1)
  
  # Defining x labels
  if(i %in% c(5)){ # vcp
    axis(side = 3, at = c(seq(0, n, 52)), labels = c(seq(2009, 2020, 1)))
  }
  if(i %in% c(3,4)){ # wjs, mpj
    axis(side = 3, at = c(seq(0, n, 52)), labels = c(seq(2009, 2021, 1)))
  }
  if(i == 6){ # vcm1
    axis(side = 3, at = c(seq(0, n, 52)), labels = c(2009,2010,2011,2012,2013))
  }
  if(i == 7){ # vcm2
    axis(side = 3, at = c(seq(0, n, 52)), labels = c(2014,2015,2016,2017,2018,2019,2020))
  }
  if(i == 8){ # vcs
    axis(side = 3, at = c(seq(0, n, 52)), labels = c(seq(2017, 2021, 1)))
  }
  if(i %in% c(1,2)){ # seg, ses
    axis(side = 3, at = c(seq(0, n, 52)), labels = c(seq(2008, 2021, 1)))
  }
  
  dev.off()

#################### with geb data
  
  load("./data_clean/gebdat.Rdata")
  
  df <- gebdat
  
  
  t1 = cbind(c(1:nrow(df)), df$Gs)
  t2 = cbind(c(1:nrow(df)), df$GPP)
  
  nrands = 10
  
  wtc.AB = wtc(t1, t2, nrands = nrands)
  
  temp.period <- wtc.AB$period # take the longest period
  temp.rsq <- as.data.frame(wtc.AB$rsq)
  temp.rsq.mat <- as.matrix(wtc.AB$rsq)
  temp.phase <- as.data.frame(wtc.AB$phase)
  temp.signif <- as.data.frame(wtc.AB$signif)
  #test <- as.data.frame(wtc.AB$coi)
  for(m in c(1:nrow(temp.rsq))){
    for(n in c(1:ncol(temp.rsq))){
      if(temp.signif[m,n] < 0.7){
        temp.rsq[m,n] = NA
        temp.phase[m,n] = NA
      }
    }
  }
  temp.avg.rsq <- rowMeans(temp.rsq, na.rm = T)
  temp.sd.rsq <- rowSds(wtc.AB$rsq, na.rm = T)
  temp.avg.phase <- rowAngleMeans(temp.phase)
  for(m in c(1:nrow(temp.rsq))){ # if it's antiphase, make the rsq negative to reflect that
    temp.avg.rsq[m] <- ifelse(temp.avg.phase[m] < 0, temp.avg.rsq[m] * -1, temp.avg.rsq[m])
  }
  #temp.avg.phase <- calclag(wtc.AB)
  list_var_rsq[[k]] <- temp.avg.rsq # will save list of avg rsq
  list_sd_rsq[[k]] <- temp.sd.rsq # will save list of avg rsq
  list_var_phase[[k]] <- temp.avg.phase # will take avg of phases
  
  # Plotting a graph
  
  png(paste("./plots/p_", "DE-Geb", "_GsvsGPP", ".png",sep=""), width = 900, height = 700)
  
  par(oma = c(0, 0, 0, 1), mar = c(5, 4, 5, 5) + 0.1)
  plot(wtc.AB, plot.phase = TRUE, lty.coi = 1, col.coi = "grey", lwd.coi = 2,
       lwd.sig = 2, arrow.lwd = 0.06, arrow.len = 0.11,
       ylab = "Scale", xlab = "Half-hour",
       plot.cb = TRUE, main = "DE-Geb Wavelet Coherence: Gs vs GPP")
  
  # Adding grid lines
  n = length(t1[, 1])
  abline(v = seq(52, n, 52), h = 1:16, col = "brown", lty = 1, lwd = 1)
  
  # Defining x labels
  if(i %in% c(5)){ # vcp
    axis(side = 3, at = c(seq(0, n, 52)), labels = c(seq(2009, 2020, 1)))
  }
  if(i %in% c(3,4)){ # wjs, mpj
    axis(side = 3, at = c(seq(0, n, 52)), labels = c(seq(2009, 2021, 1)))
  }
  if(i == 6){ # vcm1
    axis(side = 3, at = c(seq(0, n, 52)), labels = c(2009,2010,2011,2012,2013))
  }
  if(i == 7){ # vcm2
    axis(side = 3, at = c(seq(0, n, 52)), labels = c(2014,2015,2016,2017,2018,2019,2020))
  }
  if(i == 8){ # vcs
    axis(side = 3, at = c(seq(0, n, 52)), labels = c(seq(2017, 2021, 1)))
  }
  if(i %in% c(1,2)){ # seg, ses
    axis(side = 3, at = c(seq(0, n, 52)), labels = c(seq(2008, 2021, 1)))
  }
  
  dev.off()

  
 p <-  ggplot(data = gebdat) +
    geom_line(aes(x=c(1:nrow(gebdat)), y = Tair)) +
    theme_bw()
  

  






