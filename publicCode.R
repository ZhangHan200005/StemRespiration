##---------
#*2024-12-13
#*by ZHANG Han
#*descriptions:code used in analysis of 'Thermal acclimation of stem respiration implies a weaker carbon-climate feedback'

# package
{
  library(grid)
  library(data.table)
  library(dplyr)
  library(lme4)
  library(MuMIn)
  library(eivtools)
  library(ggplot2)
  library(ggpubr)
  library(ggpmisc)
  library(ggsci)
  library(rnaturalearthdata)
  library(rnaturalearth)
  library(s2)
  library(visreg)
  library(cowplot)
  library(car)
  library(purrr)
  library(progress)
  library(raster)
  library(ncdf4)
  library(tidyverse)
  library(sf)
  library(RColorBrewer)
  library(reshape2)
  library(caret)
  library(raster)
  library(ncdf4)
  library(stringr)
  library(grid)
  library(gridExtra)
  library(egg)
  library(plotrix)
}

# READ DATA
setwd("G:/proj_StemRespiration/Manuscript/revision2409/submission1215/Code")
data_global <- read.csv("GSRD_global.csv",sep=";")
data_seasonal <- read.csv("GSRD_seasonal.csv",sep=";")
w_data <- read.csv("GSRD_warmingExp.csv",sep=";")
# other data URL:https://github.com/ZhangHan200005/StemRespiration.git

# Figure1
{
  data_global$LNRSGT <- log(data_global$RS_MASS25*2.2^((data_global$T5-25)/10))
  #Figure 1A: lnrs25 ~ T5 (sapwood mass-based stem respiration to temperature)
  {
    #------------boxplot------------
    data_box <- cut(data_global$T5,breaks=c(0,5,10,15,20,25,30),
                    labels =c(2.5,7.5,12.5,17.5,22.5,27.5))
    data_box <- cbind(data_global,data_box)
    
    figure_boxplotA <- ggplot(data_box, aes(x=data_box, y=LNRS25))+ 
      geom_abline(intercept=0.7, slope=-0.39, colour="#A8021B", linetype=1,size=1)+ #just as a sketch
      geom_boxplot(width=0.5,notchwidth=0.8,lwd=0.5,
                   outlier.shape=NA)+
      # Style
      theme_bw()+
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
            panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
            legend.position = 'none')+
      theme(text = element_text(family="A"))+
      # Axis
      theme(axis.text.x = element_text(size = 10,color=rgb(0.2,0.2,0.2),family="ARL"),
            axis.text.y = element_text(size = 10,color=rgb(0.2,0.2,0.2),family="ARL"),
            axis.line.x = element_line(linetype=1,color="black",size=1))+
      xlab(NULL)+
      ylab(NULL)+
      ylim(-4.9,5)
    
    
    #------------scatter plot------------
    figure_scatterA <- ggplot(data_global,aes(x=T5,y=LNRS25))+
      # Style
      theme_bw()+
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
            panel.border = element_rect(fill=NA,color="black", linewidth =1.3, linetype="solid"))+
      theme(text = element_text(family ="HEL"))+
      # Axis
      labs(x=expression(paste(italic(T[5])," (°C)")),
           y=expression(paste("ln ",italic(r[s25]),"  (nmol ",CO[2]," ",g^-1, s^-1,")")))+
      theme(axis.title.x = element_text(vjust = 0, size=17,face="bold"),
            axis.title.y = element_text(vjust = 0,size=17,face="bold"))+
      theme(axis.text.x = element_text(size = 12,color=rgb(0.2,0.2,0.2)),
            axis.text.y = element_text(size = 12,color=rgb(0.2,0.2,0.2)),
            axis.line.x = element_line(linetype=1,color="black",linewidth=1))+
      ylim(-4.9,5)+
      xlim(0,25)+
      # Line
      stat_smooth(aes(color=MEASURMENT,fill=MEASURMENT),method = "lm",formula = y~x,se=T,linewidth=1.5,lty=1)+
      scale_fill_manual(values=c("#002FA7","#FF8000"))+
      guides(color="none",fill="none")+
      
      # Ponit
      geom_point(aes(color=MEASURMENT),shape=16,size=3.2,alpha=0.1)+
      scale_color_manual(values= c("#002FA7","#FF8000"),
                         name=" ",
                         labels = c("In Field","In lab"))+
      # Legend
      theme(legend.position = c(0.15,0.15),
            legend.background = element_rect(fill="transparent"),
            legend.title=element_text(color="black",size=15),
            legend.text = element_text(color="black",size = 14),
            legend.key = element_blank(),
            legend.key.size = unit(1.1, "cm"))+
      guides(
        color = guide_legend(
          direction = "vertical", color = "white",
          override.aes = list(fill = NA)  # 设置置信区间的填充为透明
        )
      )+
      annotate("text", x = 0.5, y = 4.8, label = "A",size=7,fontface = "bold")+
      geom_segment(aes(x = 0.5, xend =24,y = 0.6, yend = (-0.1*(24-0.5)+0.6)),
                   color = "#A8021B", linewidth = 1.5,lty=1)
    
    vie<-viewport(width=0.385,height=0.363,x=0.78,y=0.8)
    figure_scatterA
    print(figure_boxplotA,vp=vie)
    
  }
  #Figure 1B: lnrs.gt ~ T5
  {
    #------------boxplot------------ 
    figure_boxplotB <- ggplot(data_box, aes(x=data_box, y=LNRSGT))+ 
      geom_abline(intercept=-1.2, slope=-0.06, colour="#A8021B", linetype=1,size=1)+  #just as a sketch
      geom_abline(intercept=-1.2, slope=0.23, colour="grey", linetype=1,size=1)+
      geom_boxplot(width=0.5,notchwidth=0.8,lwd=0.5,
                   #outlier.size=0.5,
                   #outlier.fill="white"
                   outlier.shape = NA)+
      #scale_color_manual(values = c("#002FA7","#FF8000"))+
      # Style
      theme_bw()+
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
            panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
            legend.position = 'none')+
      theme(text = element_text(family="ARL"))+
      
      # Axis
      theme(axis.text.x = element_text(size = 10,color=rgb(0.2,0.2,0.2)),
            axis.text.y = element_text(size = 10,color=rgb(0.2,0.2,0.2)),
            axis.line.x = element_line(linetype=1,color="black",size=1))+
      xlab(NULL)+
      ylab(NULL)+
      ylim(-4.9,5)
    
    #------------scatter plot------------
    figure_scatterB <- ggplot(data_global,aes(x=T5,y=LNRSGT))+
      # Style
      theme_bw()+
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
            panel.border = element_rect(fill=NA,color="black", linewidth =1.3, linetype="solid"))+
      theme(text = element_text(family ="HEL"))+
      
      # Axis
      labs(x=expression(paste(italic(T[5])," (°C)")),
           y=expression(paste("ln ",italic(r[s.gt]),"  (nmol ",CO[2]," ",g^-1, s^-1,")")))+
      theme(axis.title.x = element_text(vjust = 0, size=17,face="bold"),
            axis.title.y = element_text(vjust = 0,size=17,face="bold"))+
      theme(axis.text.x = element_text(size = 12,color=rgb(0.2,0.2,0.2)),
            axis.text.y = element_text(size = 12,color=rgb(0.2,0.2,0.2)),
            axis.line.x = element_line(linetype=1,color="black",linewidth=1))+
      ylim(-4.9,5)+
      xlim(0,25)+
      # Line
      stat_smooth(aes(color=MEASURMENT,fill=MEASURMENT),method = "lm",formula = y~x,se=T,linewidth=1.5,lty=1)+
      scale_fill_manual(values=c("#002FA7","#FF8000"))+
      guides(color="none",fill="none")+
      
      # Ponit
      geom_point(aes(color=MEASURMENT),shape=16,size=3.2,alpha=0.1)+
      scale_color_manual(values= c("#002FA7","#FF8000"),
                         name=" ",
                         labels = c("In Field","In lab"))+
      # Legend
      theme(legend.position = 'none')+
      annotate("text", x = 0.5, y = 4.8, label = "B",size=7,fontface = "bold")+
      
      geom_segment(aes(x = 0.5, xend =24,y = -1.2, yend = (-0.023*(24-0.5)-1.2)),
                   color = "#A8021B", size = 2,lty=1)+
      geom_segment(aes(x = 0.5, xend =24,y = -1.2, yend = (0.079*(24-0.5)-1.2)),
                   color = "grey", size = 1.5,lty=1)
    
    vie<-viewport(width=0.385,height=0.363,x=0.78,y=0.8)
    figure_scatterB
    print(figure_boxplotB,vp=vie)
  }
}

# Table 1
{
  #EIV regression analysis
  reliability <- 0.96
  names(reliability) <- 'T5'
  eiv.rs25 <- eivreg(LNRS25~T5,data=data_global,reliability = reliability)
  summary(eiv.rs25)
  confint(eiv.rs25)
  
  eiv.rsgt <- eivreg(LNRSGT~T5,data=data_global,reliability = reliability)
  summary(eiv.rsgt)
  confint(eiv.rsgt)
  
  #In field/ In lab data
  exsitu <- data.table(data_global[which(data_global$MEASURMENT=="LAB"),])
  insitu <- data.table(data_global[which(data_global$MEASURMENT=="FIELD"),])
  
  eiv.rs25 <- eivreg(LNRS25~T5,data=exsitu,reliability = reliability)
  summary(eiv.rs25)
  confint(eiv.rs25)
  
  eiv.rs25 <- eivreg(LNRS25~T5,data=insitu,reliability = reliability)
  summary(eiv.rs25)
  confint(eiv.rs25)
  
  eiv.rsgt <- eivreg(LNRSGT~T5,data=exsitu,reliability = reliability)
  summary(eiv.rsgt)
  confint(eiv.rsgt)
  
  eiv.rsgt <- eivreg(LNRSGT~T5,data=insitu,reliability = reliability)
  summary(eiv.rsgt)
  confint(eiv.rsgt)
}

# fig.S2-4 and table S3-6,S9
{
  #@fig.S2
  {
    data_site <- data_global %>%  
      group_by(LAT,LON,MEASURMENT) %>%  
      summarise(lnrs25_sitemean = log(mean(exp(LNRS25), na.rm = TRUE)),
                T5_sitemean = mean(T5, na.rm = TRUE),
                lnrsgt_sitemean = log(mean(exp(LNRSGT), na.rm = TRUE)),
                .groups = 'drop')  
    
    #fig.S2A: lnrs25 ~ T5(site mean,sapwood mass-based)
    ggplot(data_site,aes(x=T5_sitemean,y=lnrs25_sitemean))+
      # Style
      theme_bw()+
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
            panel.border = element_rect(fill=NA,color="black", size=1.3, linetype="solid"))+
      theme(text = element_text(family ="A",face = "bold"))+
      
      # Ponit
      geom_point(shape="circle",size=5,color=rgb(247/256,154/256,89/256,0.7))+
      
      # Line
      geom_smooth(method = "lm",formula = y~x,se=T,linetype=2,size=1.5,color="black",
                  fill=rgb(247/256,154/256,89/256,0.2))+
      
      # Axis
      labs(x=expression(paste(italic(T[5])," (°C)")),
           y=expression(paste("ln ",italic(r[s25]),"  (nmol ",CO[2]," ",g^-1, s^-1,")")))+
      theme(axis.title.x = element_text(vjust = 0, size=17,face="bold"),
            axis.title.y = element_text(vjust = 0,size=17,face="bold"))+
      theme(axis.text.x = element_text(size = 12,color=rgb(0.2,0.2,0.2)),
            axis.text.y = element_text(size = 12,color=rgb(0.2,0.2,0.2)))+
      ylim(-4,4)+
      xlim(0,25)+
      annotate("text", x = 0.5, y = 3.9, label = "A",size=7,fontface = "bold")+
      geom_segment(aes(x = 1, xend =22,y = 0.7, yend = (-0.1*(22-1)+0.7)),
                   color = "#A8021B", linewidth = 1.5,lty=1)
    
    #fig.S2B:lnrsgt ~ T5 (site-mean,sapwood mass)
    ggplot(data_site,aes(x=T5_sitemean,y=lnrsgt_sitemean))+
      # Style
      theme_bw()+
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
            panel.border = element_rect(fill=NA,color="black", size=1.3, linetype="solid"))+
      theme(text = element_text(family ="A",face = "bold"))+
      
      # Ponit
      geom_point(shape="circle",size=5,color=rgb(247/256,154/256,89/256,0.7))+
      
      # Line
      geom_smooth(method = "lm",formula = y~x,se=T,linetype=2,size=1.5,color="darkgray",
                  fill=rgb(247/256,154/256,89/256,0.2))+
      
      # Axis
      labs(x=expression(paste(italic(T[5])," (°C)")),
           y=expression(paste("ln ",italic(r[s.gt]),"  (nmol ",CO[2]," ",g^-1, s^-1,")")))+
      theme(axis.title.x = element_text(vjust = 0, size=17,face="bold"),
            axis.title.y = element_text(vjust = 0,size=17,face="bold"))+
      theme(axis.text.x = element_text(size = 12,color=rgb(0.2,0.2,0.2)),
            axis.text.y = element_text(size = 12,color=rgb(0.2,0.2,0.2)))+
      ylim(-4,4)+
      xlim(0,25)+
      annotate("text", x = 0.5, y = 3.9, label = "B",size=7,fontface = "bold")+
      geom_segment(aes(x = 1, xend =22,y = -1.1, yend = (-0.0226*21-1.1)),
                   color = "#A8021B", linewidth = 1.5,lty=1)
  }
  
  #@fig.S3
  {
    data_tem <- data.table(data_global[which((data_global$MTEM>=24)&(data_global$MTEM<=26)),])
    data_tem$lnrs.mtem <- log(data_tem$RS_MASS)
    
    #fig.S3A:individual result
    ggplot(data_tem,aes(x=T5,y=lnrs.mtem))+
      # Ponit
      geom_point(shape="circle",size=5,color="#CC6CE7",alpha=0.5)+
      
      # Line
      geom_smooth(method = "lm",formula = y~x,se=T,linetype=2,size=1.1,color="black",fill="#CC6CE7")+
      
      # Style
      theme_bw()+
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
            panel.border = element_rect(fill=NA,color="black", linewidth =1.3, linetype="solid"))+
      theme(text = element_text(family ="HEL"))+
      
      # Axis
      labs(x=expression(paste(italic(T[5])," (°C)")),
           y=expression(paste("ln ",italic(r[s]),"  (nmol ",CO[2]," ",g^-1, s^-1,")")))+
      theme(axis.title.x = element_text(vjust = 0, size=17,face="bold"),
            axis.title.y = element_text(vjust = 0,size=17,face="bold"))+
      theme(axis.text.x = element_text(size = 12,color=rgb(0.2,0.2,0.2)),
            axis.text.y = element_text(size = 12,color=rgb(0.2,0.2,0.2)))+
      ylim(-4,2)+
      xlim(5,25)+
      
      # Legend
      theme(legend.position = c(0.2,0.15),
            legend.background = element_rect(fill="transparent"),
            legend.title=element_text(color="black",size=15),
            legend.text = element_text(color="black",size = 14),
            legend.key = element_blank(),
            legend.key.size = unit(1.1, "cm"))+
      guides(
        color = guide_legend(
          direction = "vertical", color = "white",
          override.aes = list(fill = NA)  # 设置置信区间的填充为透明
        )
      )+
      annotate("text", x = 5.5, y = 2, label = "A",size=7,fontface = "bold")+
      geom_segment(aes(x = 6, xend =22,y = 0.05, yend = (-0.1*(22-6)+0.05)),
                   color = "#A8021B", linewidth = 1.5,lty=1)
    
    #fig.S3B: group by mgdd5 bins
    bin <- 0:25
    bin_mean <- array(dim=c(25,3))
    colnames(bin_mean) <- c("T5","rsmass","lnrs")
    for(i in 1:25){
      data_tem_bin <- subset(data_tem, T5 > bin[i] & T5 < bin[i+1])
      bin_mean[i,1] <- ((bin[i]+bin[i+1])/2)
      bin_mean[i,2] <- mean(data_tem_bin$RS_MASS)
    }
    bin_mean[,3] <- log(bin_mean[,2])
    #plot(bin_mean[,1],bin_mean[,3])
    bin_mean <- data.table(bin_mean)
    bin_mean <- na.omit(bin_mean)
    
    ggplot(bin_mean,aes(x=T5,y=lnrs))+
      # Style
      theme_bw()+
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
            panel.border = element_rect(fill=NA,color="black", size=1.3, linetype="solid"))+
      theme(text = element_text(family ="A",face = "bold"))+
      
      # Ponit
      geom_point(shape="circle",size=5,color="black")+
      
      # Line
      geom_smooth(method = "lm",formula = y~x,se=T,linetype=2,size=1.1,color="black",fill="grey",alpha=0.3)+
      
      # Axis
      labs(x=expression(paste(italic(T[5])," (°C)")),
           y=expression(paste(ln," ",italic(r[s]),"  (nmol ",CO[2]," ",g^-1, s^-1,")")))+
      theme(axis.title.x = element_text(vjust = 0,size=16),
            axis.title.y = element_text(vjust = 0,size=16))+
      theme(axis.text.x = element_text(size = 12,color="black"),
            axis.text.y = element_text(size = 12,color="black"))+
      guides(color="none",fill="none")+
      geom_segment(aes(x = 6, xend =22,y = 0.38, yend = (-0.1*(22-6)+0.38)),
                   color = "#A8021B", linewidth = 1.5,lty=1)+
      annotate("text", x = 5.5, y = 2, label = "B",size=7,fontface = "bold")+
      ylim(-4,2)+
      xlim(5,25)
  }
  
  #@fig.S4
  {
    pinus <- data.table(data_global[which(data_global$GENUS=="Pinus"),])
    
    ggplot(pinus,aes(x=T5,y=LNRS25,color=FINAL_SPNAME))+
      # Ponit
      geom_point(shape="circle",size=5,alpha=0.4)+
      
      # Line
      geom_smooth(method = "lm",formula = y~x,se=T,linetype=2,size=1.5,color="black")+
      
      # Style
      theme_bw()+
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
            panel.border = element_rect(fill=NA,color="black", linewidth =1.3, linetype="solid"))+
      theme(text = element_text(family ="HEL"))+
      
      # Axis
      labs(x=expression(paste(italic(T[5])," (°C)")),
           y=expression(paste("ln ",italic(r[s25]),"  (nmol ",CO[2]," ",g^-1, s^-1,")")))+
      theme(axis.title.x = element_text(vjust = 0, size=17,face="bold"),
            axis.title.y = element_text(vjust = 0,size=17,face="bold"))+
      theme(axis.text.x = element_text(size = 12,color=rgb(0.2,0.2,0.2)),
            axis.text.y = element_text(size = 12,color=rgb(0.2,0.2,0.2)),
            legend.position = "right")+
      # ylim(-2.5,4.3)+
      xlim(0,18)+
      
      # Legend
      theme(#legend.position = c(0.2,0.15),
        #legend.background = element_rect(fill="transparent"),
        legend.title=element_text(color="black",size=12,family="HEL"),
        legend.text = element_text(color="black",size = 10),
        #legend.key = element_blank(),
        #legend.key.size = unit(1.1, "cm")
      )+
      scale_colour_discrete(name="Species")+
      guides(
        color = guide_legend(
          direction = "vertical", color = "white",
          override.aes = list(fill = NA)  # 设置置信区间的填充为透明
        )
      )+
      geom_segment(aes(x = 1, xend =17,y = 0.6, yend = (-0.1*(17-1)+0.6)),
                   color = "#A8021B", linewidth = 1.5,lty=1)
  }
  
  #@table S3
  {
    #rs25
    model <- lmer(LNRS25 ~ T5 + (1 | MEASURMENT), data = data_global) 
    summary(model)
    
    #rs.gt
    model <- lmer(LNRSGT ~ T5 + (1 | MEASURMENT), data = data_global) 
    summary(model)
  }
  
  #@table S4
  {
    anova_model <- aov(LNRS25 ~ MEASURMENT, data = data_global)
    summary(anova_model)
    anova_meas <- Anova(anova_model, type="II")
    print(anova_meas)
  }
  
  #@table S5
  {
    # site-mean value(rs25 and rs.gt)
    lm_site <- lm(lnrs25_sitemean~T5_sitemean + factor(MEASURMENT),data=data_site)
    summary(lm_site)
    lm_site <- lm(lnrsgt_sitemean~T5_sitemean + factor(MEASURMENT),data=data_site)
    summary(lm_site)
    model <- lmer(lnrs25_sitemean~T5_sitemean + (1 | MEASURMENT), data = data_site) 
    summary(model)
  }
  
  #@table S6
  {
    # without temperature standardization,with individual data
    names(reliability) <- 'T5'
    eiv.rs <- eivreg(lnrs.mtem~T5,data=data_tem,reliability = reliability)
    summary(eiv.rs)
    
    # without temperature standardization,with site-mean data
    eiv.rs <- eivreg(lnrs~T5,data=bin_mean,reliability = reliability)
    summary(eiv.rs)
  }
  
  #@table S9
  {
    eiv.rs <- eivreg(LNRS25~T5,data=pinus,reliability = reliability)
    summary(eiv.rs)
    confint(eiv.rs)
  }
}

# Figure 2
{
  fit_new <- lm(LNRS25 ~ tc6+lnE13, data = data_seasonal)
  summary(fit_new)
  vif(fit_new)
  
  par(mfrow = c(1, 2))
  visreg(fit_new)
  
  #Figure 2A: partial plot of temperature
  visreg(fit_new, "tc6", gg = TRUE,
         points = list(size = 3.6, pch = 16, alpha = 0),  
         fill = list(fill = "grey", alpha = 0.5),
         line = list(col = "black", size = 2.5, lty = 1)) +
    geom_point(color = "#002FA7",alpha=0.4, shape = 16, size = 3.6)+
    
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill = NA, color = "black", size = 1.3, linetype = "solid"),
          text = element_text(face = "bold"),
          axis.title.x = element_text(vjust = 0, size = 16),
          axis.title.y = element_text(vjust = 0, size = 16),
          axis.text.x = element_text(size = 12, color = "black"),
          axis.text.y = element_text(size = 12, color = "black")) +
    coord_cartesian(ylim = c(-4.2, 2), xlim = c(8, 19))+# 设置坐标轴范围
    # Style
    theme_bw()+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          panel.border = element_rect(fill=NA,color="black", size=1.3, linetype="solid"))+
    theme(text = element_text(family ="A"))+
    
    # Axis
    labs(x=expression(paste(italic(T)[g]," (°C)")),
         y=expression(paste(ln," ",italic(r[s25]),"  (nmol ",CO[2]," ",g^-1, s^-1,")")))+
    theme(axis.title.x = element_text(vjust = 0,size=17,family="Arial",face="bold"),
          axis.title.y = element_text(vjust = 0,size=17,family="Arial",face="bold"))+
    theme(axis.text.x = element_text(size = 12,color=rgb(0.2,0.2,0.2),family="Arial"),
          axis.text.y = element_text(size = 12,color=rgb(0.2,0.2,0.2),family="Arial"))+
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "black", size = 2.5,lty=1)+
    geom_segment(aes(x = 9, xend =18.5,y = 0.295, yend = (-0.1*(18.5-9)+0.295)),
                 color = "#A8021B",alpha=0.5, size = 1.5)+
    annotate("text", x = 8.5, y = 1.9, label = "A",size=6,fontface="bold")
  
  
  # Figure 2B: partial plot of transpiration
  visreg(fit_new, "lnE13", gg = TRUE,
         points = list(size = 3.6, pch = 16, alpha = 0),  # 设置透明度
         fill = list(fill = "gray", alpha = 0.5),
         line = list(col = "black", size = 1.1, lty = 1)) +
    geom_point(color = rgb(101/255,68/255,150/255), alpha = 0.5, shape = 16, size = 3.6)+ 
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill = NA, color = "black", size = 1.3, linetype = "solid"),
          text = element_text(family = "A", face = "bold"),
          axis.title.x = element_text(vjust = 0, size = 16),
          axis.title.y = element_text(vjust = 0, size = 16),
          axis.text.x = element_text(size = 12, color = "black"),
          axis.text.y = element_text(size = 12, color = "black")) +
    coord_cartesian(ylim = c(-4.2,2),xlim = c(3.9,6.2))+# 设置坐标轴范围
    # Style
    theme_bw()+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          panel.border = element_rect(fill=NA,color="black", size=1.3, linetype="solid"))+
    theme(text = element_text(family ="Arial"))+
    
    # Axis
    #labs(x=expression(paste("ln ", italic(E), " (mm"," ",day^-1,")")),
    labs(x=expression(paste("ln ", italic(E), " (nmol"," ",m^-2, s^-1,")")),
         y=expression(paste(ln," ",italic(r[s25]),"  (nmol ",CO[2]," ",g^-1, s^-1,")")))+
    theme(axis.title.x = element_text(vjust = 0,size=17,family="Arial",face="bold"),
          axis.title.y = element_text(vjust = 0,size=17,family="Arial",face="bold"))+
    theme(axis.text.x = element_text(size = 12,color=rgb(0.2,0.2,0.2),family="Arial"),
          axis.text.y = element_text(size = 12,color=rgb(0.2,0.2,0.2),family="Arial"))+
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "black", size = 2.5,lty=1)+
    geom_segment(aes(x = 4, xend =6.1,y = -1.9, yend = (6.1-4)-1.9),
                 color = "#C82015",alpha=0.5, size = 1.5)+
    annotate("text", x = 4, y = 1.9, label = "B",size=6,fontface="bold")
}

# fig. S5-6 and table S7
{
  #------------------------------------------------
  #VIF analysis and/or RMSE and R2
  #------------------------------------------------
  data_seasonal_timescale <- data_seasonal
  
  results_df <- data.frame(Model = character(),
                           tc = character(), 
                           lnE = character(),  
                           Slope_tc_e = numeric(), 
                           Slope_lne = numeric(),
                           R2_e = numeric(), 
                           RMSE_e = numeric(), 
                           VIF_e = numeric(),
                           AIC_e = numeric(),
                           stringsAsFactors = FALSE)
  
  # Combinations
  temperature_features <- paste0("tc", seq(1, 21))
  transpiration_features <- paste0("lnE",seq(1,21))
  combinations_e <- expand.grid(temperature_features, transpiration_features)
  
  # regression for lnrs25 ~ tc + lnE
  pb <- progress_bar$new(total = 441)
  for (i in 1:nrow(combinations_e)) {
    formula <- paste("LNRS25 ~", combinations_e[i, 1], "+", combinations_e[i, 2])
    fit <- lm(formula, data = data_seasonal_timescale)
    
    # fitting result
    summary_fit <- summary(fit)
    
    # check the significance, F-statistic p-value 
    pvalue1 <- summary_fit$coefficients[2, 4]
    pvalue2 <- summary_fit$coefficients[3, 4]
    
    if (pvalue1 > 0.05 || pvalue2>0.05) {  # p value 0.05  
      
      next
      
    } else{
      slope_tc <- coef(fit)[2]
      slope_lne <- coef(fit)[3]
      r_squared <- summary_fit$r.squared
      
      predicted_values <- predict(fit, data_seasonal_timescale)
      # RMSE
      rmse <- sqrt(mean((data_seasonal_timescale$LNRS25 - predicted_values)^2))
      
      # VIF
      vif <- vif(fit)[[1]]
      
      # AIC value
      aic_value <- AIC(fit)
      
      # print result
      results_df <- rbind(results_df, data.frame(Model = formula, 
                                                 tc = combinations_e[i,1],
                                                 lnE = combinations_e[i,2],
                                                 Slope_tc_e = slope_tc, 
                                                 Slope_lne = slope_lne, 
                                                 R2_e = r_squared, 
                                                 RMSE_e = rmse, VIF_e=vif,
                                                 AIC_e = aic_value,
                                                 stringsAsFactors = FALSE))
    }
    pb$tick()
    Sys.sleep(1/441)
  }
  
  #@fig.S5
  {
    #------------------------------------------------
    ####Different time scale of tc & lnE
    MR_e <- results_df
    r2_df_e <- subset(MR_e, select = c("tc", "lnE", "R2_e"))
    rmse_df_e <- subset(MR_e, select = c("tc", "lnE", "RMSE_e"))
    vif_df_e <- subset(MR_e, select = c("tc", "lnE", "VIF_e"))
    
    # FIGURE1:R-square
    r2_df_e$tc <- factor(r2_df_e$tc, levels = paste0("tc", 1:21))
    r2_df_e$lnE <- factor(r2_df_e$lnE, levels = paste0("lnE", 1:21))
    r2_df_e <- na.omit(r2_df_e)
    
    ggplot(r2_df_e, aes(x = tc, y = lnE)) +
      geom_point(aes(color = R2_e), size = 8, shape = 16) +
      scale_color_gradient2(low = "blue", mid = "white", high = "red",
                            midpoint = median(r2_df_e$R2_e),
                            limits = range(r2_df_e$R2_e),
                            #guide = guide_legend(orientation = "horizontal")
                            name="R-square") +
      # Style
      theme_bw()+
      theme(
        #panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(fill=NA,color="black", linewidth =0.5, linetype="solid"))+
      theme(text = element_text(family ="A"))+
      
      # Axis
      labs(x=expression(paste("Time window of ",italic(T[g]))),
           y=expression(paste("Time window of ln ",italic(E))))+
      theme(legend.position = "bottom",
            axis.text.x = element_blank(),  
            axis.text.y = element_blank(),
            axis.title.x = element_text(vjust = 0, size=17,face="bold"),
            axis.title.y = element_text(vjust = 0,size=17,face="bold"))
    
    # FIGURE2: RMSE
    rmse_df_e$tc <- factor(rmse_df_e$tc, levels = paste0("tc", 1:21))
    rmse_df_e$lnE <- factor(rmse_df_e$lnE, levels = paste0("lnE", 1:21))
    
    ggplot(rmse_df_e, aes(x = tc, y = lnE)) +
      geom_point(aes(color=RMSE_e),size=8,shape=16) +
      scale_color_gradient2(low = "red", mid = "white", high = "yellow", 
                            midpoint = median(rmse_df_e$RMSE_e),
                            limits = c(min(rmse_df_e$RMSE_e), 
                                       max(rmse_df_e$RMSE_e)),
                            name="RMSE")+
      # Style
      theme_bw()+
      theme(
        #panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(fill=NA,color="black", linewidth =0.5, linetype="solid"))+
      theme(text = element_text(family ="A"))+
      
      # Axis
      labs(x=expression(paste("Time window of ",italic(T[g]))),
           y=expression(paste("Time window of ln ",italic(E))))+
      theme(legend.position = "bottom",
            axis.text.x = element_blank(),  
            axis.text.y = element_blank(),
            axis.title.x = element_text(vjust = 0, size=17,face="bold"),
            axis.title.y = element_text(vjust = 0,size=17,face="bold"))
  }
  
  #@table S7
  {
    #table S7 is the top 10 of MR_e
  }
  
  #@fig.S6
  {
    data_seasonal$RS_MASS.GT <- data_seasonal$RS_MASS*2.2^((data_seasonal$tc6-data_seasonal$MTEM)/10)
    data_seasonal$LNRSGT <- log(data_seasonal$RS_MASS.GT)
    fit_new <- lm(LNRSGT ~ tc6+lnE13, data = data_seasonal)
    summary(fit_new)
    vif(fit_new)
    
    #fig.S6A: partial plot of temperature
    visreg(fit_new, "tc6", gg = TRUE,
           points = list(size = 3.6, pch = 16, alpha = 0),  
           fill = list(fill = "grey", alpha = 0.5),
           line = list(col = "black", size = 2.5, lty = 1)) +
      geom_point(color = "#002FA7",alpha=0.4, shape = 16, size = 3.6)+
      
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(fill = NA, color = "black", size = 1.3, linetype = "solid"),
            text = element_text(face = "bold"),
            axis.title.x = element_text(vjust = 0, size = 16),
            axis.title.y = element_text(vjust = 0, size = 16),
            axis.text.x = element_text(size = 12, color = "black"),
            axis.text.y = element_text(size = 12, color = "black")) +
      coord_cartesian(ylim = c(-5, 1), xlim = c(8, 19))+# 设置坐标轴范围
      # Style
      theme_bw()+
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
            panel.border = element_rect(fill=NA,color="black", size=1.3, linetype="solid"))+
      theme(text = element_text(family ="A"))+
      
      # Axis
      labs(x=expression(paste(italic(T)[g]," (°C)")),
           y=expression(paste(ln," ",italic(r[s.gt]),"  (nmol ",CO[2]," ",g^-1, s^-1,")")))+
      theme(axis.title.x = element_text(vjust = 0,size=17,family="Arial",face="bold"),
            axis.title.y = element_text(vjust = 0,size=17,family="Arial",face="bold"))+
      theme(axis.text.x = element_text(size = 12,color=rgb(0.2,0.2,0.2),family="Arial"),
            axis.text.y = element_text(size = 12,color=rgb(0.2,0.2,0.2),family="Arial"))+
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "black", size = 2.5,lty=1)+
      geom_segment(aes(x = 9, xend =18.5,y = -0.95, yend = (-0.023*(18.5-9)-0.95)),
                   color = "#A8021B",alpha=0.5, size = 1.5)+
      annotate("text", x = 8.5, y = 0.9, label = "A",size=6,fontface="bold")
    
    
    #fig.S6B: partial plot of transpiration
    visreg(fit_new, "lnE13", gg = TRUE,
           points = list(size = 3.6, pch = 16, alpha = 0),  # 设置透明度
           fill = list(fill = "gray", alpha = 0.5),
           line = list(col = "black", size = 1.1, lty = 1)) +
      geom_point(color = rgb(101/255,68/255,150/255), alpha = 0.5, shape = 16, size = 3.6)+ 
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(fill = NA, color = "black", size = 1.3, linetype = "solid"),
            text = element_text(family = "A", face = "bold"),
            axis.title.x = element_text(vjust = 0, size = 16),
            axis.title.y = element_text(vjust = 0, size = 16),
            axis.text.x = element_text(size = 12, color = "black"),
            axis.text.y = element_text(size = 12, color = "black")) +
      coord_cartesian(ylim = c(-5,1),xlim = c(3.9,6.2))+# 设置坐标轴范围
      # Style
      theme_bw()+
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
            panel.border = element_rect(fill=NA,color="black", size=1.3, linetype="solid"))+
      theme(text = element_text(family ="Arial"))+
      
      # Axis
      #labs(x=expression(paste("ln ", italic(E), " (mm"," ",day^-1,")")),
      labs(x=expression(paste("ln ", italic(E), " (nmol"," ",m^-2, s^-1,")")),
           y=expression(paste(ln," ",italic(r[s.gt]),"  (nmol ",CO[2]," ",g^-1, s^-1,")")))+
      theme(axis.title.x = element_text(vjust = 0,size=17,family="Arial",face="bold"),
            axis.title.y = element_text(vjust = 0,size=17,family="Arial",face="bold"))+
      theme(axis.text.x = element_text(size = 12,color=rgb(0.2,0.2,0.2),family="Arial"),
            axis.text.y = element_text(size = 12,color=rgb(0.2,0.2,0.2),family="Arial"))+
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "black", size = 2.5,lty=1)+
      geom_segment(aes(x = 4, xend =6.1,y = -2.8, yend = (6.1-4)-2.8),
                   color = "#C82015",alpha=0.5, size = 1.5)+
      annotate("text", x = 4, y = 0.9, label = "B",size=6,fontface="bold")
  }
}

# Figure 3
{
  species <- unique(w_data$FINAL_SPNAME)
  lm_species <- lm(LNRS25~GROWTH_T,data=w_data)
  summary(lm_species)
  
  #scatter plot
  scatter_plot <- ggplot(w_data,aes(x=GROWTH_T,y=LNRS25,color=FINAL_SPNAME))+
    # Style
    theme_bw()+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          panel.border = element_rect(fill=NA,color="black", size=1.3, linetype="solid"))+
    theme(text = element_text(family ="A",face = "bold"))+
    
    # Ponit
    geom_point(shape="circle",size=4)+
    scale_color_manual(values = c(rgb(46/256,85/256,150/256),
                                  rgb(159/256,189/256,215/256),
                                  "black","#30C9E8","grey"))+
    
    # Line
    geom_smooth(method = "lm",formula = y~x,se=T,linetype=2,size=1,
                color="black")+
    scale_fill_manual(values=c(rgb(243/256,141/256,147/256,0.1)))+
    
    # Axis
    labs(x=expression(paste(italic(T[g])," (°C)")),
         y=expression(paste(ln," ",italic(r[s25]),"  (nmol ",CO[2]," ",g^-1, s^-1,")")))+
    theme(axis.title.x = element_text(vjust = 0,size=16),
          axis.title.y = element_text(vjust = 0,size=16))+
    theme(axis.text.x = element_text(size = 12,color="black"),
          axis.text.y = element_text(size = 12,color="black"))+
    
    ylim(1,6)+
    xlim(15,35)+
    geom_segment(aes(x = 15, xend =35,y = 3.7, yend = (-0.1*(35-15)+3.7)),
                 color = "#A8021B",alpha=0.5, size = 1.5)+
    
    theme( legend.position =" none ")+
    annotate("text", x = 16, y = 6, label = "A",size=6,fontface="bold")
  
  Ba <- data.table(w_data[which(w_data$FINAL_SPNAME==species[1]),])
  Pnigra <- data.table(w_data[which(w_data$FINAL_SPNAME==species[2]),])
  Ppinaster <- data.table(w_data[which(w_data$FINAL_SPNAME==species[3]),])
  Ppinea <- data.table(w_data[which(w_data$FINAL_SPNAME==species[4]),])
  Ps <- data.table(w_data[which(w_data$FINAL_SPNAME==species[5]),])
  
  color_list <- c(rgb(46/256,85/256,150/256),
                  rgb(159/256,189/256,215/256),
                  "black","#30C9E8","grey")
  
  species_data <- list(Ba,Pnigra,Ppinaster,Ppinea,Ps)
  plots <- list()
  label_list <- c("B","C","D","E","F")
  for(i in 1:5){
    species_individual <- data.table(species_data[[i]])
    summary_stats <- species_individual %>%  
      group_by(GROWTH_T) %>%  
      summarise(  
        min_val = min(LNRS25),  
        max_val = max(LNRS25),  
        mean_val = mean(LNRS25)  
      )  
    p <- ggplot(species_individual, aes(x = GROWTH_T, y = LNRS25)) +  
      # 绘制表示分布范围的直线（从最小值到最大值）  
      geom_segment(data = summary_stats, aes(xend = GROWTH_T, y = min_val, yend = max_val),   
                   color = color_list[i], size = 2.5) +  # 调整size参数来加粗线条
      # 标注最小值点  
      geom_point(data = summary_stats, aes(y = min_val), color = color_list[i], shape = "-", size = 18) +  
      # 标注最大值点  
      geom_point(data = summary_stats, aes(y = max_val), color = color_list[i], shape = "-", size = 18) +  
      # 标注平均值点  
      geom_point(data = summary_stats, aes(y = mean_val), color = "#C82015", shape = "-", size = 18) +  
      theme_minimal() +  
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.border = element_rect(fill = NA, color = "black", size = 1.3, linetype = "solid")) +
      ylim(1,6)+
      xlim(12,37)+
      theme(  
        axis.title.x = element_blank(),  # 移除x轴标题  
        axis.title.y = element_blank(),  # 移除y轴标题  
        axis.text.x = element_blank(),   # 移除x轴刻度标签  
        axis.text.y = element_text(size = 14, color = "black")
      )  +
      geom_smooth(method = "lm",formula = y~x,linetype=2,size=1,
                  color="black",fill=rgb(0.8,0.8,0.8,0.1))+
      annotate("text", x = 35, y = 5.5, label = label_list[i],size=6,fontface="bold")
    
    plots[[i]] <- p
  }
  
  small_plots_grid <- plot_grid(  
    plots[[1]],  plots[[2]], plots[[3]], plots[[4]], plots[[5]],
    ncol = 1, nrow = 5,                  
    rel_widths = c(1),                   
    rel_heights = c(1, 1, 1, 1, 1)       
  )  
  final_grid <- plot_grid(  
    scatter_plot, small_plots_grid,  
    ncol = 2, nrow = 1,                 
    rel_widths = c(3, 1)                
  )  
  
  final_grid
}

# fig.S8 and table S8
{
  #@fig.S8
  {
    ggplot(w_data,aes(x=GROWTH_T,y=LNRSGT,color=FINAL_SPNAME))+
      # Style
      theme_bw()+
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
            panel.border = element_rect(fill=NA,color="black", linewidth =1.3, linetype="solid"))+
      theme(text = element_text(family ="HEL"))+
      
      # Ponit
      geom_point(shape="circle",size=4)+
      scale_color_manual(values = c(rgb(46/256,85/256,150/256),
                                    rgb(159/256,189/256,215/256),
                                    "black","#30C9E8","grey"))+
      
      # Line
      geom_smooth(method = "lm",formula = y~x,se=T,linetype=2,size=1,
                  color="black")+
      scale_fill_manual(values=c(rgb(243/256,141/256,147/256,0.1)))+
      
      # Axis
      labs(x=expression(paste(italic(T[g])," (°C)")),
           y=expression(paste("ln ",italic(r[sgt]),"  (nmol ",CO[2]," ",g^-1, s^-1,")")))+
      theme(axis.title.x = element_text(vjust = 0,size=16),
            axis.title.y = element_text(vjust = 0,size=16))+
      theme(axis.text.x = element_text(size = 12,color="black"),
            axis.text.y = element_text(size = 12,color="black"))+
      geom_segment(aes(x = 15, xend =35,y = 2.8, yend = (-0.0226*(35-15)+2.8)),
                   color = "#A8021B", size = 2,lty=1)+
      theme(legend.position = 'none')
    
    
    lm_species <- lm(LNRSGT~GROWTH_T,data=w_data)
    summary(lm_species)
  }
  
  #@table S8
  {
    lm_species <- lm(LNRS25~GROWTH_T,data=w_data)
    summary(lm_species)
    species_data <- list(Ba,Pnigra,Ppinaster,Ppinea,Ps)
    for(i in 1:5){
      species_individual <- data.table(species_data[[i]])
      lm_species <- lm(LNRS25~GROWTH_T,data=species_individual)
      f <- summary(lm_species)
      lower <- f$coefficients[2,1]-f$coefficients[2,2]
      upper <- f$coefficients[2,1]+f$coefficients[2,2]
      ince <- f$coefficients
      print(lower)
      print(upper)
      print(ince)
      
      
      print("========")
      print("========")
    }
  }
}

# Figure 4 and table S10
{
  #Figure 4A
  {
    ##global simulation
    {
      #======stem biomass======
      #Unit:above-ground biomass, Mg ha−1
      #spatial resolution: 0.5 degree
      #time:2010+-1year
      data_bio <- raster(".../global.tif")
      
      #======sapwood biomass===
      #Unit:above-ground biomass, Mg ha−1
      n <- 2.5
      a <- 0.365
      b <- 0.8771
      c <- 0.00083
      d <- 0.7478
      
      msap <- data_bio^((b-1)*d+1)* n/((n-1)*b+1)*a*c^(b-1)
      
      #======mgdd5=======
      #Unit:degree
      mgdd5_csv <- read.csv(".../mgdd5.csv",sep=";")#lat lon
      cols <- c(3,2,6)
      mgdd5_csv <- mgdd5_csv[,cols]#lon lat
      mgdd5_raster <- rasterFromXYZ(mgdd5_csv,digits = 0.2)
      data_mgdd5 <- resample(mgdd5_raster,data_bio)
      fT <- 2.2^((data_mgdd5-25)/10)
      
      #======Tgrowingtime======
      #Unit:(time) s
      #GD:growing days
      cols <- c(3,2,4)
      time <- mgdd5_csv[,cols]
      time_raster <- rasterFromXYZ(time,digits = 0.2)
      data_time <- resample(time_raster,data_bio)
      time <- data_time*24*3600
      
      #======area======
      #Unit:km2
      lat <- raster(".../lat")
      data_lat <- resample(lat,data_bio)
      pixel <- raster(resolution=0.5)
      tem_lat <- cos((pi*data_lat/180))
      sita <- abs(tem_lat)
      area <- 110*0.5*110*0.5*sita
      
      #======rs25 ======
      #Unit:nmol CO2 g-1 s-1
      #convert to per sapwood mass based
      
      rs_mtem25 <- data.table(data_global[which((data_global$MTEM>=24)&(data_global$MTEM<=26)),])
      rs_mtem25$lnrssap.mtem <- log(rs_mtem25$RS_MASS)
      
      #plot of lnrs.mtem ~ T5
      ggplot(rs_mtem25, aes(x = T5, y = lnrssap.mtem)) +
        geom_point(size = 2, alpha=0.2) +  
        # Line
        stat_smooth(method = "lm",formula = y~x,se=T,lty=2,linewidth=1.2)+
        labs(title = "lnrs25 without temperature standardization",
             x = "T5",
             y = "lnrs.mtem")
      
      ##@ table S10fix mgdd5 slope of -0.1 
      lm <- lm (lnrssap.mtem ~ offset (-0.1*T5), data=rs_mtem25) 
      summary(lm)
      intercept <- summary(lm)$coefficient[1]
      
      rs25 <- exp(intercept-2.5)
      rs25_a <- exp(log(rs25)-0.1*(data_mgdd5-25))
      
      #======Rs and total emission======
      Rs_a <- rs25_a*fT*msap*time*1.2*10^(-6)  #g C/m2 year
      plot(Rs_a)
      Rs_sum_a <- rs25_a*fT*msap*time*area*1.2*10^(-15)  #Pg C/year
      Rs_year_a <- cellStats(Rs_sum_a,stat='sum')
      Rs_year_a
      plot(Rs_sum_a)
    }
    
    #uncertainty
    {
      ## calculate T5 uncertainty
      {
        # data prepare
        nc <- nc_open("cru_ts4.04.1901.2019.tmp.dat.nc")#lon lat time
        time <- ncvar_get(nc=nc,varid="time")
        lon <- ncvar_get(nc=nc,varid = "lon")
        lat <- ncvar_get(nc=nc,varid = "lat")
        tmp <- ncvar_get(nc=nc,varid = "tmp")
        
        #*1901-2020
        CRUts <- ncvar_get(nc=nc,varid = 'tmp',start=c(1,1,721),count = c(-1,-1,360))
        TSyear <- array(dim=c(720,360,30,13)) 
        tssum <- 0
        n <- 0
        pb <- progress_bar$new(total = 720, clear = FALSE)
        for(lon in 1:720){#-179.75 - 179.75
          
          for(lat in 1:360){#-89.75 - 89.75
            
            for(year in 1:30){#30
              tssum <- 0
              n <- 0
              for(mon in 1:12){
                tstmp <- CRUts[lon,lat,((year-1)*12+mon)]
                TSyear[lon,lat,year,mon] <- tstmp
                if(is.na(tstmp)){
                  TSyear[lon,lat,year,mon] <- NA
                }else if(tstmp>=5){
                  tssum <- tssum+tstmp
                  n <- n+1
                }
              }
              TSyear[lon,lat,year,13] <- tssum/n
            }
          }
          pb$tick()
          Sys.sleep(1 /720)
        }
        
        #se 1
        se1 <- array(dim = c(720,360))
        mgdd5_oneyear <- 0
        pb <- progress_bar$new(total = 720, clear = FALSE)
        for(lon in 1:720){
          for(lat in 1:360){
            mgdd5_onepixel <- TSyear[lon,lat,,13]
            sdTS <- sd(mgdd5_onepixel)
            nTS <- 30 #n
            se1[lon,lat] <- sdTS/sqrt(nTS)
          }
          pb$tick()
          Sys.sleep(1 /720)
        }
        
        #se 2
        clv2 <- read.csv('.../CRUCLv2_mgdd5.csv', header = TRUE, sep = ",", quote = "\"",
                         dec = ".", fill = TRUE, comment.char = "")
        se2 <- array(dim = c(720,360))
        pb <- progress_bar$new(total = 720, clear = FALSE)
        for(lon in 1:720){
          for(lat in 1:360){
            lonreal <- -179.75+(lon-1)*0.5
            latreal <- -89.75+(lat-1)*0.5
            lonmax <- lonreal+0.25
            lonmin <- lonreal-0.25
            latmax <- latreal+0.25
            latmin <- latreal-0.25
            
            CLpixel <- data.table(clv2[which((clv2$V1>latmin)&(clv2$V1<latmax)&(clv2$V2<lonmax)
                                             &(clv2$V2>lonmin)),])
            sdCL <- sd(CLpixel$V15)
            nCL <- 9
            se2[lon,lat] <- sdCL/sqrt(nCL)
          }
          pb$tick()
          Sys.sleep(1 /720)
        }
        
        #SE_T5
        se_t5 <- array(dim = c(720,360))
        pb <- progress_bar$new(total = 720, clear = FALSE)
        for(lon in 1:720){
          for(lat in 1:360){
            se_t5[lon,lat] <- sqrt((se1[lon,lat]^2+se2[lon,lat]^2)/2)
          }
          pb$tick()
          Sys.sleep(1 /720)
        }
        global <- mean(se_t5,na.rm = T)#global reliability = 0.9628
        print(global)
        #global se of T5 = 0.14
      }

      ##Variables Uncertainty
      se_rs25 <- 0.03485
      se_q10 <- 0.2
      se_mgdd5 <- 0.14
      se_a <- 0.051
      se_b <- 0.029
      se_c <- 0.000125
      se_d <- 0.0189
      
      santoroAGB_SE <- raster(".../global_err.tif")
      
      #SE
      SE <- sqrt((se_rs25)^2+(se_q10)^2+(se_mgdd5)^2+(se_a)^2+(se_b)^2+(se_c)^2+(se_d)^2+(santoroAGB_SE)^2)
      SE_NA <- calc(SE, fun = function(x) { x[x == 0] <- NA; return(x) })
      SE_mean <- cellStats(SE_NA,stat='mean')
      SE_mean
    }
  }
 
  #Figure 4B
  {
    ## future prediction
    {
      setwd(".../CMIP6")
      d_ws<-".../CMIP6/"
      #文件列表
      lsfCABLE<-tibble(nmf=dir(d_ws,pattern = "tas_day_ACCESS-ESM1-5"))
      lsfORCHIDEE<-tibble(nmf=dir(d_ws,pattern = "tas_day_IPSL-CM6A-LR"))
      lsfJULES<-tibble(nmf=dir(d_ws,pattern = "tas_day_UKESM1-0-LL"))
      lsfCLM5<-tibble(nmf=dir(d_ws,pattern = "tas_day_CESM2"))
      
      #CONSTANT (without considering thermal acclimation)
      {
        Rs <- rs25*fT*msap*time*1.2*10^(-6)  #g C/m2 year
        Rs_sum <- rs25*fT*msap*time*area*1.2*10^(-15)  #Pg C/year
        Rs_year <- cellStats(Rs_sum,stat='sum')
        Rs_year
        
        Rs_a <- rs25_a*fT*msap*time*1.2*10^(-6)  #g C/m2 year
        Rs_sum_a <- rs25_a*fT*msap*time*area*1.2*10^(-15)  #Pg C/year
        Rs_year_a <- cellStats(Rs_sum_a,stat='sum')
        Rs_year_a
        
        rs25_a_mean <- cellStats(rs25_a,stat='mean')
      }
      
      # calculate number of days in one year  
      # for cabel,orchidee, consider leap year
      # for jules, 360 per year
      # for clm, 365 per year
      days_in_year <- function(year) {  
        # leap year or not
        if ((year%%4==0 & year%%100 != 0) | ((year%%400) == 0)) {  
          return(366)  # leap year 366 
        } else {  
          return(365)  # common year 365
        }  
      }  
      
      #================================================
      #cable
      #================================================
      ###----SSP126------------------------------------
      data_cable126 <- array(dim=c(86,6))
      data_cable126[,1] <- 2015:2100
      #2015-2064
      {
        nc <- nc_open(lsfCABLE[[1,1]])
        lon <- ncvar_get(nc, "lon")
        lon_len <- length(lon)
        lat <- ncvar_get(nc, "lat")
        lat_len <- length(lat)
        
        point_mat <- expand.grid(lon = lon-180, lat = lat)
        coordinates(point_mat) <- ~lon+lat
        proj4string(point_mat) <- proj4string(data_bio)

        F_bio <- raster::extract(data_bio, point_mat) 
        F_fT <- raster::extract(fT, point_mat)
        F_area <- raster::extract(area, point_mat)
        F_msap <- raster::extract(msap,point_mat)
        F_rs25a <- raster::extract(rs25_a,point_mat)

        result <- data.frame(point_mat@coords, F_bio = F_bio, F_msap = F_msap, F_area = F_area, F_rs25a=F_rs25a)
        
        pb <- progress_bar$new(total = 50, clear = FALSE)
        days_sum <- 0
        for(y in 2015:2064){
          day_num_year <- days_in_year(y)
          day_start <- days_sum+1
          days_sum <- days_sum+day_num_year

          temp_val <- ncvar_get(nc, "tas",start = c(1,1,day_start),count=c(-1,-1,day_num_year))

          temp_val[temp_val < 278.15] <- NA

          output <- lapply(seq_along(lat), function(j) {
            lapply(seq_along(lon), function(i) {
              
              num <- (j-1)*lon_len+i
              
              tas_values <- temp_val[i, j, ]
              new_time <- (365-sum(is.na(tas_values)))*24*3600 #new_time
              mean_above_5 <- mean(tas_values)-278.15 #new_mgdd5
              new_fT <- 2.2^((mean_above_5-25)/10) #new_fT
              
              F_bio <- result[num,3]
              F_msap <- result[num,4]
              F_area <- result[num,5]
              F_rs25a <- result[num,6]
              
              #not consider thermal acclimation
              new_Rs <- F_rs25a*new_fT*F_msap*new_time*F_area*1.2*10^(-15)  #Pg C/year
              #new_Rs_year <- cellStats(new_Rs,stat='sum')
              
              #consider thermal acclimation
              new_rs25a <- exp(log(rs25)-0.1*(mean_above_5-25))
              new_Rsa <- new_rs25a*new_fT*F_msap*new_time*F_area*1.2*10^(-15)  #Pg C/year
              #new_Rsa_year <- cellStats(new_Rsa,stat='sum')
              
              # recudtion of carbon emission
              reduction <- new_Rsa-new_Rs
              
              #fraction
              fraction <- -reduction/new_Rs
              return(list(lon = lon[i]-180,lat = lat[j], mean_above_5 = mean_above_5,time=new_time,new_Rs_year=new_Rs,new_Rsa_year=new_Rsa,
                          reduction=reduction,F_bio=F_bio,F_msap=F_msap,new_fT=new_fT,F_area=F_area,new_time=new_time,fraction=fraction))
            })
          })
          output <- data.frame(do.call(rbind, unlist(output, recursive = FALSE)))
          colnames(output) <- c("lon","lat","mgdd5","time","Rs","Rsa","Reduction","bio","msap","ft","area","time","fraction")
          output_raster <- rasterFromXYZ(output)
          mgdd5_raster <- resample(output_raster$mgdd5,data_bio)
          reduction_raster <- resample(output_raster$Reduction,data_bio)
          #calculate fraction
          Rs_raster <-  resample(output_raster$Rs,data_bio)
          Rsa_raster <-  resample(output_raster$Rsa,data_bio)
          fraction_raster <-  resample(output_raster$fraction,data_bio)
          
          data_cable126[(y-2014),2] <- cellStats(mgdd5_raster,stat='mean')
          data_cable126[(y-2014),3] <- cellStats(reduction_raster,stat='sum')
          data_cable126[(y-2014),4] <- cellStats(Rs_raster,stat='sum')
          data_cable126[(y-2014),5] <- cellStats(Rsa_raster,stat='sum')
          data_cable126[(y-2014),6] <- cellStats(fraction_raster,stat='mean')
          pb$tick()
          Sys.sleep(1 / 50)
        }
      }
      #2065-2100
      {
        nc <- nc_open(lsfCABLE[[2,1]])

        lon <- ncvar_get(nc, "lon")
        lon_len <- length(lon)
        lat <- ncvar_get(nc, "lat")
        lat_len <- length(lat)
        
        point_mat <- expand.grid(lon = lon-180, lat = lat)
        coordinates(point_mat) <- ~lon+lat
        proj4string(point_mat) <- proj4string(data_bio)
        
        F_bio <- raster::extract(data_bio, point_mat) 
        F_fT <- raster::extract(fT, point_mat)
        F_area <- raster::extract(area, point_mat)
        F_msap <- raster::extract(msap,point_mat)
        F_rs25a <- raster::extract(rs25_a,point_mat)

        result <- data.frame(point_mat@coords, F_bio = F_bio, F_msap = F_msap, F_area = F_area, F_rs25a=F_rs25a)
        
        pb <- progress_bar$new(total = 36, clear = FALSE)
        days_sum <- 0
        for(y in 2065:2100){

          day_num_year <- days_in_year(y)
          day_start <- days_sum+1
          days_sum <- days_sum+day_num_year

          temp_val <- ncvar_get(nc, "tas",start = c(1,1,day_start),count=c(-1,-1,day_num_year))
          
          temp_val[temp_val < 278.15] <- NA
 
          output <- lapply(seq_along(lat), function(j) {
            lapply(seq_along(lon), function(i) {
              
              num <- (j-1)*lon_len+i
              
              tas_values <- temp_val[i, j, ]
              new_time <- (365-sum(is.na(tas_values)))*24*3600 #new_time
              mean_above_5 <- mean(tas_values)-278.15 #new_mgdd5
              new_fT <- 2.2^((mean_above_5-25)/10) #new_fT
              
              F_bio <- result[num,3]
              F_msap <- result[num,4]
              F_area <- result[num,5]
              F_rs25a <- result[num,6]
              
              #not consider thermal acclimation
              new_Rs <- F_rs25a*new_fT*F_msap*new_time*F_area*1.2*10^(-15)  #Pg C/year
              #new_Rs_year <- cellStats(new_Rs,stat='sum')
              
              #consider thermal acclimation
              new_rs25a <- exp(log(rs25)-0.1*(mean_above_5-25))
              new_Rsa <- new_rs25a*new_fT*F_msap*new_time*F_area*1.2*10^(-15)  #Pg C/year
              #new_Rsa_year <- cellStats(new_Rsa,stat='sum')
              
              # recudtion of carbon emission
              reduction <- new_Rsa-new_Rs
              
              #fraction
              fraction <- -reduction/new_Rs
              
              return(list(lon = lon[i]-180,lat = lat[j], mean_above_5 = mean_above_5,time=new_time,new_Rs_year=new_Rs,new_Rsa_year=new_Rsa,
                          reduction=reduction,F_bio=F_bio,F_msap=F_msap,new_fT=new_fT,F_area=F_area,new_time=new_time,fraction=fraction))
            })
          })
          output <- data.frame(do.call(rbind, unlist(output, recursive = FALSE)))
          colnames(output) <- c("lon","lat","mgdd5","time","Rs","Rsa","Reduction","bio","msap","ft","area","time","fraction")
          output_raster <- rasterFromXYZ(output)
          mgdd5_raster <- resample(output_raster$mgdd5,data_bio)
          reduction_raster <- resample(output_raster$Reduction,data_bio)
          #calculate fraction
          Rs_raster <-  resample(output_raster$Rs,data_bio)
          Rsa_raster <-  resample(output_raster$Rsa,data_bio)
          fraction_raster <-  resample(output_raster$fraction,data_bio)
          
          data_cable126[(y-2014),2] <- cellStats(mgdd5_raster,stat='mean')
          data_cable126[(y-2014),3] <- cellStats(reduction_raster,stat='sum')
          data_cable126[(y-2014),4] <- cellStats(Rs_raster,stat='sum')
          data_cable126[(y-2014),5] <- cellStats(Rsa_raster,stat='sum')
          data_cable126[(y-2014),6] <- cellStats(fraction_raster,stat='mean')
          pb$tick()
          Sys.sleep(1 / 36)
        }
      }
      colnames(data_cable126) <- c("year","T5","reduction","Rs","Rsa","fraction")
      plot(data_cable126[,1],data_cable126[,2])#plot mgdd5
      plot(data_cable126[,1],data_cable126[,3])#plot carbon reduction
      plot(data_cable126[,1],data_cable126[,4])#plot Rs
      plot(data_cable126[,1],data_cable126[,5])#plot Rsa
      plot(data_cable126[,1],data_cable126[,6])#plot fraction
      plot(data_cable126[,1],-data_cable126[,3]/data_cable126[,4])#plot fraction(sum calculate)
      ###SSP585----------------------------------------
      data_cable585 <- array(dim=c(86,6))
      data_cable585[,1] <- 2015:2100
      #2015-2064
      {
        nc <- nc_open(lsfCABLE[[3,1]])
        lon <- ncvar_get(nc, "lon")
        lon_len <- length(lon)
        lat <- ncvar_get(nc, "lat")
        lat_len <- length(lat)
        
        point_mat <- expand.grid(lon = lon-180, lat = lat)
        coordinates(point_mat) <- ~lon+lat
        proj4string(point_mat) <- proj4string(data_bio)

        F_bio <- raster::extract(data_bio, point_mat) 
        F_fT <- raster::extract(fT, point_mat)
        F_area <- raster::extract(area, point_mat)
        F_msap <- raster::extract(msap,point_mat)
        F_rs25a <- raster::extract(rs25_a,point_mat)

        result <- data.frame(point_mat@coords, F_bio = F_bio, F_msap = F_msap, F_area = F_area, F_rs25a=F_rs25a)
        
        pb <- progress_bar$new(total = 50, clear = FALSE)
        days_sum <- 0
        for(y in 2015:2064){
          day_num_year <- days_in_year(y)
          day_start <- days_sum+1
          days_sum <- days_sum+day_num_year

          temp_val <- ncvar_get(nc, "tas",start = c(1,1,day_start),count=c(-1,-1,day_num_year))

          temp_val[temp_val < 278.15] <- NA

          output <- lapply(seq_along(lat), function(j) {
            lapply(seq_along(lon), function(i) {
              
              num <- (j-1)*lon_len+i
              
              tas_values <- temp_val[i, j, ]
              new_time <- (365-sum(is.na(tas_values)))*24*3600 #new_time
              mean_above_5 <- mean(tas_values)-278.15 #new_mgdd5
              new_fT <- 2.2^((mean_above_5-25)/10) #new_fT
              
              F_bio <- result[num,3]
              F_msap <- result[num,4]
              F_area <- result[num,5]
              F_rs25a <- result[num,6]
              
              #not consider thermal acclimation
              new_Rs <- F_rs25a*new_fT*F_msap*new_time*F_area*1.2*10^(-15)  #Pg C/year
              #new_Rs_year <- cellStats(new_Rs,stat='sum')
              
              #consider thermal acclimation
              new_rs25a <- exp(log(rs25)-0.1*(mean_above_5-25))
              new_Rsa <- new_rs25a*new_fT*F_msap*new_time*F_area*1.2*10^(-15)  #Pg C/year
              #new_Rsa_year <- cellStats(new_Rsa,stat='sum')
              
              # recudtion of carbon emission
              reduction <- new_Rsa-new_Rs
              
              #fraction
              fraction <- -reduction/new_Rs
              return(list(lon = lon[i]-180,lat = lat[j], mean_above_5 = mean_above_5,time=new_time,new_Rs_year=new_Rs,new_Rsa_year=new_Rsa,
                          reduction=reduction,F_bio=F_bio,F_msap=F_msap,new_fT=new_fT,F_area=F_area,new_time=new_time,fraction=fraction))
            })
          })
          output <- data.frame(do.call(rbind, unlist(output, recursive = FALSE)))
          colnames(output) <- c("lon","lat","mgdd5","time","Rs","Rsa","Reduction","bio","msap","ft","area","time","fraction")
          output_raster <- rasterFromXYZ(output)
          mgdd5_raster <- resample(output_raster$mgdd5,data_bio)
          reduction_raster <- resample(output_raster$Reduction,data_bio)
          #calculate fraction
          Rs_raster <-  resample(output_raster$Rs,data_bio)
          Rsa_raster <-  resample(output_raster$Rsa,data_bio)
          fraction_raster <-  resample(output_raster$fraction,data_bio)
          
          data_cable585[(y-2014),2] <- cellStats(mgdd5_raster,stat='mean')
          data_cable585[(y-2014),3] <- cellStats(reduction_raster,stat='sum')
          data_cable585[(y-2014),4] <- cellStats(Rs_raster,stat='sum')
          data_cable585[(y-2014),5] <- cellStats(Rsa_raster,stat='sum')
          data_cable585[(y-2014),6] <- cellStats(fraction_raster,stat='mean')
          pb$tick()
          Sys.sleep(1 / 50)
        }
      }
      #2065-2100
      {
        nc <- nc_open(lsfCABLE[[4,1]])

        lon <- ncvar_get(nc, "lon")
        lon_len <- length(lon)
        lat <- ncvar_get(nc, "lat")
        lat_len <- length(lat)

        point_mat <- expand.grid(lon = lon-180, lat = lat)
        coordinates(point_mat) <- ~lon+lat
        proj4string(point_mat) <- proj4string(data_bio)

        F_bio <- raster::extract(data_bio, point_mat) 
        F_fT <- raster::extract(fT, point_mat)
        F_area <- raster::extract(area, point_mat)
        F_msap <- raster::extract(msap,point_mat)
        F_rs25a <- raster::extract(rs25_a,point_mat)

        result <- data.frame(point_mat@coords, F_bio = F_bio, F_msap = F_msap, F_area = F_area, F_rs25a=F_rs25a)
        
        pb <- progress_bar$new(total = 36, clear = FALSE)
        days_sum <- 0
        for(y in 2065:2100){

          day_num_year <- days_in_year(y)
          day_start <- days_sum+1
          days_sum <- days_sum+day_num_year
          

          temp_val <- ncvar_get(nc, "tas",start = c(1,1,day_start),count=c(-1,-1,day_num_year))

          temp_val[temp_val < 278.15] <- NA
          
          output <- lapply(seq_along(lat), function(j) {
            lapply(seq_along(lon), function(i) {
              
              num <- (j-1)*lon_len+i
              
              tas_values <- temp_val[i, j, ]
              new_time <- (365-sum(is.na(tas_values)))*24*3600 #new_time
              mean_above_5 <- mean(tas_values)-278.15 #new_mgdd5
              new_fT <- 2.2^((mean_above_5-25)/10) #new_fT
              
              F_bio <- result[num,3]
              F_msap <- result[num,4]
              F_area <- result[num,5]
              F_rs25a <- result[num,6]
              
              #not consider thermal acclimation
              new_Rs <- F_rs25a*new_fT*F_msap*new_time*F_area*1.2*10^(-15)  #Pg C/year
              #new_Rs_year <- cellStats(new_Rs,stat='sum')
              
              #consider thermal acclimation
              new_rs25a <- exp(log(rs25)-0.1*(mean_above_5-25))
              new_Rsa <- new_rs25a*new_fT*F_msap*new_time*F_area*1.2*10^(-15)  #Pg C/year
              #new_Rsa_year <- cellStats(new_Rsa,stat='sum')
              
              # recudtion of carbon emission
              reduction <- new_Rsa-new_Rs
              
              #fraction
              fraction <- -reduction/new_Rs
              
              return(list(lon = lon[i]-180,lat = lat[j], mean_above_5 = mean_above_5,time=new_time,new_Rs_year=new_Rs,new_Rsa_year=new_Rsa,
                          reduction=reduction,F_bio=F_bio,F_msap=F_msap,new_fT=new_fT,F_area=F_area,new_time=new_time,fraction=fraction))
            })
          })
          output <- data.frame(do.call(rbind, unlist(output, recursive = FALSE)))
          colnames(output) <- c("lon","lat","mgdd5","time","Rs","Rsa","Reduction","bio","msap","ft","area","time","fraction")
          output_raster <- rasterFromXYZ(output)
          mgdd5_raster <- resample(output_raster$mgdd5,data_bio)
          reduction_raster <- resample(output_raster$Reduction,data_bio)
          #calculate fraction
          Rs_raster <-  resample(output_raster$Rs,data_bio)
          Rsa_raster <-  resample(output_raster$Rsa,data_bio)
          fraction_raster <-  resample(output_raster$fraction,data_bio)
          
          data_cable585[(y-2014),2] <- cellStats(mgdd5_raster,stat='mean')
          data_cable585[(y-2014),3] <- cellStats(reduction_raster,stat='sum')
          data_cable585[(y-2014),4] <- cellStats(Rs_raster,stat='sum')
          data_cable585[(y-2014),5] <- cellStats(Rsa_raster,stat='sum')
          data_cable585[(y-2014),6] <- cellStats(fraction_raster,stat='mean')
          pb$tick()
          Sys.sleep(1 / 36)
        }
      }
      
      plot(data_cable585[,1],data_cable585[,2])#plot mgdd5
      plot(data_cable585[,1],data_cable585[,3])#plot carbon reduction
      plot(data_cable585[,1],data_cable585[,4])#plot Rs
      plot(data_cable585[,1],data_cable585[,5])#plot Rsa
      plot(data_cable585[,1],data_cable585[,6])#plot fraction
      plot(data_cable585[,1],-data_cable585[,3]/data_cable585[,4])#plot fraction(sum calculate)
      
      #Plot single model
      data_cable_tem <- cbind(data_cable126[,1:2],data_cable585[,2])
      data_cable_red <- cbind(data_cable126[,c(1,3)],data_cable585[,3])
      data_cable_Rs <- cbind(data_cable126[,c(1,4)],data_cable585[,4])
      data_cable_Rsa <- cbind(data_cable126[,c(1,5)],data_cable585[,5])
      data_cable_fraction <- cbind(data_cable126[,c(1,6)],data_cable585[,6])
      
      
      colnames(data_cable_tem) <- c("Year","SSP126","SSP585")
      colnames(data_cable_red) <- c("Year","SSP126","SSP585")
      colnames(data_cable_Rs) <- c("Year","SSP126","SSP585")
      colnames(data_cable_Rsa) <- c("Year","SSP126","SSP585")
      colnames(data_cable_fraction) <- c("Year","SSP126","SSP585")
      
      data_cable_tem <- data.frame(data_cable_tem)
      data_cable_red <- data.frame(data_cable_red)
      data_cable_Rs <- data.frame(data_cable_Rs)
      data_cable_Rsa <- data.frame(data_cable_Rsa)
      data_cable_fraction <- data.frame(data_cable_fraction)
      
      data_cable_m1 <- reshape2::melt(id="Year",data_cable_tem)
      data_cable_m2 <- reshape2::melt(id="Year",data_cable_red)
      data_cable_m3 <- reshape2::melt(id="Year",data_cable_Rs)
      data_cable_m4 <- reshape2::melt(id="Year",data_cable_Rsa)
      data_cable_m5 <- reshape2::melt(id="Year",data_cable_fraction)
      data_cable <- cbind(data_cable_m1[,],data_cable_m2[,3],data_cable_m3[,3],data_cable_m4[,3],data_cable_m5[,3])
      colnames(data_cable) <- c("Year","Scenario","T5","reduction","Rs","Rsa","fraction")
      
      #plot
      {
        #Future reduction
        ggplot()+
          geom_point(data = data_cable,aes(x=Year,y=reduction,color=Scenario))+
          geom_line(data = data_cable,aes(x=Year,y=reduction,color=Scenario))+
          
          # Style
          theme_bw()+
          theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                panel.border = element_rect(fill=NA,color="black", linewidth =1.3, linetype="solid"))+
          theme(text = element_text(face = "bold"))+
          scale_color_manual(values =c("#FC4E07","#00AFBB"))+
          
          # Axis
          labs(x=expression(paste("Year")),
               y=expression(paste("Reduce of carbon emmision","  Pg C/year")))+
          theme(axis.title.x = element_text(vjust = 0,size=16),
                axis.title.y = element_text(vjust = 0,size=16))+
          theme(axis.text.x = element_text(size = 12,color="black"),
                axis.text.y = element_text(size = 12,color="black"))
        
        #Future Rs and Rsa (fraction)
        ggplot()+
          #geom_point(data = data_cable,aes(x=Year,y=Rs,color=Scenario))+
          #geom_line(data = data_cable,aes(x=Year,y=Rs,color=Scenario))+
          #geom_point(data = data_cable,aes(x=Year,y=Rsa,color=Scenario))+
          #geom_line(data = data_cable,aes(x=Year,y=Rsa,color=Scenario))+
          
          geom_point(data = data_cable,aes(x=Year,y=fraction,color=Scenario))+
          geom_line(data = data_cable,aes(x=Year,y=fraction,color=Scenario))+
          
          # Style
          theme_bw()+
          theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                panel.border = element_rect(fill=NA,color="black", linewidth =1.3, linetype="solid"))+
          theme(text = element_text(face = "bold"))+
          scale_color_manual(values =c("#FC4E07","#00AFBB"))+
          
          # Axis
          labs(x=expression(paste("Year")),
               y=expression(paste("Reduce of carbon emmision","  Pg C/year")))+
          theme(axis.title.x = element_text(vjust = 0,size=16),
                axis.title.y = element_text(vjust = 0,size=16))+
          theme(axis.text.x = element_text(size = 12,color="black"),
                axis.text.y = element_text(size = 12,color="black"))
        
        
      }
      #================================================
      #================================================
      
      #ORCHIDEE
      #================================================
      ###SSP126----------------------------------------
      data_orchidee126 <- array(dim=c(86,6))
      data_orchidee126[,1] <- 2015:2100
      {
        nc <- nc_open(lsfORCHIDEE[[1,1]])

        lon <- ncvar_get(nc, "lon")
        lon_len <- length(lon)
        lat <- ncvar_get(nc, "lat")
        lat_len <- length(lat)
        
        point_mat <- expand.grid(lon = lon-180, lat = lat)
        coordinates(point_mat) <- ~lon+lat
        proj4string(point_mat) <- proj4string(data_bio)
        
        F_bio <- raster::extract(data_bio, point_mat) 
        F_fT <- raster::extract(fT, point_mat)
        F_area <- raster::extract(area, point_mat)
        F_msap <- raster::extract(msap,point_mat)
        F_rs25a <- raster::extract(rs25_a,point_mat)

        result <- data.frame(point_mat@coords, F_bio = F_bio, F_msap = F_msap, F_area = F_area, F_rs25a=F_rs25a)
        
        pb <- progress_bar$new(total = 86, clear = FALSE)
        days_sum <- 0
        for(y in 2015:2100){

          day_num_year <- days_in_year(y)
          day_start <- days_sum+1
          days_sum <- days_sum+day_num_year

          temp_val <- ncvar_get(nc, "tas",start = c(1,1,day_start),count=c(-1,-1,day_num_year))

          temp_val[temp_val < 278.15] <- NA

          output <- lapply(seq_along(lat), function(j) {
            lapply(seq_along(lon), function(i) {
              
              num <- (j-1)*lon_len+i
              
              tas_values <- temp_val[i, j, ]
              new_time <- (365-sum(is.na(tas_values)))*24*3600 #new_time
              mean_above_5 <- mean(tas_values)-278.15 #new_mgdd5
              new_fT <- 2.2^((mean_above_5-25)/10) #new_fT
              
              F_bio <- result[num,3]
              F_msap <- result[num,4]
              F_area <- result[num,5]
              F_rs25a <- result[num,6]
              
              #not consider thermal acclimation
              new_Rs <- F_rs25a*new_fT*F_msap*new_time*F_area*1.2*10^(-15)  #Pg C/year
              #new_Rs_year <- cellStats(new_Rs,stat='sum')
              
              #consider thermal acclimation
              new_rs25a <- exp(log(rs25)-0.1*(mean_above_5-25))
              new_Rsa <- new_rs25a*new_fT*F_msap*new_time*F_area*1.2*10^(-15)  #Pg C/year
              #new_Rsa_year <- cellStats(new_Rsa,stat='sum')
              
              # recudtion of carbon emission
              reduction <- new_Rsa-new_Rs
              
              #fraction
              fraction <- -reduction/new_Rs
              
              return(list(lon = lon[i]-180,lat = lat[j], mean_above_5 = mean_above_5,time=new_time,new_Rs_year=new_Rs,new_Rsa_year=new_Rsa,
                          reduction=reduction,F_bio=F_bio,F_msap=F_msap,new_fT=new_fT,F_area=F_area,new_time=new_time,fraction=fraction))
            })
          })
          output <- data.frame(do.call(rbind, unlist(output, recursive = FALSE)))
          colnames(output) <- c("lon","lat","mgdd5","time","Rs","Rsa","Reduction","bio","msap","ft","area","time","fraction")
          output_raster <- rasterFromXYZ(output,digits = 4)
          mgdd5_raster <- resample(output_raster$mgdd5,data_bio)
          reduction_raster <- resample(output_raster$Reduction,data_bio)
          #calculate fraction
          Rs_raster <-  resample(output_raster$Rs,data_bio)
          Rsa_raster <-  resample(output_raster$Rsa,data_bio)
          fraction_raster <-  resample(output_raster$fraction,data_bio)
          
          data_orchidee126[(y-2014),2] <- cellStats(mgdd5_raster,stat='mean')
          data_orchidee126[(y-2014),3] <- cellStats(reduction_raster,stat='sum')
          data_orchidee126[(y-2014),4] <- cellStats(Rs_raster,stat='sum')
          data_orchidee126[(y-2014),5] <- cellStats(Rsa_raster,stat='sum')
          data_orchidee126[(y-2014),6] <- cellStats(fraction_raster,stat='mean')
          pb$tick()
          Sys.sleep(1 / 86)
        }
      }
      plot(data_orchidee126[,1],data_orchidee126[,2])#plot mgdd5
      plot(data_orchidee126[,1],data_orchidee126[,3])#plot carbon reduction
      ###SSP585----------------------------------------
      data_orchidee585 <- array(dim=c(86,6))
      data_orchidee585[,1] <- 2015:2100
      {
        nc <- nc_open(lsfORCHIDEE[[2,1]])
        lon <- ncvar_get(nc, "lon")
        lon_len <- length(lon)
        lat <- ncvar_get(nc, "lat")
        lat_len <- length(lat)

        point_mat <- expand.grid(lon = lon-180, lat = lat)
        coordinates(point_mat) <- ~lon+lat
        proj4string(point_mat) <- proj4string(data_bio)

        F_bio <- raster::extract(data_bio, point_mat) 
        F_fT <- raster::extract(fT, point_mat)
        F_area <- raster::extract(area, point_mat)
        F_msap <- raster::extract(msap,point_mat)
        F_rs25a <- raster::extract(rs25_a,point_mat)

        result <- data.frame(point_mat@coords, F_bio = F_bio, F_msap = F_msap, F_area = F_area, F_rs25a=F_rs25a)
        
        pb <- progress_bar$new(total = 86, clear = FALSE)
        days_sum <- 0
        for(y in 2015:2100){
          day_num_year <- days_in_year(y)
          day_start <- days_sum+1
          days_sum <- days_sum+day_num_year

          temp_val <- ncvar_get(nc, "tas",start = c(1,1,day_start),count=c(-1,-1,day_num_year))

          temp_val[temp_val < 278.15] <- NA

          output <- lapply(seq_along(lat), function(j) {
            lapply(seq_along(lon), function(i) {
              
              num <- (j-1)*lon_len+i
              
              tas_values <- temp_val[i, j, ]
              new_time <- (365-sum(is.na(tas_values)))*24*3600 #new_time
              mean_above_5 <- mean(tas_values)-278.15 #new_mgdd5
              new_fT <- 2.2^((mean_above_5-25)/10) #new_fT
              
              F_bio <- result[num,3]
              F_msap <- result[num,4]
              F_area <- result[num,5]
              F_rs25a <- result[num,6]
              
              #not consider thermal acclimation
              new_Rs <- F_rs25a*new_fT*F_msap*new_time*F_area*1.2*10^(-15)  #Pg C/year
              #new_Rs_year <- cellStats(new_Rs,stat='sum')
              
              #consider thermal acclimation
              new_rs25a <- exp(log(rs25)-0.1*(mean_above_5-25))
              new_Rsa <- new_rs25a*new_fT*F_msap*new_time*F_area*1.2*10^(-15)  #Pg C/year
              #new_Rsa_year <- cellStats(new_Rsa,stat='sum')
              
              # recudtion of carbon emission
              reduction <- new_Rsa-new_Rs
              
              #fraction
              fraction <- -reduction/new_Rs
              return(list(lon = lon[i]-180,lat = lat[j], mean_above_5 = mean_above_5,time=new_time,new_Rs_year=new_Rs,new_Rsa_year=new_Rsa,
                          reduction=reduction,F_bio=F_bio,F_msap=F_msap,new_fT=new_fT,F_area=F_area,new_time=new_time,fraction=fraction))
            })
          })
          output <- data.frame(do.call(rbind, unlist(output, recursive = FALSE)))
          colnames(output) <- c("lon","lat","mgdd5","time","Rs","Rsa","Reduction","bio","msap","ft","area","time","fraction")
          output_raster <- rasterFromXYZ(output,digits = 4)
          mgdd5_raster <- resample(output_raster$mgdd5,data_bio)
          reduction_raster <- resample(output_raster$Reduction,data_bio)
          #calculate fraction
          Rs_raster <-  resample(output_raster$Rs,data_bio)
          Rsa_raster <-  resample(output_raster$Rsa,data_bio)
          fraction_raster <-  resample(output_raster$fraction,data_bio)
          
          data_orchidee585[(y-2014),2] <- cellStats(mgdd5_raster,stat='mean')
          data_orchidee585[(y-2014),3] <- cellStats(reduction_raster,stat='sum')
          data_orchidee585[(y-2014),4] <- cellStats(Rs_raster,stat='sum')
          data_orchidee585[(y-2014),5] <- cellStats(Rsa_raster,stat='sum')
          data_orchidee585[(y-2014),6] <- cellStats(fraction_raster,stat='mean')
          pb$tick()
          Sys.sleep(1 / 86)
        }
      }
      plot(data_orchidee585[,1],data_orchidee585[,2])#plot mgdd5
      plot(data_orchidee585[,1],data_orchidee585[,3])#plot carbon reduction
      #Plot single model
      data_orchidee_tem <- cbind(data_orchidee126[,1:2],data_orchidee585[,2])
      data_orchidee_red <- cbind(data_orchidee126[,c(1,3)],data_orchidee585[,3])
      data_orchidee_Rs <- cbind(data_orchidee126[,c(1,4)],data_orchidee585[,4])
      data_orchidee_Rsa <- cbind(data_orchidee126[,c(1,5)],data_orchidee585[,5])
      data_orchidee_fraction <- cbind(data_orchidee126[,c(1,6)],data_orchidee585[,6])
      
      
      colnames(data_orchidee_tem) <- c("Year","SSP126","SSP585")
      colnames(data_orchidee_red) <- c("Year","SSP126","SSP585")
      colnames(data_orchidee_Rs) <- c("Year","SSP126","SSP585")
      colnames(data_orchidee_Rsa) <- c("Year","SSP126","SSP585")
      colnames(data_orchidee_fraction) <- c("Year","SSP126","SSP585")
      
      data_orchidee_tem <- data.frame(data_orchidee_tem)
      data_orchidee_red <- data.frame(data_orchidee_red)
      data_orchidee_Rs <- data.frame(data_orchidee_Rs)
      data_orchidee_Rsa <- data.frame(data_orchidee_Rsa)
      data_orchidee_fraction <- data.frame(data_orchidee_fraction)
      
      data_orchidee_m1 <- reshape2::melt(id="Year",data_orchidee_tem)
      data_orchidee_m2 <- reshape2::melt(id="Year",data_orchidee_red)
      data_orchidee_m3 <- reshape2::melt(id="Year",data_orchidee_Rs)
      data_orchidee_m4 <- reshape2::melt(id="Year",data_orchidee_Rsa)
      data_orchidee_m5 <- reshape2::melt(id="Year",data_orchidee_fraction)
      
      data_orchidee <- cbind(data_orchidee_m1[,],data_orchidee_m2[,3],data_orchidee_m3[,3],data_orchidee_m4[,3],data_orchidee_m5[,3])
      colnames(data_orchidee) <- c("Year","Scenario","T5","reduction","Rs","Rsa","fraction")
      #plot
      {
        ggplot()+
          geom_point(data = data_orchidee,aes(x=Year,y=reduction,color=Scenario))+
          geom_line(data = data_orchidee,aes(x=Year,y=reduction,color=Scenario))+
          
          # Style
          theme_bw()+
          theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                panel.border = element_rect(fill=NA,color="black", linewidth=1.3, linetype="solid"))+
          theme(text = element_text(face = "bold"))+
          scale_color_manual(values =c("#FC4E07","#00AFBB"))+
          
          # Axis
          labs(x=expression(paste("Year")),
               y=expression(paste("Reduce of carbon emmision","  Pg C/year")))+
          theme(axis.title.x = element_text(vjust = 0,size=16),
                axis.title.y = element_text(vjust = 0,size=16))+
          theme(axis.text.x = element_text(size = 12,color="black"),
                axis.text.y = element_text(size = 12,color="black"))
        
      }
      
      #JULES
      #================================================
      ###SSP126----------------------------------------
      data_jules126 <- array(dim=c(86,6))
      data_jules126[,1] <- 2015:2100
      #2015-2049
      {
        nc <- nc_open(lsfJULES[[1,1]])
        lon <- ncvar_get(nc, "lon")
        lon_len <- length(lon)
        lat <- ncvar_get(nc, "lat")
        lat_len <- length(lat)

        point_mat <- expand.grid(lon = lon-180, lat = lat)
        coordinates(point_mat) <- ~lon+lat
        proj4string(point_mat) <- proj4string(data_bio)

        F_bio <- raster::extract(data_bio, point_mat) 
        F_fT <- raster::extract(fT, point_mat)
        F_area <- raster::extract(area, point_mat)
        F_msap <- raster::extract(msap,point_mat)
        F_rs25a <- raster::extract(rs25_a,point_mat)

        result <- data.frame(point_mat@coords, F_bio = F_bio, F_msap = F_msap, F_area = F_area, F_rs25a=F_rs25a)
        
        pb <- progress_bar$new(total = 35, clear = FALSE)
        for(y in 2015:2049){

          day_start <- (y-2015)*360+1
          
          temp_val <- ncvar_get(nc, "tas",start = c(1,1,day_start),count=c(-1,-1,360))
          temp_val[temp_val < 278.15] <- NA

          output <- lapply(seq_along(lat), function(j) {
            lapply(seq_along(lon), function(i) {
              
              num <- (j-1)*lon_len+i
              
              tas_values <- temp_val[i, j, ]
              new_time <- (365-sum(is.na(tas_values)))*24*3600 #new_time
              mean_above_5 <- mean(tas_values)-278.15 #new_mgdd5
              new_fT <- 2.2^((mean_above_5-25)/10) #new_fT
              
              F_bio <- result[num,3]
              F_msap <- result[num,4]
              F_area <- result[num,5]
              F_rs25a <- result[num,6]
              
              #not consider thermal acclimation
              new_Rs <- F_rs25a*new_fT*F_msap*new_time*F_area*1.2*10^(-15)  #Pg C/year
              #new_Rs_year <- cellStats(new_Rs,stat='sum')
              
              #consider thermal acclimation
              new_rs25a <- exp(log(rs25)-0.1*(mean_above_5-25))
              new_Rsa <- new_rs25a*new_fT*F_msap*new_time*F_area*1.2*10^(-15)  #Pg C/year
              #new_Rsa_year <- cellStats(new_Rsa,stat='sum')
              
              # recudtion of carbon emission
              reduction <- new_Rsa-new_Rs
              
              #fraction
              fraction <- -reduction/new_Rs
              return(list(lon = lon[i]-180,lat = lat[j], mean_above_5 = mean_above_5,time=new_time,new_Rs_year=new_Rs,new_Rsa_year=new_Rsa,
                          reduction=reduction,F_bio=F_bio,F_msap=F_msap,new_fT=new_fT,F_area=F_area,new_time=new_time,fraction=fraction))
            })
          })
          output <- data.frame(do.call(rbind, unlist(output, recursive = FALSE)))
          colnames(output) <- c("lon","lat","mgdd5","time","Rs","Rsa","Reduction","bio","msap","ft","area","time","fraction")
          output_raster <- rasterFromXYZ(output)
          mgdd5_raster <- resample(output_raster$mgdd5,data_bio)
          reduction_raster <- resample(output_raster$Reduction,data_bio)
          #calculate fraction
          Rs_raster <-  resample(output_raster$Rs,data_bio)
          Rsa_raster <-  resample(output_raster$Rsa,data_bio)
          fraction_raster <-  resample(output_raster$fraction,data_bio)
          
          data_jules126[(y-2014),2] <- cellStats(mgdd5_raster,stat='mean')
          data_jules126[(y-2014),3] <- cellStats(reduction_raster,stat='sum')
          data_jules126[(y-2014),4] <- cellStats(Rs_raster,stat='sum')
          data_jules126[(y-2014),5] <- cellStats(Rsa_raster,stat='sum')
          data_jules126[(y-2014),6] <- cellStats(fraction_raster,stat='mean')
          
          pb$tick()
          Sys.sleep(1 / 35)
        }
      }
      #2050-2100
      {
        nc <- nc_open(lsfJULES[[2,1]])
        lon <- ncvar_get(nc, "lon")
        lon_len <- length(lon)
        lat <- ncvar_get(nc, "lat")
        lat_len <- length(lat)

        point_mat <- expand.grid(lon = lon-180, lat = lat)
        coordinates(point_mat) <- ~lon+lat
        proj4string(point_mat) <- proj4string(data_bio)

        F_bio <- raster::extract(data_bio, point_mat) 
        F_fT <- raster::extract(fT, point_mat)
        F_area <- raster::extract(area, point_mat)
        F_msap <- raster::extract(msap,point_mat)
        F_rs25a <- raster::extract(rs25_a,point_mat)

        result <- data.frame(point_mat@coords, F_bio = F_bio, F_msap = F_msap, F_area = F_area, F_rs25a=F_rs25a)
        
        pb <- progress_bar$new(total = 51, clear = FALSE)
        for(y in 2050:2100){
          day_start <- (y-2050)*360+1

          temp_val <- ncvar_get(nc, "tas",start = c(1,1,day_start),count=c(-1,-1,360))

          temp_val[temp_val < 278.15] <- NA

          output <- lapply(seq_along(lat), function(j) {
            lapply(seq_along(lon), function(i) {
              
              num <- (j-1)*lon_len+i
              
              tas_values <- temp_val[i, j, ]
              new_time <- (365-sum(is.na(tas_values)))*24*3600 #new_time
              mean_above_5 <- mean(tas_values)-278.15 #new_mgdd5
              new_fT <- 2.2^((mean_above_5-25)/10) #new_fT
              
              F_bio <- result[num,3]
              F_msap <- result[num,4]
              F_area <- result[num,5]
              F_rs25a <- result[num,6]
              
              #not consider thermal acclimation
              new_Rs <- F_rs25a*new_fT*F_msap*new_time*F_area*1.2*10^(-15)  #Pg C/year
              #new_Rs_year <- cellStats(new_Rs,stat='sum')
              
              #consider thermal acclimation
              new_rs25a <- exp(log(rs25)-0.1*(mean_above_5-25))
              new_Rsa <- new_rs25a*new_fT*F_msap*new_time*F_area*1.2*10^(-15)  #Pg C/year
              #new_Rsa_year <- cellStats(new_Rsa,stat='sum')
              
              # recudtion of carbon emission
              reduction <- new_Rsa-new_Rs
              
              #fraction
              fraction <- -reduction/new_Rs
              return(list(lon = lon[i]-180,lat = lat[j], mean_above_5 = mean_above_5,time=new_time,new_Rs_year=new_Rs,new_Rsa_year=new_Rsa,
                          reduction=reduction,F_bio=F_bio,F_msap=F_msap,new_fT=new_fT,F_area=F_area,new_time=new_time,fraction=fraction))
            })
          })
          output <- data.frame(do.call(rbind, unlist(output, recursive = FALSE)))
          colnames(output) <- c("lon","lat","mgdd5","time","Rs","Rsa","Reduction","bio","msap","ft","area","time","fraction")
          output_raster <- rasterFromXYZ(output)
          mgdd5_raster <- resample(output_raster$mgdd5,data_bio)
          reduction_raster <- resample(output_raster$Reduction,data_bio)
          #calculate fraction
          Rs_raster <-  resample(output_raster$Rs,data_bio)
          Rsa_raster <-  resample(output_raster$Rsa,data_bio)
          fraction_raster <-  resample(output_raster$fraction,data_bio)
          
          data_jules126[(y-2014),2] <- cellStats(mgdd5_raster,stat='mean')
          data_jules126[(y-2014),3] <- cellStats(reduction_raster,stat='sum')
          data_jules126[(y-2014),4] <- cellStats(Rs_raster,stat='sum')
          data_jules126[(y-2014),5] <- cellStats(Rsa_raster,stat='sum')
          data_jules126[(y-2014),6] <- cellStats(fraction_raster,stat='mean')
          pb$tick()
          Sys.sleep(1 / 51)
        }
      }
      plot(data_jules126[,1],data_jules126[,2])#plot mgdd5
      plot(data_jules126[,1],data_jules126[,3])#plot carbon reduction
      ###SSP585----------------------------------------
      data_jules585 <- array(dim=c(86,6))
      data_jules585[,1] <- 2015:2100
      #2015-2049
      {
        nc <- nc_open(lsfJULES[[3,1]])
        lon <- ncvar_get(nc, "lon")
        lon_len <- length(lon)
        lat <- ncvar_get(nc, "lat")
        lat_len <- length(lat)

        point_mat <- expand.grid(lon = lon-180, lat = lat)
        coordinates(point_mat) <- ~lon+lat
        proj4string(point_mat) <- proj4string(data_bio)
        
        F_bio <- raster::extract(data_bio, point_mat) 
        F_fT <- raster::extract(fT, point_mat)
        F_area <- raster::extract(area, point_mat)
        F_msap <- raster::extract(msap,point_mat)
        F_rs25a <- raster::extract(rs25_a,point_mat)

        result <- data.frame(point_mat@coords, F_bio = F_bio, F_msap = F_msap, F_area = F_area, F_rs25a=F_rs25a)
        
        pb <- progress_bar$new(total = 35, clear = FALSE)
        for(y in 2015:2049){
          day_start <- (y-2015)*360+1

          temp_val <- ncvar_get(nc, "tas",start = c(1,1,day_start),count=c(-1,-1,360))

          temp_val[temp_val < 278.15] <- NA

          output <- lapply(seq_along(lat), function(j) {
            lapply(seq_along(lon), function(i) {
              
              num <- (j-1)*lon_len+i
              
              tas_values <- temp_val[i, j, ]
              new_time <- (365-sum(is.na(tas_values)))*24*3600 #new_time
              mean_above_5 <- mean(tas_values)-278.15 #new_mgdd5
              new_fT <- 2.2^((mean_above_5-25)/10) #new_fT
              
              F_bio <- result[num,3]
              F_msap <- result[num,4]
              F_area <- result[num,5]
              F_rs25a <- result[num,6]
              
              #not consider thermal acclimation
              new_Rs <- F_rs25a*new_fT*F_msap*new_time*F_area*1.2*10^(-15)  #Pg C/year
              #new_Rs_year <- cellStats(new_Rs,stat='sum')
              
              #consider thermal acclimation
              new_rs25a <- exp(log(rs25)-0.1*(mean_above_5-25))
              new_Rsa <- new_rs25a*new_fT*F_msap*new_time*F_area*1.2*10^(-15)  #Pg C/year
              #new_Rsa_year <- cellStats(new_Rsa,stat='sum')
              
              # recudtion of carbon emission
              reduction <- new_Rsa-new_Rs
              
              #fraction
              fraction <- -reduction/new_Rs
              return(list(lon = lon[i]-180,lat = lat[j], mean_above_5 = mean_above_5,time=new_time,new_Rs_year=new_Rs,new_Rsa_year=new_Rsa,
                          reduction=reduction,F_bio=F_bio,F_msap=F_msap,new_fT=new_fT,F_area=F_area,new_time=new_time,fraction=fraction))
            })
          })
          output <- data.frame(do.call(rbind, unlist(output, recursive = FALSE)))
          colnames(output) <- c("lon","lat","mgdd5","time","Rs","Rsa","Reduction","bio","msap","ft","area","time","fraction")
          output_raster <- rasterFromXYZ(output)
          mgdd5_raster <- resample(output_raster$mgdd5,data_bio)
          reduction_raster <- resample(output_raster$Reduction,data_bio)
          #calculate fraction
          Rs_raster <-  resample(output_raster$Rs,data_bio)
          Rsa_raster <-  resample(output_raster$Rsa,data_bio)
          fraction_raster <-  resample(output_raster$fraction,data_bio)
          
          data_jules585[(y-2014),2] <- cellStats(mgdd5_raster,stat='mean')
          data_jules585[(y-2014),3] <- cellStats(reduction_raster,stat='sum')
          data_jules585[(y-2014),4] <- cellStats(Rs_raster,stat='sum')
          data_jules585[(y-2014),5] <- cellStats(Rsa_raster,stat='sum')
          data_jules585[(y-2014),6] <- cellStats(fraction_raster,stat='mean')
          
          pb$tick()
          Sys.sleep(1 / 35)
        }
      }
      #2050-2100
      {
        nc <- nc_open(lsfJULES[[4,1]])
        lon <- ncvar_get(nc, "lon")
        lon_len <- length(lon)
        lat <- ncvar_get(nc, "lat")
        lat_len <- length(lat)

        point_mat <- expand.grid(lon = lon-180, lat = lat)
        coordinates(point_mat) <- ~lon+lat
        proj4string(point_mat) <- proj4string(data_bio)
        
        F_bio <- raster::extract(data_bio, point_mat) 
        F_fT <- raster::extract(fT, point_mat)
        F_area <- raster::extract(area, point_mat)
        F_msap <- raster::extract(msap,point_mat)
        F_rs25a <- raster::extract(rs25_a,point_mat)
        
        result <- data.frame(point_mat@coords, F_bio = F_bio, F_msap = F_msap, F_area = F_area, F_rs25a=F_rs25a)
        
        pb <- progress_bar$new(total = 51, clear = FALSE)
        for(y in 2050:2100){
          day_start <- (y-2050)*360+1
          
          temp_val <- ncvar_get(nc, "tas",start = c(1,1,day_start),count=c(-1,-1,360))
          temp_val[temp_val < 278.15] <- NA

          output <- lapply(seq_along(lat), function(j) {
            lapply(seq_along(lon), function(i) {
              
              num <- (j-1)*lon_len+i
              
              tas_values <- temp_val[i, j, ]
              new_time <- (365-sum(is.na(tas_values)))*24*3600 #new_time
              mean_above_5 <- mean(tas_values)-278.15 #new_mgdd5
              new_fT <- 2.2^((mean_above_5-25)/10) #new_fT
              
              F_bio <- result[num,3]
              F_msap <- result[num,4]
              F_area <- result[num,5]
              F_rs25a <- result[num,6]
              
              #not consider thermal acclimation
              new_Rs <- F_rs25a*new_fT*F_msap*new_time*F_area*1.2*10^(-15)  #Pg C/year
              #new_Rs_year <- cellStats(new_Rs,stat='sum')
              
              #consider thermal acclimation
              new_rs25a <- exp(log(rs25)-0.1*(mean_above_5-25))
              new_Rsa <- new_rs25a*new_fT*F_msap*new_time*F_area*1.2*10^(-15)  #Pg C/year
              #new_Rsa_year <- cellStats(new_Rsa,stat='sum')
              
              # recudtion of carbon emission
              reduction <- new_Rsa-new_Rs
              
              #fraction
              fraction <- -reduction/new_Rs
              return(list(lon = lon[i]-180,lat = lat[j], mean_above_5 = mean_above_5,time=new_time,new_Rs_year=new_Rs,new_Rsa_year=new_Rsa,
                          reduction=reduction,F_bio=F_bio,F_msap=F_msap,new_fT=new_fT,F_area=F_area,new_time=new_time,fraction=fraction))
            })
          })
          output <- data.frame(do.call(rbind, unlist(output, recursive = FALSE)))
          colnames(output) <- c("lon","lat","mgdd5","time","Rs","Rsa","Reduction","bio","msap","ft","area","time","fraction")
          output_raster <- rasterFromXYZ(output)
          mgdd5_raster <- resample(output_raster$mgdd5,data_bio)
          reduction_raster <- resample(output_raster$Reduction,data_bio)
          #calculate fraction
          Rs_raster <-  resample(output_raster$Rs,data_bio)
          Rsa_raster <-  resample(output_raster$Rsa,data_bio)
          fraction_raster <-  resample(output_raster$fraction,data_bio)
          
          data_jules585[(y-2014),2] <- cellStats(mgdd5_raster,stat='mean')
          data_jules585[(y-2014),3] <- cellStats(reduction_raster,stat='sum')
          data_jules585[(y-2014),4] <- cellStats(Rs_raster,stat='sum')
          data_jules585[(y-2014),5] <- cellStats(Rsa_raster,stat='sum')
          data_jules585[(y-2014),6] <- cellStats(fraction_raster,stat='mean')
          pb$tick()
          Sys.sleep(1 / 51)
        }
      }
      plot(data_jules585[,1],data_jules585[,2])#plot mgdd5
      plot(data_jules585[,1],data_jules585[,3])#plot carbon reduction
      
      #Plot single model
      data_jules_tem <- cbind(data_jules126[,1:2],data_jules585[,2])
      data_jules_red <- cbind(data_jules126[,c(1,3)],data_jules585[,3])
      data_jules_Rs <- cbind(data_jules126[,c(1,4)],data_jules585[,4])
      data_jules_Rsa <- cbind(data_jules126[,c(1,5)],data_jules585[,5])
      data_jules_fraction <- cbind(data_jules126[,c(1,6)],data_jules585[,6])
      
      
      colnames(data_jules_tem) <- c("Year","SSP126","SSP585")
      colnames(data_jules_red) <- c("Year","SSP126","SSP585")
      colnames(data_jules_Rs) <- c("Year","SSP126","SSP585")
      colnames(data_jules_Rsa) <- c("Year","SSP126","SSP585")
      colnames(data_jules_fraction) <- c("Year","SSP126","SSP585")
      
      
      data_jules_tem <- data.frame(data_jules_tem)
      data_jules_red <- data.frame(data_jules_red)
      data_jules_Rs <- data.frame(data_jules_Rs)
      data_jules_Rsa <- data.frame(data_jules_Rsa)
      data_jules_fraction <- data.frame(data_jules_fraction)
      
      
      data_jules_m1 <- reshape2::melt(id="Year",data_jules_tem)
      data_jules_m2 <- reshape2::melt(id="Year",data_jules_red)
      data_jules_m3 <- reshape2::melt(id="Year",data_jules_Rs)
      data_jules_m4 <- reshape2::melt(id="Year",data_jules_Rsa)
      data_jules_m5 <- reshape2::melt(id="Year",data_jules_fraction)
      
      data_jules <- cbind(data_jules_m1[,],data_jules_m2[,3],data_jules_m3[,3],data_jules_m4[,3],data_jules_m5[,3])
      colnames(data_jules) <- c("Year","Scenario","T5","reduction","Rs","Rsa","fraction")
      #plot
      {
        ggplot()+
          geom_point(data = data_jules,aes(x=Year,y=reduction,color=Scenario))+
          geom_line(data = data_jules,aes(x=Year,y=reduction,color=Scenario))+
          
          # Style
          theme_bw()+
          theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                panel.border = element_rect(fill=NA,color="black", linewidth=1.3, linetype="solid"))+
          theme(text = element_text(face = "bold"))+
          scale_color_manual(values =c("#FC4E07","#00AFBB"))+
          
          # Axis
          labs(x=expression(paste("Year")),
               y=expression(paste("Reduce of carbon emmision","  Pg C/year")))+
          theme(axis.title.x = element_text(vjust = 0,size=16),
                axis.title.y = element_text(vjust = 0,size=16))+
          theme(axis.text.x = element_text(size = 12,color="black"),
                axis.text.y = element_text(size = 12,color="black"))
        
      }
      
      
      ######
      #CLM5
      #================================================
      ###SSP126----------------------------------------
      data_clm126 <- array(dim=c(86,6))
      data_clm126[,1] <- 2015:2100
      #2015-2094
      for(ls in 1:8){
        nc <- nc_open(lsfCLM5[[ls,1]])
        # 获取经度和纬度的长度和值
        lon <- ncvar_get(nc, "lon")
        lon_len <- length(lon)
        lat <- ncvar_get(nc, "lat")
        lat_len <- length(lat)

        point_mat <- expand.grid(lon = lon-180, lat = lat)
        coordinates(point_mat) <- ~lon+lat
        proj4string(point_mat) <- proj4string(data_bio)
        
        F_bio <- raster::extract(data_bio, point_mat) 
        F_fT <- raster::extract(fT, point_mat)
        F_area <- raster::extract(area, point_mat)
        F_msap <- raster::extract(msap,point_mat)
        F_rs25a <- raster::extract(rs25_a,point_mat)
        
        result <- data.frame(point_mat@coords, F_bio = F_bio, F_msap = F_msap, F_area = F_area, F_rs25a=F_rs25a)
        
        pb <- progress_bar$new(total = 10, clear = FALSE)
        start_year <- (ls-1)*10+2015
        end_year <- 2015+ls*10-1
        for(y in start_year:end_year){
          day_start <- (y-start_year)*365+1
          
          temp_val <- ncvar_get(nc, "tas",start = c(1,1,day_start),count=c(-1,-1,365))

          temp_val[temp_val < 278.15] <- NA

          output <- lapply(seq_along(lat), function(j) {
            lapply(seq_along(lon), function(i) {
              
              num <- (j-1)*lon_len+i
              
              tas_values <- temp_val[i, j, ]
              new_time <- (365-sum(is.na(tas_values)))*24*3600 #new_time
              mean_above_5 <- mean(tas_values)-278.15 #new_mgdd5
              new_fT <- 2.2^((mean_above_5-25)/10) #new_fT
              
              F_bio <- result[num,3]
              F_msap <- result[num,4]
              F_area <- result[num,5]
              F_rs25a <- result[num,6]
              
              #not consider thermal acclimation
              new_Rs <- F_rs25a*new_fT*F_msap*new_time*F_area*1.2*10^(-15)  #Pg C/year
              #new_Rs_year <- cellStats(new_Rs,stat='sum')
              
              #consider thermal acclimation
              new_rs25a <- exp(log(rs25)-0.1*(mean_above_5-25))
              new_Rsa <- new_rs25a*new_fT*F_msap*new_time*F_area*1.2*10^(-15)  #Pg C/year
              #new_Rsa_year <- cellStats(new_Rsa,stat='sum')
              
              # recudtion of carbon emission
              reduction <- new_Rsa-new_Rs
              
              #fraction
              fraction <- -reduction/new_Rs
              return(list(lon = lon[i]-180,lat = lat[j], mean_above_5 = mean_above_5,time=new_time,new_Rs_year=new_Rs,new_Rsa_year=new_Rsa,
                          reduction=reduction,F_bio=F_bio,F_msap=F_msap,new_fT=new_fT,F_area=F_area,new_time=new_time,fraction=fraction))
            })
          })
          output <- data.frame(do.call(rbind, unlist(output, recursive = FALSE)))
          colnames(output) <- c("lon","lat","mgdd5","time","Rs","Rsa","Reduction","bio","msap","ft","area","time","fraction")
          output_raster <- rasterFromXYZ(output)
          mgdd5_raster <- resample(output_raster$mgdd5,data_bio)
          reduction_raster <- resample(output_raster$Reduction,data_bio)
          #calculate fraction
          Rs_raster <-  resample(output_raster$Rs,data_bio)
          Rsa_raster <-  resample(output_raster$Rsa,data_bio)
          fraction_raster <-  resample(output_raster$fraction,data_bio)
          
          data_clm126[(y-2014),2] <- cellStats(mgdd5_raster,stat='mean')
          data_clm126[(y-2014),3] <- cellStats(reduction_raster,stat='sum')
          data_clm126[(y-2014),4] <- cellStats(Rs_raster,stat='sum')
          data_clm126[(y-2014),5] <- cellStats(Rsa_raster,stat='sum')
          data_clm126[(y-2014),6] <- cellStats(fraction_raster,stat='mean')
          pb$tick()
          Sys.sleep(1 / 10)
        }
      }
      #2095-2100
      {
        nc <- nc_open(lsfCLM5[[9,1]])
        lon <- ncvar_get(nc, "lon")
        lon_len <- length(lon)
        lat <- ncvar_get(nc, "lat")
        lat_len <- length(lat)

        point_mat <- expand.grid(lon = lon-180, lat = lat)
        coordinates(point_mat) <- ~lon+lat
        proj4string(point_mat) <- proj4string(data_bio)

        F_bio <- raster::extract(data_bio, point_mat) 
        F_fT <- raster::extract(fT, point_mat)
        F_area <- raster::extract(area, point_mat)
        F_msap <- raster::extract(msap,point_mat)
        F_rs25a <- raster::extract(rs25_a,point_mat)

        result <- data.frame(point_mat@coords, F_bio = F_bio, F_msap = F_msap, F_area = F_area, F_rs25a=F_rs25a)
        
        pb <- progress_bar$new(total = 6, clear = FALSE)
        for(y in 2095:2100){
          day_start <- (y-2095)*365+1

          temp_val <- ncvar_get(nc, "tas",start = c(1,1,day_start),count=c(-1,-1,365))

          temp_val[temp_val < 278.15] <- NA

          output <- lapply(seq_along(lat), function(j) {
            lapply(seq_along(lon), function(i) {
              
              num <- (j-1)*lon_len+i
              
              tas_values <- temp_val[i, j, ]
              new_time <- (365-sum(is.na(tas_values)))*24*3600 #new_time
              mean_above_5 <- mean(tas_values)-278.15 #new_mgdd5
              new_fT <- 2.2^((mean_above_5-25)/10) #new_fT
              
              F_bio <- result[num,3]
              F_msap <- result[num,4]
              F_area <- result[num,5]
              F_rs25a <- result[num,6]
              
              #not consider thermal acclimation
              new_Rs <- F_rs25a*new_fT*F_msap*new_time*F_area*1.2*10^(-15)  #Pg C/year
              #new_Rs_year <- cellStats(new_Rs,stat='sum')
              
              #consider thermal acclimation
              new_rs25a <- exp(log(rs25)-0.1*(mean_above_5-25))
              new_Rsa <- new_rs25a*new_fT*F_msap*new_time*F_area*1.2*10^(-15)  #Pg C/year
              #new_Rsa_year <- cellStats(new_Rsa,stat='sum')
              
              # recudtion of carbon emission
              reduction <- new_Rsa-new_Rs
              
              #fraction
              fraction <- -reduction/new_Rs
              return(list(lon = lon[i]-180,lat = lat[j], mean_above_5 = mean_above_5,time=new_time,new_Rs_year=new_Rs,new_Rsa_year=new_Rsa,
                          reduction=reduction,F_bio=F_bio,F_msap=F_msap,new_fT=new_fT,F_area=F_area,new_time=new_time,fraction=fraction))
            })
          })
          output <- data.frame(do.call(rbind, unlist(output, recursive = FALSE)))
          colnames(output) <- c("lon","lat","mgdd5","time","Rs","Rsa","Reduction","bio","msap","ft","area","time","fraction")
          output_raster <- rasterFromXYZ(output)
          mgdd5_raster <- resample(output_raster$mgdd5,data_bio)
          reduction_raster <- resample(output_raster$Reduction,data_bio)
          #calculate fraction
          Rs_raster <-  resample(output_raster$Rs,data_bio)
          Rsa_raster <-  resample(output_raster$Rsa,data_bio)
          fraction_raster <-  resample(output_raster$fraction,data_bio)
          
          data_clm126[(y-2014),2] <- cellStats(mgdd5_raster,stat='mean')
          data_clm126[(y-2014),3] <- cellStats(reduction_raster,stat='sum')
          data_clm126[(y-2014),4] <- cellStats(Rs_raster,stat='sum')
          data_clm126[(y-2014),5] <- cellStats(Rsa_raster,stat='sum')
          data_clm126[(y-2014),6] <- cellStats(fraction_raster,stat='mean')
          pb$tick()
          Sys.sleep(1 / 6)
        }
      }
      plot(data_clm126[,1],data_clm126[,2])#plot mgdd5
      plot(data_clm126[,1],data_clm126[,3])#plot carbon reduction
      ###SSP585----------------------------------------
      data_clm585 <- array(dim=c(86,6))
      data_clm585[,1] <- 2015:2100
      #2015-2094
      for(ls in 1:8){
        nc <- nc_open(lsfCLM5[[(ls+9),1]])
        lon <- ncvar_get(nc, "lon")
        lon_len <- length(lon)
        lat <- ncvar_get(nc, "lat")
        lat_len <- length(lat)

        point_mat <- expand.grid(lon = lon-180, lat = lat)
        coordinates(point_mat) <- ~lon+lat
        proj4string(point_mat) <- proj4string(data_bio)
        
        F_bio <- raster::extract(data_bio, point_mat) 
        F_fT <- raster::extract(fT, point_mat)
        F_area <- raster::extract(area, point_mat)
        F_msap <- raster::extract(msap,point_mat)
        F_rs25a <- raster::extract(rs25_a,point_mat)

        result <- data.frame(point_mat@coords, F_bio = F_bio, F_msap = F_msap, F_area = F_area, F_rs25a=F_rs25a)
        
        pb <- progress_bar$new(total = 10, clear = FALSE)
        start_year <- (ls-1)*10+2015
        end_year <- 2015+ls*10-1
        for(y in start_year:end_year){
          day_start <- (y-start_year)*365+1

          temp_val <- ncvar_get(nc, "tas",start = c(1,1,day_start),count=c(-1,-1,365))

          temp_val[temp_val < 278.15] <- NA

          output <- lapply(seq_along(lat), function(j) {
            lapply(seq_along(lon), function(i) {
              
              num <- (j-1)*lon_len+i
              
              tas_values <- temp_val[i, j, ]
              new_time <- (365-sum(is.na(tas_values)))*24*3600 #new_time
              mean_above_5 <- mean(tas_values)-278.15 #new_mgdd5
              new_fT <- 2.2^((mean_above_5-25)/10) #new_fT
              
              F_bio <- result[num,3]
              F_msap <- result[num,4]
              F_area <- result[num,5]
              F_rs25a <- result[num,6]
              
              #not consider thermal acclimation
              new_Rs <- F_rs25a*new_fT*F_msap*new_time*F_area*1.2*10^(-15)  #Pg C/year
              #new_Rs_year <- cellStats(new_Rs,stat='sum')
              
              #consider thermal acclimation
              new_rs25a <- exp(log(rs25)-0.1*(mean_above_5-25))
              new_Rsa <- new_rs25a*new_fT*F_msap*new_time*F_area*1.2*10^(-15)  #Pg C/year
              #new_Rsa_year <- cellStats(new_Rsa,stat='sum')
              
              # recudtion of carbon emission
              reduction <- new_Rsa-new_Rs
              
              #fraction
              fraction <- -reduction/new_Rs
              return(list(lon = lon[i]-180,lat = lat[j], mean_above_5 = mean_above_5,time=new_time,new_Rs_year=new_Rs,new_Rsa_year=new_Rsa,
                          reduction=reduction,F_bio=F_bio,F_msap=F_msap,new_fT=new_fT,F_area=F_area,new_time=new_time,fraction=fraction))
            })
          })
          output <- data.frame(do.call(rbind, unlist(output, recursive = FALSE)))
          colnames(output) <- c("lon","lat","mgdd5","time","Rs","Rsa","Reduction","bio","msap","ft","area","time","fraction")
          output_raster <- rasterFromXYZ(output)
          mgdd5_raster <- resample(output_raster$mgdd5,data_bio)
          reduction_raster <- resample(output_raster$Reduction,data_bio)
          #calculate fraction
          Rs_raster <-  resample(output_raster$Rs,data_bio)
          Rsa_raster <-  resample(output_raster$Rsa,data_bio)
          fraction_raster <-  resample(output_raster$fraction,data_bio)
          
          data_clm585[(y-2014),2] <- cellStats(mgdd5_raster,stat='mean')
          data_clm585[(y-2014),3] <- cellStats(reduction_raster,stat='sum')
          data_clm585[(y-2014),4] <- cellStats(Rs_raster,stat='sum')
          data_clm585[(y-2014),5] <- cellStats(Rsa_raster,stat='sum')
          data_clm585[(y-2014),6] <- cellStats(fraction_raster,stat='mean')
          pb$tick()
          Sys.sleep(1 / 10)
        }
      }
      #2095-2100
      {
        nc <- nc_open(lsfCLM5[[18,1]])
        lon <- ncvar_get(nc, "lon")
        lon_len <- length(lon)
        lat <- ncvar_get(nc, "lat")
        lat_len <- length(lat)
        
        point_mat <- expand.grid(lon = lon-180, lat = lat)
        coordinates(point_mat) <- ~lon+lat
        proj4string(point_mat) <- proj4string(data_bio)
        
        F_bio <- raster::extract(data_bio, point_mat) 
        F_fT <- raster::extract(fT, point_mat)
        F_area <- raster::extract(area, point_mat)
        F_msap <- raster::extract(msap,point_mat)
        F_rs25a <- raster::extract(rs25_a,point_mat)

        result <- data.frame(point_mat@coords, F_bio = F_bio, F_msap = F_msap, F_area = F_area, F_rs25a=F_rs25a)
        
        pb <- progress_bar$new(total = 6, clear = FALSE)
        for(y in 2095:2100){
          day_start <- (y-2095)*365+1

          temp_val <- ncvar_get(nc, "tas",start = c(1,1,day_start),count=c(-1,-1,365))

          temp_val[temp_val < 278.15] <- NA

          output <- lapply(seq_along(lat), function(j) {
            lapply(seq_along(lon), function(i) {
              
              num <- (j-1)*lon_len+i
              
              tas_values <- temp_val[i, j, ]
              new_time <- (365-sum(is.na(tas_values)))*24*3600 #new_time
              mean_above_5 <- mean(tas_values)-278.15 #new_mgdd5
              new_fT <- 2.2^((mean_above_5-25)/10) #new_fT
              
              F_bio <- result[num,3]
              F_msap <- result[num,4]
              F_area <- result[num,5]
              F_rs25a <- result[num,6]
              
              #not consider thermal acclimation
              new_Rs <- F_rs25a*new_fT*F_msap*new_time*F_area*1.2*10^(-15)  #Pg C/year
              #new_Rs_year <- cellStats(new_Rs,stat='sum')
              
              #consider thermal acclimation
              new_rs25a <- exp(log(rs25)-0.1*(mean_above_5-25))
              new_Rsa <- new_rs25a*new_fT*F_msap*new_time*F_area*1.2*10^(-15)  #Pg C/year
              #new_Rsa_year <- cellStats(new_Rsa,stat='sum')
              
              # recudtion of carbon emission
              reduction <- new_Rsa-new_Rs
              
              #fraction
              fraction <- -reduction/new_Rs
              return(list(lon = lon[i]-180,lat = lat[j], mean_above_5 = mean_above_5,time=new_time,new_Rs_year=new_Rs,new_Rsa_year=new_Rsa,
                          reduction=reduction,F_bio=F_bio,F_msap=F_msap,new_fT=new_fT,F_area=F_area,new_time=new_time,fraction=fraction))
            })
          })
          output <- data.frame(do.call(rbind, unlist(output, recursive = FALSE)))
          colnames(output) <- c("lon","lat","mgdd5","time","Rs","Rsa","Reduction","bio","msap","ft","area","time","fraction")
          output_raster <- rasterFromXYZ(output)
          mgdd5_raster <- resample(output_raster$mgdd5,data_bio)
          reduction_raster <- resample(output_raster$Reduction,data_bio)
          #calculate fraction
          Rs_raster <-  resample(output_raster$Rs,data_bio)
          Rsa_raster <-  resample(output_raster$Rsa,data_bio)
          fraction_raster <-  resample(output_raster$fraction,data_bio)
          
          data_clm585[(y-2014),2] <- cellStats(mgdd5_raster,stat='mean')
          data_clm585[(y-2014),3] <- cellStats(reduction_raster,stat='sum')
          data_clm585[(y-2014),4] <- cellStats(Rs_raster,stat='sum')
          data_clm585[(y-2014),5] <- cellStats(Rsa_raster,stat='sum')
          data_clm585[(y-2014),6] <- cellStats(fraction_raster,stat='mean')
          
          pb$tick()
          Sys.sleep(1 / 6)
        }
      }
      plot(data_clm585[,1],data_clm585[,2])#plot mgdd5
      plot(data_clm585[,1],data_clm585[,3])#plot carbon reduction
      
      #Plot single model
      data_clm_tem <- cbind(data_clm126[,1:2],data_clm585[,2])
      data_clm_red <- cbind(data_clm126[,c(1,3)],data_clm585[,3])
      data_clm_Rs <- cbind(data_clm126[,c(1,4)],data_clm585[,4])
      data_clm_Rsa <- cbind(data_clm126[,c(1,5)],data_clm585[,5])
      data_clm_fraction <- cbind(data_clm126[,c(1,6)],data_clm585[,6])
      
      colnames(data_clm_tem) <- c("Year","SSP126","SSP585")
      colnames(data_clm_red) <- c("Year","SSP126","SSP585")
      colnames(data_clm_Rs) <- c("Year","SSP126","SSP585")
      colnames(data_clm_Rsa) <- c("Year","SSP126","SSP585")
      colnames(data_clm_fraction) <- c("Year","SSP126","SSP585")
      
      data_clm_tem <- data.frame(data_clm_tem)
      data_clm_red <- data.frame(data_clm_red)
      data_clm_Rs <- data.frame(data_clm_Rs)
      data_clm_Rsa <- data.frame(data_clm_Rsa)
      data_clm_fraction <- data.frame(data_clm_fraction)
      
      data_clm_m1 <- reshape2::melt(id="Year",data_clm_tem)
      data_clm_m2 <- reshape2::melt(id="Year",data_clm_red)
      data_clm_m3 <- reshape2::melt(id="Year",data_clm_Rs)
      data_clm_m4 <- reshape2::melt(id="Year",data_clm_Rsa)
      data_clm_m5 <- reshape2::melt(id="Year",data_clm_fraction)
      
      data_clm <- cbind(data_clm_m1[,],data_clm_m2[,3],data_clm_m3[,3],data_clm_m4[,3],data_clm_m5[,3])
      colnames(data_clm) <- c("Year","Scenario","T5","reduction","Rs","Rsa","fraction")
      {
        ggplot()+
          geom_point(data = data_cable,aes(x=Year,y=reduction,color=Scenario))+
          geom_line(data = data_cable,aes(x=Year,y=reduction,color=Scenario))+
          
          # Style
          theme_bw()+
          theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                panel.border = element_rect(fill=NA,color="black", linewidth=1.3, linetype="solid"))+
          theme(text = element_text(face = "bold"))+
          scale_color_manual(values =c("#FC4E07","#00AFBB"))+
          
          # Axis
          labs(x=expression(paste("Year")),
               y=expression(paste("Reduce of carbon emmision","  Pg C/year")))+
          theme(axis.title.x = element_text(vjust = 0,size=16),
                axis.title.y = element_text(vjust = 0,size=16))+
          theme(axis.text.x = element_text(size = 12,color="black"),
                axis.text.y = element_text(size = 12,color="black"))
        
      }
      
      
      
      #calculate the ensemble
      #=========================================================
      data_ensemble_tem <- cbind(data_cable[,1:3],data_orchidee[,3],data_jules[,3],data_clm[,3])
      colnames(data_ensemble_tem) <- c("Year","Scenario","CABLE","ORCHIDEE","JULES","CLM5")
      data_ensemble_tem$mean <- rowMeans(data_ensemble_tem[, c(3:6)])
      data_ensemble_tem$max <- apply(data_ensemble_tem[, c(3:6)], 1, max)
      data_ensemble_tem$min <- apply(data_ensemble_tem[, c(3:6)], 1, min)
      
      data_ensemble_red <- cbind(data_cable[,c(1,2,4)],data_orchidee[,4],data_jules[,4],data_clm[,4])
      colnames(data_ensemble_red) <- c("Year","Scenario","CABLE","ORCHIDEE","JULES","CLM5")
      data_ensemble_red$mean <- rowMeans(data_ensemble_red[, c(3:6)])
      data_ensemble_red$max <- apply(data_ensemble_red[, c(3:6)], 1, max)
      data_ensemble_red$min <- apply(data_ensemble_red[, c(3:6)], 1, min)
      
      data_ensemble_fraction <- cbind(data_cable[,1:2],data_cable[,7],data_orchidee[,7],data_jules[,7],data_clm[,7])
      colnames(data_ensemble_fraction) <- c("Year","Scenario","CABLE","ORCHIDEE","JULES","CLM5")
      data_ensemble_fraction$meanfraction <- mean()
      write.csv(data_ensemble_fraction,"G:/proj_StemRespiration/Manuscript/revision2409simulation_future.csv")
      
      #calculate our prediction
      #=========================================================
      #Plot the ensemble
      par(mar = c(6, 6, 6, 6))  # c(bottom, left, top, right)
      ggplot()+
        geom_point(data = data_ensemble_red,aes(x=Year,y=mean,color=Scenario))+
        geom_line(data = data_ensemble_red,aes(x=Year,y=mean,color=Scenario))+
        
        # band
        geom_ribbon(data=data_ensemble_red,aes(ymin = min,ymax=max,x=Year,fill=Scenario), alpha = 0.1) +
        #geom_line(data=data_ensemble_tem,aes(y = max,x=Year,color=Scenario), linetype = "dashed") +
        #geom_line(data=data_ensemble_tem,aes(y = min,x=Year,color=Scenario), linetype = "dashed") +
        
        # Style
        theme_bw()+
        theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
              panel.border = element_rect(fill=NA,color="black", size=1.1, linetype="solid"),
              legend.position = c(0.25,0.25),legend.key.size = unit(1.5, "cm"),
              legend.text = element_text(size = 15),
              legend.title = element_text(size = 18,face="bold"))+
        theme(text = element_text(family="Arial"))+
        scale_color_manual(values =c("#FC4E07","#00AFBB"))+
        
        # Axis
        labs(x=expression(paste("Year")),
             y=expression(paste("Reduction of carbon emmision","  (Pg C/year)")))+
        theme(axis.title.x = element_text(vjust = 0,size=20),
              axis.title.y = element_text(vjust = 0,size=20))+
        theme(axis.text.x = element_text(size = 16,color="black"),
              axis.text.y = element_text(size = 16,color="black"))+
        
        #(b)
        annotate("text", x = 2018, y = -1, label = "B",size=6)
    }
    
    ## CLM5 simulation
    {
      # expression in our model
      rs25 <- exp(intercept-2.5)
      rs25_a <- exp(log(rs25)-0.1*(data_mgdd5-25))
      
      # expression in CLM5
      #CLM5
      rs25_clm5 <- exp(log(rs25)-0.00794*log(10)*(data_mgdd5-25))
      Rs_clm5 <- Rs_ensemble*1.83/10
      
      ######
      #CLM5
      #================================================
      ###SSP126----------------------------------------
      data_clm126_code <- array(dim=c(1,6))
      data_clm126_code[,1] <- "SSP126"
      {
        nc <- nc_open(lsfCLM5[[9,1]])
        lon <- ncvar_get(nc, "lon")
        lon_len <- length(lon)
        lat <- ncvar_get(nc, "lat")
        lat_len <- length(lat)
        
        point_mat <- expand.grid(lon = lon-180, lat = lat)
        coordinates(point_mat) <- ~lon+lat
        proj4string(point_mat) <- proj4string(data_bio)

        F_bio <- raster::extract(data_bio, point_mat) 
        F_fT <- raster::extract(fT, point_mat)
        F_area <- raster::extract(area, point_mat)
        F_msap <- raster::extract(msap,point_mat)
        F_rs25a <- raster::extract(rs25_clm5,point_mat)

        result <- data.frame(point_mat@coords, F_bio = F_bio, F_msap = F_msap, F_area = F_area, F_rs25a=F_rs25a)
        
        #2100
        y <- 2100
        day_start <- (y-2095)*365+1
        temp_val <- ncvar_get(nc, "tas",start = c(1,1,day_start),count=c(-1,-1,365))

        temp_val[temp_val < 278.15] <- NA

        output <- lapply(seq_along(lat), function(j) {
          lapply(seq_along(lon), function(i) {
            
            num <- (j-1)*lon_len+i
            
            tas_values <- temp_val[i, j, ]
            new_time <- (365-sum(is.na(tas_values)))*24*3600 #new_time
            mean_above_5 <- mean(tas_values)-278.15 #new_mgdd5
            new_fT <- 2.2^((mean_above_5-25)/10) #new_fT
            
            F_bio <- result[num,3]
            F_msap <- result[num,4]
            F_area <- result[num,5]
            F_rs25a <- result[num,6]
            
            #not consider thermal acclimation
            new_Rs <- F_rs25a*new_fT*F_msap*new_time*F_area*1.2*10^(-15)  #Pg C/year
            #new_Rs_year <- cellStats(new_Rs,stat='sum')
            
            #consider thermal acclimation
            new_rs25a <- exp(log(rs25)-0.1*(mean_above_5-25))
            new_Rsa <- new_rs25a*new_fT*F_msap*new_time*F_area*1.2*10^(-15)  #Pg C/year
            #new_Rsa_year <- cellStats(new_Rsa,stat='sum')
            
            # recudtion of carbon emission
            reduction <- new_Rsa-new_Rs
            
            #fraction
            fraction <- -reduction/new_Rs
            return(list(lon = lon[i]-180,lat = lat[j], mean_above_5 = mean_above_5,time=new_time,new_Rs_year=new_Rs,new_Rsa_year=new_Rsa,
                        reduction=reduction,F_bio=F_bio,F_msap=F_msap,new_fT=new_fT,F_area=F_area,new_time=new_time,fraction=fraction))
          })
        })
        output <- data.frame(do.call(rbind, unlist(output, recursive = FALSE)))
        colnames(output) <- c("lon","lat","mgdd5","time","Rs","Rsa","Reduction","bio","msap","ft","area","time","fraction")
        output_raster <- rasterFromXYZ(output)
        mgdd5_raster <- resample(output_raster$mgdd5,data_bio)
        reduction_raster <- resample(output_raster$Reduction,data_bio)
        #calculate fraction
        Rs_raster <-  resample(output_raster$Rs,data_bio)
        Rsa_raster <-  resample(output_raster$Rsa,data_bio)
        fraction_raster <-  resample(output_raster$fraction,data_bio)
        
        data_clm126_code[1,2] <- cellStats(mgdd5_raster,stat='mean')
        data_clm126_code[1,3] <- cellStats(reduction_raster,stat='sum')
        data_clm126_code[1,4] <- cellStats(Rs_raster,stat='sum')
        data_clm126_code[1,5] <- cellStats(Rsa_raster,stat='sum')
        data_clm126_code[1,6] <- cellStats(fraction_raster,stat='mean')
        
      }
      
      ###SSP585----------------------------------------
      data_clm585_code <- array(dim=c(1,6))
      data_clm585_code[,1] <- "SSP585"
      {
        #2100
        nc <- nc_open(lsfCLM5[[18,1]])
        lon <- ncvar_get(nc, "lon")
        lon_len <- length(lon)
        lat <- ncvar_get(nc, "lat")
        lat_len <- length(lat)

        point_mat <- expand.grid(lon = lon-180, lat = lat)
        coordinates(point_mat) <- ~lon+lat
        proj4string(point_mat) <- proj4string(data_bio)
        
        F_bio <- raster::extract(data_bio, point_mat) 
        F_fT <- raster::extract(fT, point_mat)
        F_area <- raster::extract(area, point_mat)
        F_msap <- raster::extract(msap,point_mat)
        F_rs25a <- raster::extract(rs25_clm5,point_mat)
        
        result <- data.frame(point_mat@coords, F_bio = F_bio, F_msap = F_msap, F_area = F_area, F_rs25a=F_rs25a)
        y <- 2100
        day_start <- (y-2095)*365+1
        temp_val <- ncvar_get(nc, "tas",start = c(1,1,day_start),count=c(-1,-1,365))
        temp_val[temp_val < 278.15] <- NA

        output <- lapply(seq_along(lat), function(j) {
          lapply(seq_along(lon), function(i) {
            
            num <- (j-1)*lon_len+i
            
            tas_values <- temp_val[i, j, ]
            new_time <- (365-sum(is.na(tas_values)))*24*3600 #new_time
            mean_above_5 <- mean(tas_values)-278.15 #new_mgdd5
            new_fT <- 2.2^((mean_above_5-25)/10) #new_fT
            
            F_bio <- result[num,3]
            F_msap <- result[num,4]
            F_area <- result[num,5]
            F_rs25a <- result[num,6]
            
            #not consider thermal acclimation
            new_Rs <- F_rs25a*new_fT*F_msap*new_time*F_area*1.2*10^(-15)  #Pg C/year
            #new_Rs_year <- cellStats(new_Rs,stat='sum')
            
            #consider thermal acclimation
            new_rs25a <- exp(log(rs25)-0.1*(mean_above_5-25))
            new_Rsa <- new_rs25a*new_fT*F_msap*new_time*F_area*1.2*10^(-15)  #Pg C/year
            #new_Rsa_year <- cellStats(new_Rsa,stat='sum')
            
            # recudtion of carbon emission
            reduction <- new_Rsa-new_Rs
            
            #fraction
            fraction <- -reduction/new_Rs
            return(list(lon = lon[i]-180,lat = lat[j], mean_above_5 = mean_above_5,time=new_time,new_Rs_year=new_Rs,new_Rsa_year=new_Rsa,
                        reduction=reduction,F_bio=F_bio,F_msap=F_msap,new_fT=new_fT,F_area=F_area,new_time=new_time,fraction=fraction))
          })
        })
        output <- data.frame(do.call(rbind, unlist(output, recursive = FALSE)))
        colnames(output) <- c("lon","lat","mgdd5","time","Rs","Rsa","Reduction","bio","msap","ft","area","time","fraction")
        output_raster <- rasterFromXYZ(output)
        mgdd5_raster <- resample(output_raster$mgdd5,data_bio)
        reduction_raster <- resample(output_raster$Reduction,data_bio)
        #calculate fraction
        Rs_raster <-  resample(output_raster$Rs,data_bio)
        Rsa_raster <-  resample(output_raster$Rsa,data_bio)
        fraction_raster <-  resample(output_raster$fraction,data_bio)
        
        data_clm585_code[1,2] <- cellStats(mgdd5_raster,stat='mean')
        data_clm585_code[1,3] <- cellStats(reduction_raster,stat='sum')
        data_clm585_code[1,4] <- cellStats(Rs_raster,stat='sum')
        data_clm585_code[1,5] <- cellStats(Rsa_raster,stat='sum')
        data_clm585_code[1,6] <- cellStats(fraction_raster,stat='mean')
        
      }
      
      #Plot single model
      colnames(data_clm126_code) <- c("SSP","tem","reduction","rs","rsa","fraction")
      colnames(data_clm585_code) <- c("SSP","tem","reduction","rs","rsa","fraction")
      
      # compare clm5 with our model
      
      
    }
  }
}

# figs.S9-11, table S11-15
{
  #@fig.S9(rap)
  {
    n <- 2.5
    a <- 0.365
    b <- 0.8771
    c <- 0.00083
    d <- 0.7478
    s <- 4/7
    data_global$rsmass_rap <- ifelse(data_global$MEASURMENT == "FIELD",  
                                     0.1*data_global$RS_AREA/data_global$WD/(0.01*data_global$RAP)*(data_global$DIAMETER)^(1-2*b)*(s/a*4^b*((n-1)*b+1)*pi^(1-b)),  
                                     data_global$RS_MASS/(0.01*data_global$RAP))
    
    data_global$rsmass25_rap <- data_global$rsmass_rap*2.2^((25-data_global$MTEM)/10)
    data_global$lnrs25 <- log(data_global$rsmass25_rap)
    
    ggplot(data_global,aes(x=T5,y=lnrs25))+
      # Style
      theme_bw()+
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
            panel.border = element_rect(fill=NA,color="black", linewidth =1.3, linetype="solid"))+
      theme(text = element_text(family ="HEL"))+
      # Axis
      labs(x=expression(paste(italic(T[5])," (°C)")),
           y=expression(paste("ln ",italic(r[s25_parenchyma.mass-based]),"  (nmol ",CO[2]," ",g^-1, s^-1,")")))+
      theme(axis.title.x = element_text(vjust = 0, size=17,face="bold"),
            axis.title.y = element_text(vjust = 0,size=17,face="bold"))+
      theme(axis.text.x = element_text(size = 12,color=rgb(0.2,0.2,0.2)),
            axis.text.y = element_text(size = 12,color=rgb(0.2,0.2,0.2)),
            axis.line.x = element_line(linetype=1,color="black",linewidth=1))+
      ylim(-3.5,6.2)+
      xlim(0,25)+
      # Line
      stat_smooth(aes(color=group_ex,fill=group_ex),method = "lm",formula = y~x,se=T,linewidth=1.5,lty=1)+
      scale_fill_manual(values=c("darkgreen","orange","#ecced0"))+
      guides(color="none",fill="none")+
      
      # Ponit
      geom_point(aes(color=group_ex),shape=16,size=3.6,alpha=0.1)+
      scale_color_manual(values= c("darkgreen","orange","#FFB6C1"),
                         name=" ",
                         labels = c("Angiosperms","Gymnosperms","In lab"))+
      # Legend
      theme(legend.position = c(0.2,0.15),
            legend.background = element_rect(fill="transparent"),
            legend.title=element_text(color="black",size=15),
            legend.text = element_text(color="black",size = 14),
            legend.key = element_blank(),
            legend.key.size = unit(1.1, "cm"))+
      guides(
        color = guide_legend(
          direction = "vertical", color = "white",
          override.aes = list(fill = NA)  # 设置置信区间的填充为透明
        )
      )+
      #annotate("text", x = 0.5, y = 6, label = "A",size=7,fontface = "bold")+
      geom_segment(aes(x = 0.5, xend =24,y = 1.95, yend = (-0.1*(24-0.5)+1.95)),
                   color = "#A8021B", linewidth = 1.5,lty=1)
  }
  
  #@fig.S10-11
  {
    #ENVI
    #VPD value_avg (unit:pa)
    vpd_site <- read.csv(".../site_vpd.csv")
    
    data_global_vpd_rs$lnvpd <- log(data_global_vpd_rs$Avg_VPD)
    fit_glm <- lm(LNRS25 ~ T5 +lnvpd, data = data_global_vpd_rs)
    summary(fit_glm)
    
    ###===========Figure===========
    ## s1: visreg vpd and temperature
    {
      visreg(fit_glm, "T5", gg = TRUE,
             points = list(size = 3.2, pch = 16, alpha = 0),  # 设置透明度
             fill = list(fill = "grey", alpha = 0.5),
             line = list(col = "black", size = 1.1, lty = 2)) +
        geom_point(color = "#002FA7", alpha = 0.2, shape = 16, size = 3.6) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(fill = NA, color = "black", size = 1.3, linetype = "solid"),
              text = element_text(family = "A", face = "bold"),
              axis.title.x = element_text(vjust = 0, size = 16),
              axis.title.y = element_text(vjust = 0, size = 16),
              axis.text.x = element_text(size = 12, color = "black"),
              axis.text.y = element_text(size = 12, color = "black")) +
        coord_cartesian(ylim = c(-4.5, 4.5))+# 设置坐标轴范围
        # Style
        theme_bw()+
        theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
              panel.border = element_rect(fill=NA,color="black", size=1.3, linetype="solid"))+
        theme(text = element_text(family ="Arial"))+
        
        # Axis
        labs(x=expression(paste(italic(T[5])," (°C)")),
             y=expression(paste("ln ",italic(r[s25]),"  (nmol ",CO[2]," ",g^-1, s^-1,")")))+
        theme(axis.title.x = element_text(vjust = 0, size=17,face="bold"),
              axis.title.y = element_text(vjust = 0,size=17,face="bold"))+
        theme(axis.text.x = element_text(size = 12,color=rgb(0.2,0.2,0.2)),
              axis.text.y = element_text(size = 12,color=rgb(0.2,0.2,0.2)))+
        geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "black", size = 1.5,lty=2)+
        annotate("text", x = 0.5, y = 4, label = "A",size=6,fontface="bold")
      
      
      visreg(fit_glm, "lnvpd", gg = TRUE,
             points = list(size = 3.2, pch = 16, alpha = 0),  # 设置透明度
             fill = list(fill = "grey", alpha = 0.5),
             line = list(col = "black", size = 1.5, lty = 2)) +
        geom_point(color = "#2F7FC1",alpha=0.2, shape = 16, size = 3.6) +
        
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(fill = NA, color = "black", size = 1.3, linetype = "solid"),
              text = element_text(face = "bold"),
              axis.title.x = element_text(vjust = 0, size = 16),
              axis.title.y = element_text(vjust = 0, size = 16),
              axis.text.x = element_text(size = 12, color = "black"),
              axis.text.y = element_text(size = 12, color = "black"),
              axis.line.x = element_line(linetype = 1, color = "black", size = 1)) +
        coord_cartesian(ylim = c(-4.5, 4.5))+
        # Style
        theme_bw()+
        theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
              panel.border = element_rect(fill=NA,color="black", size=1.3, linetype="solid"))+
        theme(text = element_text(family ="A"))+
        
        # Axis
        labs(x=expression(paste("ln ", italic(VPD), " (Pa)")),
             y=expression(paste(ln," ",italic(r[s25]),"  (nmol ",CO[2]," ",g^-1, s^-1,")")))+
        theme(axis.title.x = element_text(vjust = 0,size=17,family="Arial",face="bold"),
              axis.title.y = element_text(vjust = 0,size=17,family="Arial",face="bold"))+
        theme(axis.text.x = element_text(size = 12,color=rgb(0.2,0.2,0.2),family="Arial"),
              axis.text.y = element_text(size = 12,color=rgb(0.2,0.2,0.2),family="Arial"))+
        geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "black", size = 1.5,lty=2) +
        annotate("text", x = 5, y = 4, label = "B",size=6,fontface="bold")
      
      
    }
    
    ## s2: theoretical prediciton
    #===========================================
    #--------vpd theory examination
    #* rs25 =a*1.6*gs*D
    #* rs25 =a*1.6*(gpp/(ca-ci))*D
    # Pmodel (Stocker et al., 2019; Wang et al., 2017)
    
    tc_c <- 25 #degree C
    atm_c <- 101325 #pa
    fapar_c <- 1 #unitless
    ppfd_c <- 300 #mol m-2 month-1
    vpd <- seq(from = 50, to = 1000, by = 1) #pa
    
    GPP <- unlist(rpmodel( tc_c, vpd, co2=400, fapar_c, ppfd_c, atm_c, elv = NA,returnvar = "gpp", verbose = FALSE ))#unit:g C m-2 month-1
    gs <- unlist(rpmodel( tc_c, vpd, co2=400, fapar_c, ppfd_c, atm_c, elv = NA,returnvar = "gs", verbose = FALSE ))
    
    lnrs25 <- log(1.6*gs*vpd/atm_c)
    lnvpd <- log(vpd)
    plot(lnvpd,lnrs25)
    sensi <- (lnrs25[950]-lnrs25[1])/(lnvpd[950]-lnvpd[1])
    
    data_vpd_predict <- cbind(lnrs25,lnvpd)
    data_vpd_predict <- data.table(data_vpd_predict)
    lm <- lm(lnrs25~lnvpd,data=data_vpd_predict)
    summary(lm)
    
    ggplot(data_vpd_predict,aes(x=lnvpd,y=lnrs25))+
      # Style
      theme_bw()+
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
            panel.border = element_rect(fill=NA,color="black", linewidth =1.3, linetype="solid"))+
      theme(text = element_text(family ="HEL"))+
      # Axis
      labs(x=expression(paste("ln ",italic(VPD)," (Pa)")),
           y=expression(paste("Predicted ln ",italic(r[s25]))))+
      theme(axis.title.x = element_text(vjust = 0, size=17,face="bold"),
            axis.title.y = element_text(vjust = 0,size=17,face="bold"))+
      theme(axis.text.x = element_text(size = 12,color=rgb(0.2,0.2,0.2)),
            axis.text.y = element_text(size = 12,color=rgb(0.2,0.2,0.2)))+
      
      # Ponit
      geom_point(color="black",shape=16,size=3.2,alpha=0.2)+
      guides(color="none",fill="none")+
      # Legend
      theme(legend.position = c(0.15,0.15),
            legend.background = element_rect(fill="transparent"),
            legend.title=element_text(color="black",size=15),
            legend.text = element_text(color="black",size = 14),
            legend.key = element_blank(),
            legend.key.size = unit(1.1, "cm"))+
      guides(
        color = guide_legend(
          direction = "vertical", color = "white",
          override.aes = list(fill = NA)  # 设置置信区间的填充为透明
        )
      )+
      annotate("text", x = 0.5, y = 4.8, label = "A",size=7,fontface = "bold")
  }
  
  #@table S10
  {
    rsdata_new <- data_global
    rsdata_new$lnrs25_q1.6<- log(rsdata_new$rsmass_sap*1.6^((25-rsdata_new$MTEM)/10))
    rsdata_new$lnrs25_q1.8<- log(rsdata_new$rsmass_sap*1.8^((25-rsdata_new$MTEM)/10))
    rsdata_new$lnrs25_q2.0<- log(rsdata_new$rsmass_sap*2^((25-rsdata_new$MTEM)/10))
    rsdata_new$lnrs25_q2.2<- log(rsdata_new$rsmass_sap*2.2^((25-rsdata_new$MTEM)/10))
    rsdata_new$lnrs25_q2.4<- log(rsdata_new$rsmass_sap*2.4^((25-rsdata_new$MTEM)/10))
    rsdata_new$lnrs25_q2.6<- log(rsdata_new$rsmass_sap*2.6^((25-rsdata_new$MTEM)/10))
    
    
    eiv.rs <- eivreg(lnrs25_q1.6 ~ T5 , data=rsdata_new,reliability = reliability)
    summary(eiv.rs)
    eiv.rs <- eivreg(lnrs25_q1.8 ~ T5 , data=rsdata_new,reliability = reliability)
    summary(eiv.rs)
    eiv.rs <- eivreg(lnrs25_q2.0 ~ T5 , data=rsdata_new,reliability = reliability)
    summary(eiv.rs)
    eiv.rs <- eivreg(lnrs25_q2.2 ~ T5 , data=rsdata_new,reliability = reliability)
    summary(eiv.rs)
    eiv.rs <- eivreg(lnrs25_q2.4 ~ T5 , data=rsdata_new,reliability = reliability)
    summary(eiv.rs)
    eiv.rs <- eivreg(lnrs25_q2.6 ~ T5 , data=rsdata_new,reliability = reliability)
    summary(eiv.rs)
    
    ggplot(rsdata_new, aes(x = mgdd5, y = lnrs25_q1.6,color=group_ex)) +
      geom_point(size = 2, alpha=0.2) +  
      # Line
      stat_smooth(aes(color=group_ex,fill=group_ex),method = "lm",formula = y~x,se=T,lty=2,linewidth=1.2)+
      labs(title = "Q10 = 1.6",
           x = "mgdd5",
           y = "lnrs25")
  }
  
  #@table S12.13
  {
    colnames(data_global)
    eiv.rs <- eivreg(DIAMETER ~ T5 , data=data_global,reliability = reliability)
    summary(eiv.rs)
    
    eiv.rs <- eivreg(LNRS25 ~ T5+DIAMETER , data=data_global,reliability = reliability)
    summary(eiv.rs)
  }
  
  #@table S14
  {
    insitu <- data.table(data_global[which((data_global$MEASURMENT=="FIELD")),])
    exsitu <- data.table(data_global[which((data_global$MEASURMENT=="LAB")),])
    angio <- data.table(insitu[which(insitu$GROUP=="Angiosperms"),])
    conifer <- data.table(insitu[which(insitu$GROUP=="Gymnosperms"),])

    eiv.rs <- eivreg(LNRS25~T5,data=exsitu,reliability = reliability)
    summary(eiv.rs)
    eiv.rs <- eivreg(LNRS25 ~ T5 , data=angio,reliability = reliability)
    summary(eiv.rs)
    eiv.rs <- eivreg(LNRS25 ~ T5 , data=conifer,reliability = reliability)
    summary(eiv.rs)
  }
  
}
