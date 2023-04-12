#Clear 
rm(list=ls())
#Set Test Run
test.run=FALSE
#Libraries
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(ggcorrplot)
library(vegan)
library(Hmsc)
library(snow)
library(corrplot)
#Data
p1<-read.csv("allsites_algae.csv")
p2<-read.csv("Traits_algae.csv")
#Set Seed
set.seed(123)

#Environmental Features
cov=p1[,c(5,6)]
#can add and manipulate variables to cov
#cov and X are the same for now
X=as.data.frame(cov[c("Temperature","day")])
XFormula= ~Temperature+day

#Traits
TrData=select(p2, 3:5, 12:16, 18:21)
rownames(TrData)<-p2[,1]
#colnames(TrData)<-"algae_or_animal"+"Algae Class"
TrFormula=~Rhodophyta+Chlorophyta+Phaeophyta+Erect+Prostrate+Small+Medium+Large+Perennial+Pseudo_perennial+Aseasonal_annual+Seasonal_annual

#Species Occurrence Box
Y = as.matrix(p1[,7:38])
Y<-Y/100
#Study Design
# Variable selection 
studyDesign=p1[1]
#Sites and Seasons need to be converted to numbers
#CapDan=1,DallasRock=2,DanSpot=3, mile23=4, reef200200=5)
#Fall 2011=1, Fall 2012=2, 
#Spring 2012=3, Summer 2011=4,
#Summer 2012=5, Winter 2011-12=6,
#Winter 2012-13=7, 
studyDesign$Site=as.numeric(as.factor(studyDesign$Site))
studyDesign=data.frame(apply(studyDesign,2,as.factor))
studyDesign$Site=as.factor(studyDesign$Site)
# Variable structuring 
Site=HmscRandomLevel(units=studyDesign$Site)
ranlevels=list(Site=Site) 
ranlevels

#Creating Model
simul <- Hmsc(Y=Y, XData = X,
              XFormula = XFormula,
              TrData = TrData,
              TrFormula = TrFormula,
              studyDesign = studyDesign,
              ranLevels  = ranlevels,
              distr = "probit")  

#Running Model
nChains = 2
if(test.run){
  thin=10
  samples=50
  transient=ceiling(thin*samples*.5)
}else{
  thin=10
  samples=1000
  transient=ceiling(thin*samples*.5)
}
m=sampleMcmc(simul,
             samples = samples,
             thin = thin,
             transient = transient,
             nChains = nChains, 
             nParallel = nChains)
mpost <- convertToCodaObject(m)
preds = computePredictedValues(m)
MF = evaluateModelFit(hM=m, predY=preds)
VP <- computeVariancePartitioning(m)

## Convergence tests
par(mar=rep(2,4))
#Visual chain tests for different coefficients of interest 
plot(mpost$Beta[,1:2])
plot(mpost$Gamma[,1:2])

#Organization for Diagnostics
psrf.V = gelman.diag(mpost$Beta[,1:32])$psrf%>%
  as_tibble() %>% dplyr::rename(psrf_v = `Point est.`)
ess.beta <- effectiveSize(mpost$Beta) %>%
  as_tibble() %>% dplyr::rename(ess_beta = value)
ess.v <- effectiveSize(mpost$V)%>%
  as_tibble() %>% dplyr::rename(ess_v = value)
psrf.beta <- gelman.diag(mpost$Beta, multivariate=FALSE)$psrf%>%
  as_tibble() %>% dplyr::rename(psrf_beta = `Point est.`)
#Effective Sample Size
ggplot(ess.beta, aes(x=ess_beta))+
  geom_histogram()+
  xlab("Effective Sample Size")
if(!test.run){
  ggsave("ess_algae.pdf", width=11.5, height=8.5)
}

#Gelman Diagnostic                      
ggplot(psrf.beta, aes(x=psrf_beta))+
  geom_histogram()+
  xlab("Gelman Diagnostic")
if(!test.run){
  ggsave("gelman_algae.pdf", width=11.5, height=8.5)
}

#Variance Partitioning
vp_df <- VP$vals%>%
  as_tibble(rownames = "variable") %>%
  pivot_longer(cols=names(.)[2:ncol(.)],
               names_to = "Species",
               values_to = "value") %>%
  na.omit()
vp_summary <- vp_df %>%
  group_by(variable) %>%
  summarise(value = mean(value)) %>%
  ungroup()
vp_order <- vp_df %>%
  filter(variable == "Temperature") %>%
  arrange(value)%>%
  mutate(Species_f = factor(Species, levels = .$Species)) %>%
  dplyr::select(Species, Species_f)
#
vp <- left_join(vp_df, vp_order)
vp %>%
  ggplot(aes(x=value,y=Species_f, fill = variable)) +
  geom_bar(stat="identity")+
  theme_classic() +
  ylab("Species") +
  xlab("Proportion of Variance Explained") +
  scale_fill_brewer(palette = "Dark2")+
  theme(legend.position = c(1,.315),
        legend.justification = c(1,0),
        legend.background = element_rect(color="black")) +
  ggtitle("Variance Partitioning, Occurrence Model")
if(!test.run){
  ggsave("vp_algae.pdf", width=12, height=7)
}

#Temperature Correlation Heatmap
postBeta <- getPostEstimate(m, parName = "Beta")
means <- postBeta$mean %>%
  as_tibble() %>%
  rowid_to_column("env_var") %>%
  mutate(env_var = c("intercept",VP$groupnames)) %>%
  pivot_longer(cols=names(.)[2:ncol(.)], names_to = "Species", values_to = "Mean")
postBeta$mean %>%
  as_tibble() %>%
  rowid_to_column("env_var")
postBeta$support %>%
  as_tibble() %>%
  rowid_to_column("env_var") %>%
  mutate(env_var = c("intercept",VP$groupnames))
VP$groupnames
c("intercept",VP$groupnames)
postBeta$support %>%
  as_tibble() %>%
  rowid_to_column("env_var")
c("intercept",VP$groupnames)
supported <- postBeta$support %>%
  as_tibble() %>%
  rowid_to_column("env_var") %>%
  mutate(env_var = c("intercept",  "Temperature", "day")) %>%
  pivot_longer(cols=names(.)[2:ncol(.)], names_to = "Species", values_to = "Support") %>%
  #filter(Support >0.90|Support<0.10, env_var != "intercept")%>%
  left_join(means, by = c("env_var", "Species"))%>%
  mutate(sign = ifelse(Mean>0, "+", "-"))
postBeta$support %>%
  as_tibble() %>%
  rowid_to_column("env_var") %>%
  mutate(env_var = c("intercept",  "Temperature", "day")) %>%
  pivot_longer(cols=names(.)[2:ncol(.)],names_to = "Species",values_to = "Support")
means <- postBeta$mean %>%
  as_tibble() %>%
  rowid_to_column("env_var") %>%
  mutate(env_var = c("intercept",VP$groupnames)) %>%
  pivot_longer(cols=names(.)[2:ncol(.)], names_to = "Species", values_to = "Mean")
means <- postBeta$mean %>%
  as_tibble() %>%
  rowid_to_column("env_var") %>%
  mutate(env_var =c("intercept",  "Temperature", "day")) %>%
  pivot_longer(cols=names(.)[2:ncol(.)], names_to = "Species", values_to = "Mean")
supported <- postBeta$support %>%
  as_tibble() %>%
  rowid_to_column("env_var") %>%
  mutate(env_var = c("intercept",  "Temperature", "day")) %>%
  pivot_longer(cols=names(.)[2:ncol(.)],names_to = "Species",values_to = "Support") %>%
  filter(Support >0.95|Support<0.05,env_var != "intercept") %>%
  left_join(means, by = c("env_var", "Species"))%>%
  mutate(sign = ifelse(Mean>0, "+", "-"))
#supported <- postBeta$support %>%
#as_tibble() %>%
#rowid_to_column("env_var") %>%
#mutate(env_var = c("intercept",  "Temperature")) %>%
#pivot_longer(cols=names(.)[2:ncol(.)], names_to = "Species", values_to = "Support") %>%
#filter(Support >0.90|Support<0.10, env_var != "intercept") %>%
#left_join(means, by = c("env_var", "Species"))%>%
#mutate(sign = ifelse(Mean>0, "+", "-")) %>%
ggplot(supported,mapping=aes(x=env_var,y=Species, fill = Mean, color = sign)) +
  geom_tile(lwd=.5) +
  theme_clean()+
  scale_fill_gradient2(mid = "grey90") +
  scale_color_manual(values = c(("red"), ("blue"))) +
  guides(color = "none")+
  theme(axis.text.x = element_text(angle=45, vjust=1,hjust = 1),
        axis.title = element_blank())
if(!test.run){
  ggsave("heatmap_algae.pdf", width=12, height=7)
}

#Species Association
OmegaCor = computeAssociations(m) 
supportLevel = 0.95 
toPlot = ((OmegaCor[[1]]$support>supportLevel)
          + (OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[1]]$mean

toPlot_p = ((OmegaCor[[1]]$support>supportLevel)
            + (OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[1]]$mean

hmdf_mean <- OmegaCor[[1]]$mean %>%
  as.matrix
hmdf_support <- OmegaCor[[1]]$support %>%
  as.matrix
omega_plot<- ggcorrplot::ggcorrplot(hmdf_mean, type="lower",
                                    hc.order=TRUE, title="Occurence")
omega_plot
if(!test.run){
  ggsave(omega_plot, filename="speciesassociation_algae.pdf", width=18, height=18)
}

#Residual Associations Between Species
OmegaCor = computeAssociations(m) 
supportLevel = 0.95 
toPlot = ((OmegaCor[[1]]$support>supportLevel)
          + (OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[1]]$mean 
corrplot(toPlot, method = "color", col=colorRampPalette(c("blue","white","red"))(200),
         tl.cex=.3, tl.col="black",
         title=paste("random effect level:", m$rLNames[1]), mar=c(0,0,1,0))

#Temperature-Species Richness Gradient
GradientTemp = constructGradient(m,
                                 focalVariable = "Temperature")
predYTemp = predict(m,
                    XData = GradientTemp$XDataNew, 
                    studyDesign = GradientTemp$studyDesignNew,
                    ranLevels = GradientTemp$rLNew,
                    expected = TRUE)

plotGradient(m,
             GradientTemp,
             pred=predYTemp,
             measure="S",
             showData = TRUE)

#Environmental Effects on Species
postBeta_b = getPostEstimate(m, parName = "Beta")

par(mar=c(5,11,2.5,0))

plotBeta(m,
         post = postBeta_b, 
         plotTree = F,
         spNamesNumbers = c(T,F))
plotBeta(m, 
         post = postBeta_b,
         param = "Mean",
         plotTree = F,  
         spNamesNumbers = c(T,F))
#Environmental Effect on Traits
par(mar=c(5,11,2.5,0))
postGamma = getPostEstimate(m, parName = "Gamma")
plotGamma(m, post = postGamma)
#No support at 0.95
par(mar=c(5,11,2.5,0))
postGamma = getPostEstimate(m, parName = "Gamma")
plotGamma(m, post = postGamma, supportLevel = 0.1)

#Temperature-Traits Gradient
plotGradient(m, 
             GradientTemp,
             pred=predYTemp,
             measure="T",
             index = 2,
             showData = TRUE)

plotGradient(m, 
             GradientTemp,
             pred=predYTemp,
             measure="T",
             index = 3,
             showData = TRUE)

plotGradient(m, 
             GradientTemp,
             pred=predYTemp,
             measure="T",
             index = 4,
             showData = TRUE)
plotGradient(m, 
             GradientTemp,
             pred=predYTemp,
             measure="T",
             index = 5,
             showData = TRUE)
plotGradient(m, 
             GradientTemp,
             pred=predYTemp,
             measure="T",
             index = 6,
             showData = TRUE)
plotGradient(m, 
             GradientTemp,
             pred=predYTemp,
             measure="T",
             index = 7,
             showData = TRUE)
plotGradient(m, 
             GradientTemp,
             pred=predYTemp,
             measure="T",
             index = 8,
             showData = TRUE)
plotGradient(m, 
             GradientTemp,
             pred=predYTemp,
             measure="T",
             index = 9,
             showData = TRUE)
plotGradient(m, 
             GradientTemp,
             pred=predYTemp,
             measure="T",
             index = 10,
             showData = TRUE)
plotGradient(m, 
             GradientTemp,
             pred=predYTemp,
             measure="T",
             index = 11,
             showData = TRUE)
plotGradient(m, 
             GradientTemp,
             pred=predYTemp,
             measure="T",
             index = 12,
             showData = TRUE)
plotGradient(m, 
             GradientTemp,
             pred=predYTemp,
             measure="T",
             index = 13,
             showData = TRUE)

#Day-Species Richness Gradient
GradientDay = constructGradient(m,
                                focalVariable = "day")
predYDay = predict(m,
                   XData = GradientDay$XDataNew, 
                   studyDesign = GradientDay$studyDesignNew,
                   ranLevels = GradientDay$rLNew,
                   expected = TRUE)
plotGradient(m,
             GradientDay,
             pred=predYDay,
             measure="S",
             showData = TRUE)
#Day-Traits Gradient
plotGradient(m, 
             GradientDay,
             pred=predYDay,
             measure="T",
             index = 2,
             showData = TRUE)

plotGradient(m, 
             GradientDay,
             pred=predYDay,
             measure="T",
             index = 3,
             showData = TRUE)

plotGradient(m, 
             GradientDay,
             pred=predYDay,
             measure="T",
             index = 4,
             showData = TRUE)
plotGradient(m, 
             GradientDay,
             pred=predYDay,
             measure="T",
             index = 5,
             showData = TRUE)
plotGradient(m, 
             GradientDay,
             pred=predYDay,
             measure="T",
             index = 6,
             showData = TRUE)
plotGradient(m, 
             GradientDay,
             pred=predYDay,
             measure="T",
             index = 7,
             showData = TRUE)
plotGradient(m, 
             GradientDay,
             pred=predYDay,
             measure="T",
             index = 8,
             showData = TRUE)
plotGradient(m, 
             GradientDay,
             pred=predYDay,
             measure="T",
             index = 9,
             showData = TRUE)
plotGradient(m, 
             GradientDay,
             pred=predYDay,
             measure="T",
             index = 10,
             showData = TRUE)
plotGradient(m, 
             GradientDay,
             pred=predYDay,
             measure="T",
             index = 11,
             showData = TRUE)
plotGradient(m, 
             GradientDay,
             pred=predYDay,
             measure="T",
             index = 12,
             showData = TRUE)
plotGradient(m, 
             GradientDay,
             pred=predYDay,
             measure="T",
             index = 13,
             showData = TRUE)
#New Multispecies GGPLOT Gradients

n_runs <- nChains*samples
pred_df <- do.call("rbind", predYTemp) %>%
  as_tibble() %>%
  mutate(Temperature = rep(GradientTemp$XDataNew$Temperature,
                      n_runs),
         run = rep(1:n_runs,each=20)) %>%
  pivot_longer(values_to = "cover", names_to = "Species", -c(Temperature,run)) %>%
  #left_join(prevalence) %>%
  #filter(prevalence > 1) %>%
  #arrange(desc(prevalence)) %>%
  mutate(Species_f = factor(Species, levels = unique(.$Species)))
temp_species <- pred_df %>%
  ggplot(aes(x=Temperature, y=cover)) +
  geom_line(alpha = 0.03, aes(group=run), key_glyph="rect")+
  # geom_line(data = pred_df_grazing,
  #           lwd=1, alpha=0.95, color="black", aes(y=mean))+
  facet_wrap(~Species_f, nrow=4)+
  xlab("Temperature") +
  ylab("Probability of Occurrence") +
  guides(color=guide_legend(override.aes = list(alpha=1)))+
  # scale_color_brewer(palette = "Dark2") +
  # scale_color_manual(values = pal_nat)+
  theme(legend.position = "right",
        #strip.text = element_markdown(),
        panel.border = element_rect(fill=NA, size=0.75),
        legend.justification = c(1,0),
        legend.title = element_blank())
temp_species
if(!test.run){
  ggsave(temp_species, filename="tempspecies.pdf", width=12, height=8)
}
temp_gradient <- pred_df %>%
  ggplot(aes(x=Temperature, y=cover)) +
  geom_line(alpha = 0.03, aes(group=run), key_glyph="rect")+
  # geom_line(data = pred_df_grazing,
  #           lwd=1, alpha=0.95, color="black", aes(y=mean))+
  #facet_wrap(~Species_f, nrow=4)+
  xlab("Temperature") +
  ylab("Probability of Occurrence") +
  guides(color=guide_legend(override.aes = list(alpha=1)))+
  # scale_color_brewer(palette = "Dark2") +
  # scale_color_manual(values = pal_nat)+
  theme(legend.position = "right",
        #strip.text = element_markdown(),
        panel.border = element_rect(fill=NA, size=0.75),
        legend.justification = c(1,0),
        legend.title = element_blank())
temp_gradient

n_runs <- nChains*samples
pred_df <- do.call("rbind", predYDay) %>%
  as_tibble() %>%
  mutate(Day = rep(GradientDay$XDataNew$day,
                           n_runs),
         run = rep(1:n_runs,each=20)) %>%
  pivot_longer(values_to = "cover", names_to = "Species", -c(day,run)) %>%
  #left_join(prevalence) %>%
  #filter(prevalence > 1) %>%
  #arrange(desc(prevalence)) %>%
  mutate(Species_f = factor(Species, levels = unique(.$Species)))
day_species <- pred_df %>%
  ggplot(aes(x=day, y=cover)) +
  geom_line(alpha = 0.03, aes(group=run), key_glyph="rect")+
  # geom_line(data = pred_df_grazing,
  #           lwd=1, alpha=0.95, color="black", aes(y=mean))+
  facet_wrap(~Species_f, nrow=4)+
  xlab("Day") +
  ylab("Probability of Occurrence") +
  guides(color=guide_legend(override.aes = list(alpha=1)))+
  # scale_color_brewer(palette = "Dark2") +
  # scale_color_manual(values = pal_nat)+
  theme(legend.position = "right",
        #strip.text = element_markdown(),
        panel.border = element_rect(fill=NA, size=0.75),
        legend.justification = c(1,0),
        legend.title = element_blank())
day_species
if(!test.run){
  ggsave(temp_species, filename="dayspecies.pdf", width=12, height=8)
}
day_gradient <- pred_df %>%
  ggplot(aes(x=day, y=cover)) +
  geom_line(alpha = 0.03, aes(group=run), key_glyph="rect")+
  # geom_line(data = pred_df_grazing,
  #           lwd=1, alpha=0.95, color="black", aes(y=mean))+
  #facet_wrap(~Species_f, nrow=4)+
  xlab("Day") +
  ylab("Probability of Occurrence") +
  guides(color=guide_legend(override.aes = list(alpha=1)))+
  # scale_color_brewer(palette = "Dark2") +
  # scale_color_manual(values = pal_nat)+
  theme(legend.position = "right",
        #strip.text = element_markdown(),
        panel.border = element_rect(fill=NA, size=0.75),
        legend.justification = c(1,0),
        legend.title = element_blank())
day_gradient
