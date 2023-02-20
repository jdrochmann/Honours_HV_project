## load the necessary libraries
library(tidyverse)
library(lme4)
library(sjPlot)
library(lmerTest)
library(webshot) #webshot::install_phantomjs()
library(emmeans)
library(factoextra)
library(ggpubr)
library(performance)
#install.packages("equatiomatic") # converts a model to LaTeX
library(equatiomatic)

## read in the cleaned dataset
source("R/HV_data_cleanup.R")
source("R/functions.R")

# create the dataframe of species names sorted alphabetically
species_names<-sort(unique(survival_data$species))
scientific_names<-sort(scientific_names$scientific_name)
scientific_prov_names <- c("Eucalyptus crebra","Eucalyptus tereticornis")
short_names <- c("A. glaucocarpa", "A. salicina", "A. floribunda/subvelutina",
                 "C. citriodora subsp. variegata", "C. intermedia", 
                 "C. tesselaris", "E. crebra", "E. crebra (dry)", "E. moluccana", 
                 "E. tereticornis", "E. tereticornis (dry)", "L. suaveolens",
                 "M. irbyana")
trait_names <- c("LA","SLA","LDMC","SSD","RTD","log(SRL)","log(Seed Mass)")

# Colour palettes  to use for plottings 
# aiming for pretty colour groups based on eucalypts, corymbias, acacias, 
# angophoras and melaleuca
line_pal <- c('#DC050C','#1965B0','#F1932D')

spec_pal_rainbow <- c('#4EB265','#90C987',
                      '#882E72',
                      '#1965B0','#5289C7','#7BAFDE',
                      '#F7F056','#F6C141','#F1932D','#E8601C','#DC050C',
                      '#D1BBD7','#AE76A3',
                      '#777777')


######## ######## ######## ######## ######## ######## ######## ######## ######## 
######## SURVIVAL MODELS ######## 
######## ######## ######## ######## ######## ######## ######## ######## ######## 
# Preliminary simple survival model to see if quadratic effects needed
surv_mod_htonly <- glmer(alive_resurvey ~ sqrt(height_baseline_cm) + 
                           height_baseline_cm + 
                           (sqrt(height_baseline_cm)+height_baseline_cm|species) + 
                           (1|block/plot),
                         data=survival_data, family=binomial(link="logit"))
summary(surv_mod_htonly)
r2(surv_mod_htonly)
## model for survival using the trait PCA
surv_mod_pc <- glmer(alive_resurvey ~ sqrt(height_baseline_cm) + innoculation + 
                       compost + mean_elev + cross_fall + trait_pc1 + trait_pc2 + 
                       trait_pc3 + date_diff_b2r + innoculation:compost + 
                       trait_pc1:sqrt(height_baseline_cm) + 
                       trait_pc2:sqrt(height_baseline_cm) + 
                       trait_pc3:sqrt(height_baseline_cm) +
                       trait_pc1:innoculation + trait_pc2:innoculation + 
                       trait_pc3:innoculation + trait_pc1:compost + 
                       trait_pc2:compost + trait_pc3:compost + (1|plot) + 
                       (sqrt(height_baseline_cm) + innoculation + compost|species), 
                     data = survival_data, family = binomial(link = "logit"),
                     control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(surv_mod_pc) 
r2(surv_mod_pc)

# A more simple survival model excluding non-significant interactions in order 
# to better understand the main effects. 
# ht:trait - traits modulate species survival relationships w size but we only 
# retain sig ht:trait interactions
surv_mod_simp <- glmer(alive_resurvey ~ sqrt(height_baseline_cm) + innoculation +
                             compost + mean_elev + cross_fall + 
                             trait_pc1 + trait_pc2 + 
                             trait_pc3 + date_diff_b2r + innoculation:compost +
                         sqrt(height_baseline_cm):trait_pc3 + (1|plot) + 
                             (sqrt(height_baseline_cm) + innoculation + compost|species), 
                           data = survival_data, family = binomial(link = "logit"), 
                           control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(surv_mod_simp)
r2(surv_mod_simp)


######## ######## ######## ######## ######## ######## ######## ######## ######## 
######## HEIGHT GROWTH MODELS ######## 
######## ######## ######## ######## ######## ######## ######## ######## ########
# Preliminary simple height-growth model to see if quadratic effects needed
growth_mod_htonly <- lmer(sqrt(growth_mm_day) ~ sqrt(height_baseline_cm) + 
                            height_baseline_cm + 
                            (sqrt(height_baseline_cm)+height_baseline_cm|species) + 
                            (1|block/plot), 
                          pos_growth_data)
summary(growth_mod_htonly)

growth_mod_pc <- lmer(sqrt(growth_mm_day) ~ sqrt(height_baseline_cm) + 
                        innoculation + compost + mean_elev + cross_fall + 
                        trait_pc1 + trait_pc2 + trait_pc3 +  
                        date_diff_b2r + innoculation:compost + 
                        trait_pc1:sqrt(height_baseline_cm) + 
                        trait_pc2:sqrt(height_baseline_cm) + 
                        trait_pc3:sqrt(height_baseline_cm) + 
                        trait_pc1:innoculation + trait_pc2:innoculation + 
                        trait_pc3:innoculation + trait_pc1:compost + 
                        trait_pc2:compost + trait_pc3:compost +(1|block/plot) + 
                        (sqrt(height_baseline_cm) + compost|species), 
                      pos_growth_data, 
                      control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(growth_mod_pc) 
r2(growth_mod_pc)

growth_mod_simp <- lmer(sqrt(growth_mm_day) ~ sqrt(height_baseline_cm)  + 
                          innoculation + compost + mean_elev + cross_fall + 
                          trait_pc1 + trait_pc2 + trait_pc3 + 
                          date_diff_b2r + innoculation:compost + 
                          (1|block/plot) + 
                          (sqrt(height_baseline_cm) + compost|species), 
                        pos_growth_data, 
                        control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(growth_mod_simp) 
r2(growth_mod_simp)
######## ######## ######## ######## ######## ######## ######## ######## ######## 
######## OUTPUT MODELS TO PRETTY TABLE ######## 
######## ######## ######## ######## ######## ######## ######## ######## ######## 
# Output is html - could save as file = "summary.doc" for word doc
htonly_model_table <- tab_model(surv_mod_htonly, growth_mod_htonly, show.re.var=TRUE, 
                                dv.labels= c("Survival Model", "Height Growth Model"),
                                pred.labels = c("(Intercept)", 
                                                "Initial height [sqrt]",
                                                "Initial height"),
                                string.pred = "Coeffcient",
                                string.ci = "CI (95%)",
                                string.p = "P-Value",
                                transform = NULL,
                                file = "Outputs/Tables/htonly_model_table.html")

htonly_model_table
# then take this html file and make .png file
webshot("Outputs/Tables/htonly_model_table.html", "Outputs/Tables/htonly_model_table.pdf") 

pc_model_table <- tab_model(surv_mod_pc, growth_mod_pc, show.re.var= TRUE, 
                            dv.labels= c("Survival Model", "Height Growth Model"),
                            pred.labels =c("(Intercept)", 
                                           "Initial height [sqrt]",
                                           "Innoculation treatment [YES]", 
                                           "Compost treatment [YES]",
                                           "Mean elevation", "Cross-fall",
                                           "Trait PC1", "Trait PC2", 
                                           "Trait PC3", "Census duration",
                                           "Innoculation treatment [YES] : 
                                           compost treatment [YES]", "Initial 
                                           height [sqrt] : trait PC1", "Initial 
                                           height [sqrt] : trait PC2", "Initial 
                                           height [sqrt] : trait PC3", 
                                           "Innoculation treatment [YES] : trait PC1", 
                                           "Innoculation treatment [YES] : trait PC2", 
                                           "Innoculation treatment [YES] : trait PC3",
                                           "Compost treatment [YES] : trait PC1", 
                                           "Compost treatment [YES] : trait PC2", 
                                           "Compost treatment [YES] : trait PC3"),
                            string.pred = "Coeffcient",
                            string.ci = "CI (95%)",
                            string.p = "P-Value",
                            transform = NULL,
                            file = "Outputs/Tables/pc_model_table.html")

pc_model_table
# then take this html file and make .png file
webshot("Outputs/Tables/pc_model_table.html", "Outputs/Tables/pc_model_table.pdf") 

simp_model_table <- tab_model(surv_mod_simp, growth_mod_simp, show.re.var= TRUE, 
                              dv.labels= c("Survival Model", "Height Growth Model"), 
                              pred.labels =c("(Intercept)", 
                                             "Initial height [sqrt]",
                                             "Innoculation treatment [YES]", 
                                             "Compost treatment [YES]",
                                             "Mean elevation", "Cross-fall",
                                             "Trait PC1", "Trait PC2", 
                                             "Trait PC3", "Census duration",
                                             "Innoculation treatment [YES] : 
                                           compost treatment [YES]", "Initial 
                                           height [sqrt] : trait PC3"),
                              string.pred = "Coeffcient",
                              string.ci = "CI (95%)",
                              string.p = "P-Value",
                              transform = NULL,
                              file = "Outputs/Tables/simp_model_table.html")

simp_model_table
# then take this html file and make .png file
webshot("Outputs/Tables/simp_model_table.html", "Outputs/Tables/simp_model_table.pdf") 


######## ######## ######## ######## ######## ######## ######## ######## ######## 
######## Trait PCA ######## 
######## ######## ######## ######## ######## ######## ######## ######## ######## 
#### PCA based on measurements of functional traits for four replicates of species;...
avg_trait_vals <- as.data.frame(avg_trait_vals)
row.names(avg_trait_vals)<-scientific_names
avg_trait_vals_w_colnames <- avg_trait_vals
colnames(avg_trait_vals_w_colnames)<- c("Species", "LA", "SLA", "LDMC", "SSD", "RTD", 
                             "log(SRL)", "log(seed_mass)", "PCA axis 1", 
                             "PCA axis 2", "PCA axis 3")
#row.names(avg_trait_vals)<-short_names
trait_pca<-princomp(avg_trait_vals_w_colnames[,c(2:8)], cor = TRUE)

traitpc <- fviz_pca_biplot(trait_pca, col.var=spec_pal_rainbow[5], 
                           labelsize=5, arrowsize = 1,
                           pointsize=3, repel = TRUE, title="")
set.seed(5)
traitpc + theme_bw() + theme(plot.background = element_blank(), 
                             panel.grid.major = element_blank(), 
                             panel.grid.minor = element_blank(),
                             axis.text = element_text(size = 15),
                             text = element_text(size = 15)) +
  labs(y = "PC2 (29%)", x = "PC1 (34.5%)") + xlim(-3.7, 3.7) + ylim (-3.7, 3.7)
ggsave("Outputs/Figures/figure_pca.pdf", width = 10, height = 10) 

summary(trait_pca) 
loadings(trait_pca) 
write.csv(unclass(loadings(trait_pca)), file = "Outputs/pca.loadings.csv")

######## ######## ######## ######## ######## ######## ######## ######## ######## 
######## FIGURE: Species survival, compost ON/OFF ######## 
######## ######## ######## ######## ######## ######## ######## ######## ######## 

species_surv_coefs <- coef(surv_mod_simp)$species
species_surv_coefs[1,] ##this is to know what to put for the cbind. 

## values for plotting - holding everything constant at the mean except for compost ON/OFF
mean_mean_elev <- mean(survival_data$mean_elev, na.rm=T) #88.91011
mean_cross_fall <- mean(survival_data$cross_fall, na.rm=T) #1.002128
mean_date_diff_b2r <- mean(survival_data$date_diff_b2r, na.rm=T) #257.0161

## plot x = sqrt(height_baseline_cm) with one curve for compost = NO and one curve for compost = Y
#Need cbind(1,sqrt(height),innoculation OFF, compost ON/OFF, mean_mean_elev, 
# mean_cross_fall, trait_pc1, trait_pc2, trait_pc3, 0.5 (avg of block B ON/OFF)?-REMOVED, 
# mean_date_diff_b2r, 0*0/1, sqrt(height)*mean_pc3)
pdf(file="Outputs/Figures/figure_2.pdf", width=12, height=12)
par(mfrow=c(4,4))
for(i in 1:length(species_names)){
  species_trait_data<-droplevels(subset(survival_data, species==species_names[i]))
  with(species_trait_data, plot(jitter(alive_resurvey, amount=0.05) ~ sqrt(height_baseline_cm),
                                col=alpha(ifelse(compost=="NO", line_pal[2], line_pal[1]),0.5), pch=16, 
                                xaxt="n", cex.lab=1.15, xlim=c(0,9),
                                xlab = expression(Initial~Height~(cm)), 
                                ylab=expression(Probability~of~Survival), main=scientific_names[i], cex.main=1.4))
  axis(1, at = 1:9, labels = (1:9)^2)
  curve(plogis(cbind(1,x,0, 0, mean_mean_elev, mean_cross_fall, 
                     avg_trait_vals$trait_pc1[i], 
                     avg_trait_vals$trait_pc2[i], avg_trait_vals$trait_pc3[i], 
                     mean_date_diff_b2r, 0*0, 
                     x*avg_trait_vals$trait_pc3[i])%*%as.numeric(species_surv_coefs[i,])), 
        add=T, col=line_pal[2], lwd=2, lty=3)
  curve(plogis(cbind(1,x,0, 1, mean_mean_elev, mean_cross_fall, 
                     avg_trait_vals$trait_pc1[i], 
                     avg_trait_vals$trait_pc2[i], avg_trait_vals$trait_pc3[i],
                     mean_date_diff_b2r, 0*1, 
                     x*avg_trait_vals$trait_pc3[i])%*%as.numeric(species_surv_coefs[i,])), 
        add=T, col=line_pal[1], lwd=2,lty=3)
  curve(plogis(cbind(1,x,0, 0, mean_mean_elev, mean_cross_fall, 
                     avg_trait_vals$trait_pc1[i], 
                     avg_trait_vals$trait_pc2[i], avg_trait_vals$trait_pc3[i], 
                     mean_date_diff_b2r, 0*0, 
                     x*avg_trait_vals$trait_pc3[i])%*%fixef(surv_mod_simp)), 
        add=T, col=line_pal[2], lwd=2)
  curve(plogis(cbind(1,x,0, 1, mean_mean_elev, mean_cross_fall, 
                     avg_trait_vals$trait_pc1[i], 
                     avg_trait_vals$trait_pc2[i], avg_trait_vals$trait_pc3[i], 
                     mean_date_diff_b2r, 0*1, 
                     x*avg_trait_vals$trait_pc3[i])%*%fixef(surv_mod_simp)), 
        add=T, col=line_pal[1], lwd=2)
}
plot.new() #blank plot to add legend to
plot.new() #blank plot to add legend to
plot.new() #blank plot to add legend to
legend('top', title = expression(Model~coefficients), legend = c("Fixed","Random"), 
       col = "black", lwd = 2:1, lty=1:3, xpd = TRUE, horiz = FALSE, 
       cex = 1.5, seg.len=2, bty='n', y.intersp=1)
legend('bottom', title = expression(Organic~amendment~applied), legend = c("YES","NO"), 
       col = c(line_pal[1],line_pal[2]), lwd = 2, xpd = TRUE, horiz = FALSE, 
       cex = 1.5, seg.len=2, bty='n',y.intersp=1)
dev.off()


######## ######## ######## ######## ######## ######## ######## ######## ######## 
######## FIGURE 3: Species survival for sig. effects ######## 
######## ######## ######## ######## ######## ######## ######## ######## ######## 
mean_sqrt_height <- mean(sqrt(survival_data$height_baseline_cm), na.rm=T) #4.16036
mean_trait_pc1 <- mean(survival_data$trait_pc1, na.rm=T) #0.05387055
mean_trait_pc2 <- mean(survival_data$trait_pc2, na.rm=T) #0.07254566
mean_trait_pc3 <- mean(survival_data$trait_pc3, na.rm=T) #0.2548202

avg_start_height <- survival_data %>%
  group_by(species) %>%
  summarise(avg_height = mean(height_baseline_cm, na.rm=TRUE))
mean_height <- avg_trait_vals %>%
  left_join(avg_start_height, by="species")

# Get error bars/ confidence intervals on random effect coeff estimates
# Error bars (arrows)
surv_error_bars <- logistic.calcs(y = survival_data$alive_resurvey, 
                                  groups = as.factor(survival_data$species), 
                                  se_mult = 1.96)
avg_trait_vals_error <- data.frame(avg_trait_vals, surv_error_bars)

# make a prediction with CIs
pc1_pred <- lmer.predict(mod=surv_mod_simp, 
                         newdat=data.frame(1,mean_sqrt_height,0, 0, mean_mean_elev, 
                                           mean_cross_fall, seq.func(survival_data$trait_pc1),
                                           mean_trait_pc2, mean_trait_pc3, 
                                           mean_date_diff_b2r, 0*0, 
                                           mean_sqrt_height*mean_trait_pc3), 
                         se.mult=1.96, binom=T, poisson=F)
pc2_pred <- lmer.predict(mod=surv_mod_simp, 
                         newdat=data.frame(1,mean_sqrt_height,0, 0, mean_mean_elev, 
                                           mean_cross_fall, mean_trait_pc1, 
                                           seq.func(survival_data$trait_pc2), 
                                           mean_trait_pc3, 
                                           mean_date_diff_b2r, 0*0, 
                                           mean_sqrt_height*mean_trait_pc3), 
                         se.mult=1.96, binom=T, poisson=F)
pc3_pred <- lmer.predict(mod=surv_mod_simp, 
                         newdat=data.frame(1,mean_sqrt_height,0, 0, mean_mean_elev, 
                                           mean_cross_fall, mean_trait_pc1, 
                                           mean_trait_pc2, seq.func(survival_data$trait_pc3), 
                                           mean_date_diff_b2r, 0*0, 
                                           mean_sqrt_height*seq.func(survival_data$trait_pc3)), 
                         se.mult=1.96, binom=T, poisson=F)

pdf(file="Outputs/Figures/figure_3.pdf", width=12, height=12)
par(mfrow=c(2,2))
#Fig. 3 (a) Survival~PC1
with(survival_data, plot(alive_resurvey ~ trait_pc1, type = "n", cex.lab=1.15,
                         xlab = expression(Trait~PC1), 
                         ylab = expression(Probability~of~Survival), 
                         xlim = c(-2.2,3.5)))
## plot the fitted line AND the CIs
plot.CI.func(x.for.plot = seq.func(survival_data$trait_pc1), pred=pc1_pred$y, 
             upper = pc1_pred$phi, lower = pc1_pred$plo, env.colour = "grey", 
             env.trans = 80, line.colour = "black", line.weight = 3, line.type = 2)
points(avg_trait_vals_error$trait_pc1, avg_trait_vals_error$mean.y, 
       col=spec_pal_rainbow[as.factor(avg_trait_vals_error$species)], 
       ylab=expression(Probability~of~Survival), pch=16, cex=2)
text(avg_trait_vals_error$trait_pc1, avg_trait_vals_error$mean.y, 
     labels=short_names, col=spec_pal_rainbow[as.factor(avg_trait_vals_error$species)], 
     cex=1, pos=c(rep(4,12),2)) ## keep tinkering with label positions
arrows(x0=avg_trait_vals_error$trait_pc1, y0=avg_trait_vals_error$lower.y, 
       x1=avg_trait_vals_error$trait_pc1, y1=avg_trait_vals_error$upper.y, code=3, 
       angle=90, length=0, col=spec_pal_rainbow[as.factor(avg_trait_vals_error$species)])
mtext(side=3, adj=0, line=0.25, text="(a)", cex=1.25)

#Fig. 3 (b) Survival~PC2 
with(survival_data, plot(alive_resurvey ~ trait_pc2, type = "n", cex.lab=1.15, 
                         xlab = expression(Trait~PC2), 
                         ylab = expression(Probability~of~Survival), 
                         xlim = c(-3.2,1.8))) 
plot.CI.func(x.for.plot = seq.func(survival_data$trait_pc2), pred=pc2_pred$y, 
             upper = pc2_pred$phi, lower = pc2_pred$plo, env.colour = "grey", 
             env.trans = 80, line.colour = "black", line.weight = 3, line.type = 2)
points(avg_trait_vals_error$trait_pc2, avg_trait_vals_error$mean.y, 
       col=spec_pal_rainbow[as.factor(avg_trait_vals_error$species)], 
       ylab=expression(Probability~of~Survival), pch=16, cex=2)
text(avg_trait_vals_error$trait_pc2, avg_trait_vals_error$mean.y, 
     labels=short_names, 
     col=spec_pal_rainbow[as.factor(avg_trait_vals_error$species)], cex=1, 
     pos=c(rep(2,5),1,2,4,3,4,rep(2,2),4))
arrows(x0=avg_trait_vals_error$trait_pc2, y0=avg_trait_vals_error$lower.y, 
       x1=avg_trait_vals_error$trait_pc2, y1=avg_trait_vals_error$upper.y, code=3, 
       angle=90, length=0, col=spec_pal_rainbow[as.factor(avg_trait_vals_error$species)])
mtext(side=3, adj=0, line=0.25, text="(b)", cex=1.25)

specs_to_plot<-c("cote", "coin", "eumo")

#Fig. 3 (c) Survival~PC3
with(subset(survival_data, species %in%specs_to_plot), plot(jitter(alive_resurvey, amount=0.02) ~ 
                                                              jitter(sqrt(height_baseline_cm), amount=0.05), 
                                                            col = spec_pal_rainbow[c(6,5,9)][as.factor(species)], pch=16, cex=1.5, cex.lab=1.15, xaxt="n",
                                                            xlab = expression(Initial~Height~(cm)), 
                                                            ylab = expression(Probability~of~Survival)))
axis(1, at = c(1:8), labels = c(1:8)^2)
height_limits<-survival_data%>%
  group_by(species)%>%
  summarise(min_sqrt_height = min(sqrt(height_baseline_cm), na.rm=T),
            max_sqrt_height = max(sqrt(height_baseline_cm), na.rm=T))
spec_names<-row.names(species_surv_coefs)

for(i in c(6,5,9)){
  curve(plogis(cbind(1,x,0, 0, mean_mean_elev, mean_cross_fall, 0,
                     0, mean_height$trait_pc3[i], mean_date_diff_b2r,
                     0*0, x*mean_height$trait_pc3[i]) %*% fixef(surv_mod_simp)), 
        add=T, from=as.numeric(height_limits[i,2]), to=as.numeric(height_limits[i,3]),
        col=spec_pal_rainbow[i], lwd = 3)
}

text(x=1, y=c(0.4, 0.3, 0.2, 0.1), labels =c("Effect of trait PC3:", "C. tessellaris (small leaves, high SLA)",
                                             "C. intermedia (average leaves, average SLA)",
                                             "E. moluccana (large leaves, low SLA)"),
     col=c("black", spec_pal_rainbow[c(6,5,9)]), pos=4)
mtext(side=3, adj=0, line=0.25, text="(c)", cex=1.25)

#Fig. 3 (d) Survival~Mean Elev
with(survival_data, plot(jitter(alive_resurvey, amount=0.02) ~ 
                           jitter(mean_elev, amount=0.05), 
                         col = alpha(spec_pal_rainbow, 0.4)[as.factor(species)], pch=16, cex=1.5, cex.lab=1.15,
                         xlab = expression(Mean~Elevation~(m)), xlim=c(83,95), xaxt="n", 
                         ylab = expression(Probability~of~Survival)))
axis(1, at = 86:95, labels = 86:95)
curve(plogis(cbind(1,mean_sqrt_height,0, 0, x, mean_cross_fall, mean_trait_pc1, 
                   mean_trait_pc2, mean_trait_pc3, mean_date_diff_b2r, 
                   0*0, mean_sqrt_height*mean_trait_pc3)%*%fixef(surv_mod_simp)), 
      add=T, from=86, col="black", lwd=3)
for(i in 1:dim(species_surv_coefs)[1]){
  curve(plogis(cbind(1,mean_sqrt_height,0, 0, x, mean_cross_fall, mean_height$trait_pc1[i],
                     mean_height$trait_pc2[i], mean_height$trait_pc3[i], mean_date_diff_b2r,
                     0*0, mean_sqrt_height*mean_height$trait_pc3[i]) %*%
                 as.numeric(species_surv_coefs[i,])), add=T, from=86, col=spec_pal_rainbow[i], 
        lwd = 1)
}
text(x=c(85.05269, 85.32729, 84.52706, 84.26031, 85.13899, 85.20961, 85.38220, 
         85.08408, 85.13115, 85.13115, 84.84087, 85.12331, 85.30375), 
     y=c(0.5292175,0.5975983,0.8544236,0.9064753,0.8047402,0.8695267,0.7075284,
         0.8344063,0.9759116,0.9592626,0.8864581,0.7940035,0.8209330), cex=0.75,
     labels=short_names,col=spec_pal_rainbow[as.factor(avg_trait_vals_error$species)])
mtext(side=3, adj=0, line=0.25, text="(d)", cex=1.25)
dev.off()


######## ######## ######## ######## ######## ######## ######## ######## ######## 
######## FIGURE: Species growth  ######## 
######## ######## ######## ######## ######## ######## ######## ######## ######## 

species_growth_coefs <- coef(growth_mod_simp)$species
species_growth_coefs[1,] ##this is to know what to put for the cbind. 

pdf(file="Outputs/Figures/figure_4.pdf", width=12, height=12)
par(mfrow=c(4,4))
for(i in 1:length(species_names)){
  species_trait_data <- droplevels(subset(pos_growth_data, 
                                          species == species_names[i]))
  with(species_trait_data, 
       plot(sqrt(growth_mm_day) ~ sqrt(height_baseline_cm), 
            col= alpha(line_pal[2], 0.5), pch = 16, xaxt = "n", yaxt="n", 
            xlab = expression(Initial~Height~(cm)), xlim = c(0, 9), ylim=c(0,2.5),
            ylab=expression(Height~Growth~"(mm/day)"), main = scientific_names[i], 
            cex.main=1.4, cex.lab=1.15))
  axis(1, at = 0:9, labels = (0:9)^2)
  axis(2, at = seq(0, 2.5, 0.25), labels = seq(0, 2.5, 0.25)^2)
  curve(cbind(1,x,0, 0.5, mean_mean_elev, mean_cross_fall,  
              avg_trait_vals$trait_pc1[i], 
              avg_trait_vals$trait_pc2[i], avg_trait_vals$trait_pc3[i], 
              mean_date_diff_b2r, 0*0)%*%as.numeric(species_growth_coefs[i,]), 
        add=T, col=line_pal[1], lwd=2, lty=3)
  curve(cbind(1,x,0, 0.5, mean_mean_elev, mean_cross_fall,  
              avg_trait_vals$trait_pc1[i], 
              avg_trait_vals$trait_pc2[i], avg_trait_vals$trait_pc3[i],
              mean_date_diff_b2r, 0*0)%*%fixef(growth_mod_simp), 
        add=T, col=line_pal[1], lwd=2)
}
plot.new() #blank plot to add legend to
plot.new() #blank plot to add legend to
plot.new() #blank plot to add legend to
legend('center', title = expression(Model~coefficients), legend = c("Fixed","Random"), 
       col = line_pal[1], lwd = 2, lty=1:3, xpd = TRUE, horiz = FALSE, 
       cex = 1.5, seg.len=3, bty='n')
dev.off()

####### ####### ####### ####### ####### ####### ####### ####### ####### #######
####### FIGURE 5: Provenance survival ########
####### ####### ####### ####### ####### ####### ####### ####### ####### #######
# eute, eute_dry, eucr, eucr_dry only
# simple model for prov as have already explored full data
# we can't let height interact with species and prov because it varies 
# systematically with species and provenance
prov_surv_mod <- glmer(alive_resurvey ~ sqrt(height_baseline_cm) + species + 
                         provenance + species:provenance + (1|plot), 
                       data = provenance_surv_data, family = binomial(link = "logit"))
summary(prov_surv_mod)

prov_growth_mod <- lmer(sqrt(growth_mm_day) ~ sqrt(height_baseline_cm) + species + 
                          provenance + species:provenance + (1|block/plot), 
                        provenance_growth_data)
summary(prov_growth_mod)

prov_table <- tab_model(prov_surv_mod, prov_growth_mod, show.re.var= TRUE, 
          dv.labels= c("Survival Model", "Height Growth Model"),
          pred.labels =c("(Intercept)", 
                         "Initial height [sqrt]",
                         "Species [E. teretecornis]",
                         "Provenance [dry]",
                         "Species [E. teretecornis] : provenance [dry]"),
          string.pred = "Coeffcient",
          string.ci = "CI (95%)",
          string.p = "P-Value",
          transform = NULL,
          file = "Outputs/Tables/prov_table.html")
prov_table
webshot("Outputs/Tables/prov_table.html", "Outputs/Tables/prov_table.pdf") 
  
## Creating the provenance figure objects... 
# include ~ provenance | species + sqrt(height_baseline_cm))  
prov_surv_mod_plot <- plot(emmeans(prov_surv_mod, ~ provenance | species), 
                           horizontal=F, type = "response", 
                           colors = spec_pal_rainbow[4], 
                           xlab = expression(Probability~of~Survival), 
                           ylab = expression(Provenance)) 

## to see plotting options for emmeans objects look at this ?plot.emmGrid
prov_growth_mod_plot<-plot(emmeans(prov_growth_mod, ~ provenance | species, 
                                   pbkrtest.limit = 4782, lmerTest.limit = 4782), 
                           horizontal=F, type = "response", 
                           colors = spec_pal_rainbow[4], 
                           xlab = expression(Height~Growth~(mm/day)), 
                           ylab = expression(Provenance))
ref_grid(prov_growth_mod)
ref_grid(prov_surv_mod)
#estimated marginal means held at average height of 20cm

## p values for tests in plots
test<-emmeans(prov_surv_mod, ~ provenance | species)
emmeans::contrast(test)

test2<-emmeans(prov_growth_mod, ~ provenance | species, 
               pbkrtest.limit = 4782, lmerTest.limit = 4782)
emmeans::contrast(test2)

# Plotting the actual figures
spp <- prov_surv_mod_plot +
  theme_bw() + theme(plot.background = element_blank(), 
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank()) + 
  theme(axis.line = element_line(color = "black")) + 
  facet_wrap(~ species, labeller = custom_labels, strip.position = "right", 
             scales = "free", ncol = 1)

gpp <- prov_growth_mod_plot + 
  theme_bw() + theme(plot.background = element_blank(), 
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank()) + 
  theme(axis.line = element_line(color = "black")) +
  facet_wrap(~ species, labeller = custom_labels, strip.position = "right", 
             scales = "free", ncol = 1)

anno <- data.frame(
  label = c("p = 0.0102", "p = 0.0001"),
  species = c("eucr","eute"),
  x = c(0.95, 1.01),
  y = c(1.5,1.5),
  x1 = c(0.93, 0.997), 
  x2 = c(0.92, 0.99), 
  y1 = c(1, 1), 
  y2 = c(2, 2))

anno2 <- data.frame(
  label = c("p = 0.5666", "p = 0.2705"),
  species   = c("eucr","eute"),
  x     = c(0.945, 1.87),
  y     = c(1.5,1.5),
  x1 = c(0.92, 1.82), 
  x2 = c(0.90, 1.8), 
  y1 = c(1, 1), 
  y2 = c(2, 2))

spp_final <- spp + geom_text(data = anno, mapping = aes(x = x, y = y, label = label)) + 
  geom_segment(data = anno, mapping = aes(x = x1, xend = x1, y = y1, yend = y2), 
               colour = "black") + 
  geom_segment(data = anno, mapping = aes(x = x1, xend = x2, y = y1, yend = y1), 
               colour = "black") +
  geom_segment(data = anno, mapping = aes(x = x1, xend = x2, y = y2, yend = y2),
               colour = "black") + 
  geom_segment(data = anno, mapping = aes(x = x1, xend = x2, y = y2, yend = y2),
               colour = "black")  +
geom_segment(data = anno, mapping = aes(x = x1, xend = x2, y = y2, yend = y2),
             colour = "black") + 
  geom_segment(data = anno, mapping = aes(x = x1, xend = x2, y = y2, yend = y2),
               colour = "black")  

gpp_final <- gpp + geom_text(data = anno2, mapping = aes(x = x, y = y, label = label)) + 
  geom_segment(data = anno2, mapping = aes(x = x1, xend = x1, y = y1, yend = y2), 
               colour = "black") + 
  geom_segment(data = anno2, mapping = aes(x = x1, xend = x2, y = y1, yend = y1), 
               colour = "black") +
  geom_segment(data = anno2, mapping = aes(x = x1, xend = x2, y = y2, yend = y2),
               colour = "black") + 
  geom_segment(data = anno2, mapping = aes(x = x1, xend = x2, y = y2, yend = y2),
               colour = "black")  +
  geom_segment(data = anno2, mapping = aes(x = x1, xend = x2, y = y2, yend = y2),
               colour = "black") + 
  geom_segment(data = anno2, mapping = aes(x = x1, xend = x2, y = y2, yend = y2),
               colour = "black")  

ggarrange(spp_final, gpp_final, font.label = list(size = 11), labels =c("(a)","(b)"), ncol = 1, nrow = 2 )
ggsave("Outputs/Figures/prov_plots.pdf", width = 3, height = 8) 

######## ######## ######## ######## ######## ######## ######## ######## ######## 
######## Simple height only models ######## 
######## ######## ######## ######## ######## ######## ######## ######## ######## 

species_surv_coefs_htonly <- coef(surv_mod_htonly)$species
species_surv_coefs_htonly[1,] ##this is to know what to put for the cbind. 

pdf(file="Outputs/Figures/figure_surv_mod_htonly.pdf", width=12, height=12)
par(mfrow=c(4,4))
for(i in 1:length(species_names)){
  species_trait_data<-droplevels(subset(survival_data, species==species_names[i]))
  with(species_trait_data, plot(jitter(alive_resurvey, amount=0.05) ~ sqrt(height_baseline_cm),
                                col=alpha(line_pal[2],0.5), pch=16, xlim = c(0, 9),
                                xaxt="n", cex.lab=1.15, 
                                xlab = expression(Initial~Height~(cm)), 
                                ylab=expression(Probability~of~Survival), 
                                main=scientific_names[i], cex.main=1.4))
  axis(1, at = 0:9, labels = (0:9)^2)
  curve(plogis(cbind(1,x,x^2)%*%as.numeric(species_surv_coefs_htonly[i,])), 
        add=T, col=line_pal[1], lwd=2, lty=3)
  curve(plogis(cbind(1,x,x^2)%*%fixef(surv_mod_htonly)), 
        add=T, col=line_pal[1], lwd=2)
}
plot.new() #blank plot to add legend to
plot.new() #blank plot to add legend to
plot.new() #blank plot to add legend to
legend('top', title = expression(Model~coefficients), legend = c("Fixed","Random"), 
       col = line_pal[1], lwd = 2:1, lty=1:3, xpd = TRUE, horiz = FALSE, 
       cex = 1.5, seg.len=2, bty='n', y.intersp=1)
dev.off()

species_growth_coefs_htonly <- coef(growth_mod_htonly)$species
species_growth_coefs_htonly[1,] ##this is to know what to put for the cbind. 

pdf(file="Outputs/Figures/figure_growth_mod_htonly.pdf", width=12, height=12)
par(mfrow=c(4,4))
for(i in 1:length(species_names)){
  species_trait_data <- droplevels(subset(pos_growth_data, 
                                          species == species_names[i]))
  with(species_trait_data, 
       plot(sqrt(growth_mm_day) ~ sqrt(height_baseline_cm), 
            col= alpha(line_pal[2], 0.5), pch = 16, xaxt = "n", yaxt = "n", 
            xlab = expression(Initial~Height~(cm)), ylim = c(0, 2.5), xlim = c(0, 9),
            ylab=expression(Height~Growth~"(mm/day)"), main = scientific_names[i], 
            cex.main=1.4, cex.lab=1.15))
  axis(1, at = 1:9, labels = (1:9)^2)
  axis(2, at = seq(0, 2.5, 0.25), labels = seq(0, 2.5, 0.25)^2)
  curve(cbind(1,x,x^2)%*%as.numeric(species_growth_coefs_htonly[i,]), 
        add=T, col=line_pal[1], lwd=2, lty=3)
  curve(cbind(1,x,x^2)%*%fixef(growth_mod_htonly), 
        add=T, col=line_pal[1], lwd=2)
}
plot.new() #blank plot to add legend to
plot.new() #blank plot to add legend to
plot.new() #blank plot to add legend to
legend('center', title = expression(Model~coefficients), legend = c("Fixed","Random"), 
       col = line_pal[1], lwd = 2, lty=1:3, xpd = TRUE, horiz = FALSE, 
       cex = 1.5, seg.len=3, bty='n')
dev.off()


##### Overall probabilities for survival/ growth etc for write up
treatment_survival <- survival_data %>%
  mutate(treatment = ifelse(compost == "YES" & innoculation == "NO", 
                            "compost", 
                            ifelse(compost == "NO" & innoculation == "YES", "innoc", 
                                   ifelse(compost == "NO" & innoculation == "NO", "neither", 
                                          "both"))))
treatment_surv_means <- treatment_survival %>%
  group_by(treatment) %>%
  summarise(survival = mean(alive_resurvey))

treatment_error_bars <- logistic.calcs(y = treatment_survival$alive_resurvey, 
                                       groups = as.factor(treatment_survival$treatment), 
                                       se_mult = 1.96)
treatment_probs <- data.frame(unique(treatment_survival$treatment), treatment_error_bars)

treatment_error_bars_species <- logistic.calcs(y = treatment_survival$alive_resurvey, 
                                       groups = as.factor(treatment_survival$species), 
                                       se_mult = 1.96)
treatment_probs_sp <- data.frame(sort(unique(treatment_survival$species)), treatment_error_bars_species)
survival_pm <- (treatment_probs_species$upper.y - treatment_probs_species$lower.y)/2
treatment_probs_species <- data.frame(sort(unique(treatment_survival$species)), treatment_error_bars_species, 
                                      survival_pm)

treatment_growth <- pos_growth_data %>%
  mutate(treatment = ifelse(compost == "YES" & innoculation == "NO", 
                            "compost", 
                            ifelse(compost == "NO" & innoculation == "YES", "innoc", 
                                   ifelse(compost == "NO" & innoculation == "NO", "neither", 
                                          "both"))))

treatment_growth_means <- treatment_growth %>%
  group_by(treatment) %>% 
  summarise(growth=mean(growth_mm_day))

treatment_growth_means_spec <- treatment_growth %>%
  group_by(species) %>% 
  summarise(growth=mean(growth_mm_day))

## The better way of estimating species average survival which was not used for thesis table
simple_surv_mod <- glmer(alive_resurvey ~ 1 + (1|plot) + (1|species), 
                         data = survival_data, family = binomial(link = "logit"))
avg_species_surv_coefs <- mean(plogis(coef(simple_surv_mod)$species[,1]), rm.na=T)
simple_species_surv_coefs <- plogis(coef(simple_surv_mod)$species[,1])
summary(simple_surv_mod)

simple_growth_mod <- lmer(sqrt(growth_mm_day) ~ 1 + (1|plot) + (1|species),
                          pos_growth_data)
summary(simple_growth_mod)
avg_species_growth_coefs <- coef(simple_growth_mod)$species
mean(avg_species_growth_coefs[,1], rm.na=T)

simple_surv_treatment<- glmer(alive_resurvey ~ 1 + (1|plot) + (1|species) + (1|treatment), 
                         data = treatment_survival, family = binomial(link = "logit"))
avg_treatment_surv_coefs <- mean(plogis(coef(simple_surv_treatment)$species[,1]), rm.na=T)
simple_treatment_surv_coefs <- plogis(coef(simple_surv_treatment)$treatment[,1])
summary(simple_surv_treatment)

simple_growth_treatment <- lmer(sqrt(growth_mm_day) ~ 1 + (1|plot) + (1|species) +(1|treatment),
                                treatment_growth)
summary(simple_growth_treatment)
avg_species_growth_coefs <- coef(simple_growth_treatment)$species
mean(avg_species_growth_coefs[,1], rm.na=T)


######## ######## ######## ######## ######## ######## ######## ######## ######## 
######## FIGURE: Growth vs Trait PC1-3 ######## 
######## ######## ######## ######## ######## ######## ######## ######## ######## 
pc1_pred_g <- lmer.predict(mod=growth_mod_simp, 
                           newdat=data.frame(1,mean_sqrt_height,0, 0, mean_mean_elev, 
                                             mean_cross_fall, seq.func(pos_growth_data$trait_pc1),
                                             mean_trait_pc2, mean_trait_pc3, 
                                             mean_date_diff_b2r, 0*0), 
                           se.mult=1.96, binom=F, poisson=F)
pc2_pred_g <- lmer.predict(mod=growth_mod_simp, 
                           newdat=data.frame(1,mean_sqrt_height,0, 0, mean_mean_elev, 
                                             mean_cross_fall, mean_trait_pc1, 
                                             seq.func(pos_growth_data$trait_pc2), 
                                             mean_trait_pc3, 
                                             mean_date_diff_b2r, 0*0), 
                           se.mult=1.96, binom=F, poisson=F)
pc3_pred_g <- lmer.predict(mod=growth_mod_simp, 
                           newdat=data.frame(1,mean_sqrt_height,0, 0, mean_mean_elev, 
                                             mean_cross_fall, mean_trait_pc1, 
                                             mean_trait_pc2, seq.func(pos_growth_data$trait_pc3), 
                                             mean_date_diff_b2r, 0*0), 
                           se.mult=1.96, binom=F, poisson=F)

pdf(file="Outputs/Figures/figure_growthvtrait.pdf", width=12, height=12)
par(mfrow=c(2,2))
#(a) 
with(species_trait_data, 
     plot(sqrt(growth_mm_day) ~ trait_pc1, type="n", 
          cex.main=1.4, cex.lab=1.15, yaxt="n", 
          xlab = expression(Trait~PC1), 
          ylab=expression(Height~Growth~"(mm/day)"), 
          xlim = c(-2.4, 3.6), ylim=c(0,1.5)))
axis(2, at = seq(0, 1.5, 0.25), labels = seq(0, 1.5, 0.25)^2)
plot.CI.func(x.for.plot = seq.func(pos_growth_data$trait_pc1), pred=pc1_pred_g$y, 
             upper = pc1_pred_g$phi, lower = pc1_pred_g$plo, env.colour = "grey", 
             env.trans = 80, line.colour = "black", line.weight = 3, line.type = 1)
points(avg_trait_vals_error$trait_pc1, sqrt(treatment_growth_means_spec$growth), 
       col=spec_pal_rainbow[as.factor(avg_trait_vals_error$species)], pch=16, cex=2) 
text(avg_trait_vals_error$trait_pc1, sqrt(treatment_growth_means_spec$growth), 
     labels=short_names, col=spec_pal_rainbow[as.factor(avg_trait_vals_error$species)], 
     cex=1, pos=c(4,1,4,3,4,4,3,3,4,2,2,4,1)) ## keep tinkering with label positions
mtext(side=3, adj=0, line=0.25, text="(a)", cex=1.25)
#(b)  
with(species_trait_data, 
     plot(sqrt(growth_mm_day) ~ trait_pc2, type="n", 
          cex.main=1.4, cex.lab=1.15, yaxt="n", 
          xlab = expression(Trait~PC2), 
          ylab=expression(Height~Growth~"(mm/day)"), 
          xlim = c(-3.3, 1.7), ylim=c(0,1.5)))
axis(2, at = seq(0, 1.5, 0.25), labels = seq(0, 1.5, 0.25)^2)
plot.CI.func(x.for.plot = seq.func(pos_growth_data$trait_pc2), pred=pc2_pred_g$y, 
             upper = pc2_pred_g$phi, lower = pc2_pred_g$plo, env.colour = "grey", 
             env.trans = 80, line.colour = "black", line.weight = 3, line.type = 1)
points(avg_trait_vals_error$trait_pc2, sqrt(treatment_growth_means_spec$growth), 
       col=spec_pal_rainbow[as.factor(avg_trait_vals_error$species)], 
       pch=16, cex=2)
text(avg_trait_vals_error$trait_pc2, sqrt(treatment_growth_means_spec$growth), 
     labels=short_names, col=spec_pal_rainbow[as.factor(avg_trait_vals_error$species)], 
     cex=1, pos=c(2,2,2,2,2,1,4,4,2,2,2,3,1)) ## keep tinkering with label positions
mtext(side=3, adj=0, line=0.25, text="(b)", cex=1.25)
#(c)  
with(species_trait_data, 
     plot(sqrt(growth_mm_day) ~ trait_pc3, type="n", 
          cex.main=1.4, cex.lab=1.15, yaxt="n", 
          xlab = expression(Trait~PC3), 
          ylab=expression(Height~Growth~"(mm/day)"), 
          xlim = c(-2.2, 1.65), ylim=c(0,1.5)))
axis(2, at = seq(0, 1.5, 0.25), labels = seq(0, 1.5, 0.25)^2)
plot.CI.func(x.for.plot = seq.func(pos_growth_data$trait_pc3), pred=pc3_pred_g$y, 
             upper = pc3_pred_g$phi, lower = pc3_pred_g$plo, env.colour = "grey", 
             env.trans = 80, line.colour = "black", line.weight = 3, line.type = 2)
points(avg_trait_vals_error$trait_pc3, sqrt(treatment_growth_means_spec$growth), 
       col=spec_pal_rainbow[as.factor(avg_trait_vals_error$species)], 
       pch=16, cex=2)
text(avg_trait_vals_error$trait_pc3, sqrt(treatment_growth_means_spec$growth), 
     labels=short_names, col=spec_pal_rainbow[as.factor(avg_trait_vals_error$species)], 
     cex=1, pos=c(3,4,4,3,1,4,4,1,2,2,2,2,1)) ## keep tinkering with label positions
mtext(side=3, adj=0, line=0.25, text="(c)", cex=1.25)
dev.off()

notseedlingheight <- survival_data %>%
  group_by(species) %>%
  filter(height_baseline_cm > 50)





