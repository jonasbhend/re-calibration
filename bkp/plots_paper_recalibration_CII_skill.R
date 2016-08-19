library(myhelpers)
library(ggplot2, lib='~/R/oldlib')
library(scales)
library(sp)
library(raster)
library(RColorBrewer)
library(dplyr)
library(tidyr)

plot_boxplot <- function(data, group1="index", group2="method", ...){
  scores <- unique(data$score)
  inits <- unique(data$initmon)
  group1s <- unique(data[[group1]])
  group2s <- unique(data[[group2]])
  data[[group2]] <- factor(as.character(data[[group2]]), 
                           unique(as.character(data[[group2]])))
  data[[group1]] <- factor(as.character(data[[group1]]), 
                           unique(as.character(data[[group1]])))
  
  ylimtmp <- data %>% group_by(score) %>% summarise(ymin=min(value, na.rm=T), ymax=max(value, na.rm = T))
  ylims <- as.matrix(ylimtmp[,-1])
  rownames(ylims) <- as.character(ylimtmp$score)
  if ("EnsCorr" %in% scores) ylims['EnsCorr',] <- c(-1,1)
  if ("FairCrpss" %in% scores) ylims['FairCrpss',1] <- -1
  
  xat <- outer(seq(along=group2s), (seq(along=group1s) - 1)*(length(group2s) + 2), '+')
  
  j <- 0
  for (tscore in scores){
    for (tinit in inits){
      j <- j + 1
      plot(0, type='n', xlim=range(xat) + c(-0.5,0.5), ylim=ylims[tscore,], axes=FALSE, 
           las=1, log=if (tscore == 'EnsRmse') 'y' else '')
      grid(nx=NA, ny=NULL, equilogs=FALSE)
      if (tscore != "EnsRmse"){
        abline(h=0, lty=2)
      } 
      if (tinit == '05') {
        axis(2, at= mean(if (tscore == 'EnsRmse') exp(par('usr')[3:4]) else par('usr')[3:4]), 
             label = relabeller("score", tscore), line=2, tick=F,
             cex.axis=par("cex.lab"))
        axis(2, las=1)
      }
      boxplot(as.formula(paste("value ~ ", group2, "*", group1)), axes=FALSE,
              data=data, at=xat, add=TRUE,
              las=1, lty=1, pch=16, cex=0.5, 
              col=ifelse(group2 == "method",methcols, vericols))
      box()
      if (tscore == scores[length(scores)]){
        axis(1, at=colMeans(xat), relabeller(group1, group1s), tick=F, cex.axis=par("cex.lab"))
      }
      axis(3, at=par("usr")[1], paste0(letters[j], ') ', relabeller("score", tscore), ' of forecasts for ', c(`11`='DJF (Nov. init.)', `05`="JJA (May init.)")[tinit]), hadj=0, line=0, cex.axis=par('cex.lab'), tick=F)
      if (tinit == inits[1] & tscore == scores[1]){
        legend('bottom', relabeller(group2, group2s), 
               fill=ifelse(group2 == 'method', methcols, vericols), 
               ncol=2, cex=par("cex.axis"), inset=0.02, bg='white')
      }
    }
  }
}

relabeller <- function(var, value){
  value <- as.character(value)
  if (var=="ccr") { 
    value[value == "TRUE"] <- "With recalibration on seasonal index"
    value[value == "FALSE"]   <- "Without recalibration"
  } else if (var == 'initmon'){
    value[value == "05"] <- "Summer (JJA, May init.)"
    value[value == '11'] <- "Winter (DJF, Nov. init.)"
  } else if (var == 'method'){
    value[grep('none', value)] <- 'no calibration'
    value <- gsub('smooth', 'debias', value)
    value <- gsub('comb', 'all', value)
    value <- gsub("fastqqmap", 'quantile mapping', value)
    value <- gsub('-forward', '', value)
  } else if (var == 'verimethod'){
    value[value == 'none'] <- 'no calibration'
    value[value == 'insample'] <- 'in sample'
    value[value == 'crossval1'] <- 'cross-validation (1 year)'
    value[value == 'crossval10'] <- 'cross-validation (10 years)'
    value[value == 'forward'] <- 'forward calibration'  
  } else if (var == 'score'){
    value[value == 'EnsCorr'] <- 'correlation'
    value[value == 'EnsRmse'] <- 'RMSE'
    value[value == 'FairCrpss'] <- 'CRPSS'
    value[value == 'FairSprErr'] <- 'spread to error ratio'
  } 
  return(value)
}


## read in scores from all indices
scores <- c("EnsCorr", "FairCrpss", "FairCrpss.sigma", "FairSprErr")
indices <- c("tas", "ITV", "PDD90", "PDD98", "NDD10", "NDD02", "HDDch", "CDD")
initmons <- c("05", "11")
methods <- c("none", paste0(outer(c("fastqqmap", "smooth", "comb"), c("", "Recal"), paste0), '-forward_1981-2014_ERA-INT'))

slong <- read_scores(ind=indices, init=initmons, 
                     scores=scores, method=methods, 
                     obs='ERA-INT', period='1981-2014',
                     detrend=FALSE, ccr=c(FALSE, TRUE), 
                     cleanup=TRUE)

## convert file names
slong <- slong %>% 
  mutate(index=factor(index, unique(slong$index)),
         index2=factor(gsub("PDD90|NDD10", "DD90/10", 
                            gsub("PDD98|NDD02", "DD98/02", 
                                 gsub("HDDch|CDD", "HDD/CDD", index))), 
                       c("tas", "ITV", "DD90/10", "DD98/02", "HDD/CDD")),
         lon=ifelse(lon > 180, lon - 360, lon))

ggplot(filter(slong, score %in% c("EnsCorr", "FairCrpss"), 
              value > -2, lsm > 0.5, lat > -60, ! ccr), 
       aes(x=index,y=value, fill=as.factor(method))) + 
  geom_boxplot(outlier.size=0.6) + 
  facet_grid(score ~ initmon, scales='free')


save(slong, file="data/all_scores.Rdata")

q(save='no')

######################################################################################################
## compare the different verification strategies
######################################################################################################

dpath <- '/store/msclim/bhendj/EUPORIAS/skill_scores/global2/seasonal'
indices <- 'tas'
## verimethods <- c('none', 'smooth', 'smooth-crossval1', 'smooth-crossval10', 'smooth-forward')
verimethods <- c('none', 'insample', 'crossval1', 'crossval10', 'forward')
methods <- c('smooth', 'comb', 'fastqqmap')
vericols <- c(grey(0.6), hcl(c(120,200,240,10), l=60, c=60))
names(vericols) <- verimethods
## methods <- gsub('-forward', '', methods)
initmonths <- c("05", "11")
scores <- c('EnsCorr', "EnsRmse", "FairCrpss")

## set up script for 
vmethods <- gsub("-insample", "", c('none', outer(methods, setdiff(verimethods, 'none'), paste, sep='-')))

## read in skill scores
skill <- read_scores(model='ecmwf-system4',
                     method=vmethods,
                     score=scores,
                     index=indices,
                     init=initmonths,
                     grid='global2',
                     lead=2,
                     ccr=FALSE,
                     detrend=FALSE)
skill$method <- factor(gsub('_1.*', '', skill$method))


## adjust longitudes
subskill <- skill %>% filter(lsm > 0.5 & lat > -60) %>%
  mutate(lon2 = ifelse(lon > 188, lon - 360, lon),
         verimethod = ifelse(method %in% methods, 'insample', gsub('.*-', "", method)),
         method = ifelse(method == "none", methods[1], gsub("-.*", "", method)), 
         value2 = ifelse(score == 'EnsRmse', log(pmax(value, 0)), value))
for (meth in methods[-1]){
  subskill <- rbind(subskill, subskill %>% filter(method == methods[1] & verimethod == "none") %>% mutate(method=meth))
}
## reorder factors
subskill$verimethod <- factor(as.character(subskill$verimethod), verimethods)
subskill$method <- factor(as.character(subskill$method), methods)

ylimtmp <- subskill %>% group_by(score) %>% summarise(ymin=min(value), ymax=max(value))
ylims <- as.matrix(ylimtmp[,-1])
rownames(ylims) <- as.character(ylimtmp$score)
ylims['EnsCorr',] <- c(-1,1)
ylims['FairCrpss',] <- c(-1, 0.6)
xat <- outer(seq(verimethods), (seq(methods) - 1)*(length(verimethods) + 2), '+')

## revert to base graphics
png('/users/bhendj/doc/figures/recalibration/recalibration-compare_verification_strategies.png', width=8, height=10, units='in', res=200)
par(mfrow=c(3,2), mar=c(0.5, 0.5, 3, 0.5), oma=c(2,4,0,0), 
    tcl=0.3, mgp=c(1.5, 0.5, 0), cex.axis=1, cex.lab=1.4)
j <- 0
for (tscore in unique(subskill$score)){
  for (tinit in unique(subskill$initmon)){
    j <- j + 1
    plot(0, type='n', xlim=range(xat) + c(-0.5,0.5), ylim=ylims[tscore,], axes=FALSE, 
         las=1, log=if (tscore == 'EnsRmse') 'y' else '')
    grid(nx=NA, ny=NULL, equilogs=FALSE)
    if (tscore != "EnsRmse"){
      abline(h=0, lty=2)
    } 
    if (tinit == '05') {
      axis(2, at= mean(if (tscore == 'EnsRmse') exp(par('usr')[3:4]) else par('usr')[3:4]), 
           label = relabeller("score", tscore), line=2, tick=F,
           cex.axis=par("cex.lab"))
      axis(2, las=1)
    }
    boxplot(value ~ verimethod*method, axes=FALSE,
            data=filter(subskill, score == tscore & initmon == tinit), 
            las=1, lty=1, pch=16, cex=0.5, col=vericols, at=xat, add=TRUE)
    box()
    if (tscore == 'FairCrpss'){
      axis(1, at=colMeans(xat), relabeller("method", methods), tick=F, cex.axis=par("cex.lab"))
    }
    axis(3, at=par("usr")[1], paste0(letters[j], ') ', relabeller("score", tscore), ' of forecasts for ', c(`11`='DJF (Nov. init.)', `05`="JJA (May init.)")[tinit]), hadj=0, line=0, cex.axis=par('cex.lab'), tick=F)
    if (tinit == '05' & tscore == 'EnsCorr'){
      legend('bottom', relabeller("verimethod", verimethods), fill=vericols, ncol=2, cex=par("cex.axis"), inset=0.02, bg='white')
    }
  }
}
dev.off()


######################################################################################################
## compare the different calibration methods
######################################################################################################
dpath <- '/store/msclim/bhendj/EUPORIAS/skill_scores/global2/seasonal'
## indices <- c('tas', 'ITV', list.files(dpath, pattern='[PN]DD'))
indices <- c('tas', 'ITV', 'NDD10', 'NDD02', 'PDD90', 'PDD98', 'HDDch', 'CDD')
methods <- c(c('none', paste(c('smooth', 'trend', 'conditional', 'comb'), 'forward', sep='-'), paste(c('smooth', 'trend', 'conditional', 'comb'), 'Recal-forward', sep=''))[c(1,2,6,3,7,4,8,5,9)], "fastqqmap-forward")
submethods <- c(methods[-grep("Recal", methods)], "conditionalRecal-forward")
methcols <- c(grey(0.6), hcl(rep(1:5, each=2)/5*360, l=c(50,70), c=c(70,50)))
## methods <- c("none", paste0(c("unbias", "unbiasRecal", "comb", "combRecal", "fastqqmap"), '-forward'))
## methcols <- c(grey(0.6), hcl(rep(1:3, each=2)/3*360, l=c(50,70), c=c(70,50)))
## methcols <- methcols[-length(methcols)]
names(methcols) <- methods
## methods <- gsub('-forward', '', methods)
initmonths <- c("05", "11")
scores <- c("EnsRmsess", "EnsRmse", "EnsMe", "Ens2AFC", "FairCrpss", 'FairSprErr', 'FairCrps', 'EnsCorr')

skill.long  <- read_scores(model='ecmwf-system4',
                           index=indices,
                           score=scores,
                           method=submethods,
                           init=initmonths,
                           ccr=c(FALSE,TRUE),
                           grid='global2',
                           lead=2,
                           detrend=FALSE,
                           cleanup=TRUE)

## remove Heating and Cooling degree days in non-heating, cooling seasons
dind <- which((skill.long$index == 'HDDch' & skill.long$initmon == "05") | 
                (skill.long$index == "CDD" & skill.long$initmon == "11"))
if (length(dind) > 0) skill.long <- skill.long[-dind, ]

skill <- skill.long %>% filter(lsm > 0.5 & lat > -60) %>%
  mutate(lon = ifelse(lon > 188, lon - 360, lon),
         index = gsub("NDD", "DD", as.character(index)),
         index = gsub("CDD", "HDD/CDD", as.character(index)),
         index = gsub("HDDch", "HDD/CDD", as.character(index)),
         index = ifelse(substr(index, 1, 3) == "PDD", paste0('DD', formatC(100 - as.numeric(substr(index, 4,5)), width=2, flag=0)), index),
         gridarea = 4*cos(lat/180*pi))
skill$index <- factor(as.character(skill$index), unique(as.character(skill$index)))


ylims <- list(EnsCorr=c(-1,1), Ens2AFC=c(0,1), EnsCrps=c(0,3), EnsRmse=c(0.005,9), 
              FairCrpss=c(-1.5,1), FairSprErr=c(0,2.5), EnsMe=c(-2,2))


## plot all methods (non-recal) for the various scores
par(mfrow=c(3,2), mar=c(0.5, 0.5, 3, 0.5), oma=c(2,4,0,0), 
    tcl=0.3, mgp=c(1.5, 0.5, 0), cex.axis=1, cex.lab=1.4)
plot_boxplot(data=filter(skill, !ccr, method != 'conditionalRecal-forward', 
                   score %in% c("EnsCorr", "EnsRmse", "FairCrpss")))



for (thisscore in scores){
  png(paste0('/users/bhendj/doc/figures/recalibration/recalibration-', thisscore, '.png'), 
      width=8, height=5, units='in', res=200)
  p <- ggplot(subset(skill, !ccr & score == thisscore), aes(x=index, y=value)) + 
    facet_grid(. ~ initmon, labeller=relabeller) + 
    geom_hline(yintercept=ifelse(thisscore %in% c("EnsRmsess", "FairCrpss", "EnsCorr"), 
                                 0, ifelse(thisscore %in% "Ens2AFC", 0.5, 0)), lty=2) + 
    geom_boxplot(aes(fill=method), outlier.size=0.6, labeller=relabeller) + 
    ggtitle(paste0(thisscore, " (land only)")) + 
    ylab(thisscore)
  if (thisscore == "EnsRmse"){
    p <- p + coord_trans(ytrans="log10", limy=ylims[[thisscore]])
  }  else{
    p <- p + coord_cartesian(ylim=ylims[[thisscore]])
  }
  p <- p + 
    scale_fill_manual(name="Calibration methods",
                      breaks=methods, values=methcols,
                      labels=relabeller('method', methods)) + 
    theme_bw()  
  print(p)
  dev.off()  
}


## plot different scores for some of the methods
submethods <- c('none', 'smooth-forward', 'comb-forward', 'fastqqmap-forward')
subscores <- c('EnsCorr',  'FairCrpss', "FairSprErr")
skill$method <- factor(as.character(skill$method), submethods)
png(paste0('/users/bhendj/doc/figures/recalibration/recalibration-submethods_subscores.png'), 
    width=12, height=8, units='in', res=200)
print({
  ggplot(filter(skill, !ccr, score %in% subscores, method %in% submethods) %>%
           mutate(score = factor(as.character(score), subscores)) %>% 
           mutate(value2 = ifelse(score == 'FairSprErr', log(value), value)), 
                  aes(x=index, y=value2)) + 
    facet_grid(initmon ~ score, labeller=relabeller) + 
    geom_hline(yintercept=0, lty=2) + 
    geom_boxplot(aes(fill=method), outlier.size=0.6, labeller=relabeller) + 
    ylab("Verification metric (land only)") + 
    coord_cartesian(ylim=c(-1.5,1))+
    scale_fill_manual(name="Calibration methods",
                      breaks=methods, values=methcols,
                      labels=relabeller('method', methods)) + 
    theme_bw()  
})
dev.off()  



par(mfrow=c(length(submethods), 2), mar=c(0.5, 0.5, 2, 0.5))
for (smeth in c("none", "smooth-forward", "smoothRecal-forward", "fastqqmap-forward")){
  for (init in c("05", "11")){
    plot_scores(filter(skill, score == 'FairCrpss', 
                       !ccr, method==smeth, 
                       initmon==init, index=='DD02'),
                ylim=c(-60,85),
                lev=c(-1e4, seq(-0.7, 0.7, 0.2), 1))
    axis(3, at=par("usr")[1], paste(init, smeth), hadj=0, tick=F, line=-0.5)
  }
}






## take differences in skill metrics with respect to different references
noland.cast <- dcast(noland, ... ~ method, value.var='value')
none.diff <- subset(noland.cast[,!names(noland.cast) %in% c('none')], score %in% c('EnsRmse', 'FairCrps'))
none.diff[,which(names(none.diff) %in% methods[-1])] <- none.diff[,which(names(none.diff) %in% methods[-1])]/subset(noland.cast, score %in% c('EnsRmse', 'FairCrps'))[,'none']
none.diff.long <- melt(none.diff, id.vars=c('index', 'initmon', 'score', 'gridID', 'area', 'ccr', 'lsm', 'lon', 'lat'), variable.name='method')




png('/users/bhendj/doc/figures/recalibration/recalibration-CRPSS.png', width=8, height=5, units='in', res=200)
print({
  ggplot(subset(noland, !ccr & score == 'FairCrpss'), aes(x=index, y=-tanh(log(1 - value)))) + 
    facet_grid(. ~ initmon, labeller=relabeller) + 
    geom_hline(yintercept=0) +
    geom_boxplot(aes(fill=method), outlier.size=0.6, labeller=relabeller) + 
    ggtitle("CRPSS (land only)") + 
    ylab("CRPSS'") + 
    scale_fill_manual(name="Calibration methods",
                      breaks=methods, values=methcols,
                      labels=relabeller('method', methods)) + 
    theme_bw()  
})
dev.off()

improvement <- subset(noland.cast, score %in% c('FairCrps', 'EnsRmse'))
improvement[,methods] <- 1 - improvement[,methods] / improvement[,'smooth-forward'] 
improvement.long <- melt(improvement, id.vars=1:9, variable.name="method")

## improvement (CRPSS) over debiased forecasts
png('/users/bhendj/doc/figures/recalibration/recalibration-comparison_of_calibration_methods.png', width=8, height=5, units='in', res=200)
print({
  ggplot(subset(improvement.long, !ccr & score == 'FairCrps' & method %in% methods[-grep('Recal', methods)]), aes(x=index, y=-tanh(log(1 - value)))) + 
    facet_grid(. ~ initmon, labeller=relabeller) + 
    geom_hline(yintercept=0) +
    geom_boxplot(aes(fill=method), outlier.size=0.6, labeller=relabeller) + 
    ggtitle("CRPSS - debias as reference forecast (land only)") + 
    ylab("CRPSS'") + 
    scale_fill_manual(name="Calibration methods",
                      breaks=methods, values=methcols,
                      labels=relabeller('method', methods)) + 
    theme_bw()  
})
dev.off()
  
improvement2 <- subset(noland.cast, score %in% c('FairCrps', 'EnsRmse'))[, -grep('none', names(noland.cast))]
improvement2[,methods[seq(2,8,2)]] <- 1 - improvement2[,methods[seq(3,9,2)]] / improvement2[,methods[seq(2,8,2)]]
improvement2.long <- melt(improvement2, id.vars=1:9, variable.name="method")

## improvement (CRPSS) of recalibration forecasts
png('/users/bhendj/doc/figures/recalibration/recalibration-effect_of_daily_recalibration.png', width=8, height=5, units='in', res=200)
print({
  ggplot(subset(improvement2.long, !ccr & score == 'FairCrps' & method %in% methods[seq(2,8,2)]), aes(x=index, y=-tanh(log(1 - value)))) + 
    facet_grid(. ~ initmon, labeller=relabeller) + 
    geom_hline(yintercept=0) +
    geom_boxplot(aes(fill=method), outlier.size=0.6, labeller=relabeller) + 
    ggtitle("CRPSS - no rescaling as reference forecast (land only)") + 
    ylab("CRPSS'") + 
    scale_fill_manual(name="Calibration methods",
                      breaks=methods, values=methcols,
                      labels=relabeller('method', methods)) + 
    theme_bw()  
})
dev.off()
  

recali <- dcast(subset(noland, score %in% c('FairCrps', 'EnsRmse')), ... ~ ccr, value.var='value')
recali$value <- 1 - recali[['TRUE']] / recali[['FALSE']]


png('/users/bhendj/doc/figures/map_recalibration-effect_of_recalibration_on_index_CRPSS_vs_combRecal.png', width=12, height=11, units='in', res=200)
print(ggplot(subset(recali, score == 'FairCrps' & method == 'combRecal-forward'), aes(x=lon, y=lat)) + 
        facet_grid(index ~ initmon, labeller = relabeller) + 
        geom_tile(aes(fill=-tanh(log(1 - value)))) + 
        coord_equal() + 
        scale_fill_gradientn(name="CRPSS'", colours=rev(brewer.pal(11, 'RdBu')), limits=c(-1,1)) + 
        geom_path(data=mm, aes(x=x, y=y), colour=grey(0.4), size=0.5) + 
        theme_bw()
); dev.off()



## improvement (CRPSS) of recalibration forecasts
png('/users/bhendj/doc/figures/recalibration/recalibration-effect_of_recalibration_on_index_CRPSS.png', width=8, height=5, units='in', res=200)
print({
  ggplot(subset(recali, score == 'FairCrps'), aes(x=index, y=-tanh(log(1 - value)))) + 
    facet_grid(. ~ initmon, labeller=relabeller) + 
    geom_hline(yintercept=0) +
    geom_boxplot(aes(fill=method), outlier.size=0.6, labeller=relabeller) + 
    ggtitle("CRPSS - no recalibration on index as reference (land only)") + 
    ylab("CRPSS'") + 
    scale_fill_manual(name="Calibration methods",
                      breaks=methods, values=methcols,
                      labels=relabeller('method', methods)) + 
    theme_bw()  
})
dev.off()

## improvement (Spread to error ratio) of recalibration forecasts
png('/users/bhendj/doc/figures/recalibration/recalibration-effect_of_recalibration_on_index_SPRERR.png', width=8, height=8, units='in', res=200)
print({
  ggplot(subset(noland, score == 'FairSprErr'), aes(x=index, y=value)) + 
    facet_grid(ccr ~ initmon, labeller=relabeller) + 
    geom_hline(yintercept=1) +
    geom_boxplot(aes(fill=method), outlier.size=0.6, labeller=relabeller) + 
    ggtitle("Spread to error ratio (land only)") + 
    ylab("Spread to error ratio") + 
    scale_fill_manual(name="Calibration methods",
                      breaks=methods, values=methcols,
                      labels=relabeller('method', methods)) + 
    scale_y_continuous(trans=log2_trans(), 
                       limits=c(0.1, 10)) + 
    theme_bw()  
})
dev.off()





















### plot global mean CRPS of HDD/CDD with and without bias correction


skill.cast <- dcast(skill.long, ... ~ method, value.var='value')

for (index in c('HDD/CDD', 'DD10', 'DD02')){
  indname <- gsub("HDD/CDD", "HDDCDD", index)
  
  ski <- subset(skill.long[skill.long$index == index, ], score == "FairCrpss" & method %in% c('none', 'smooth-forward', 'comb-forward') & ! ccr  & lsm > 0.5 & lat > -60)
  ski2 <- subset(skill.cast[skill.cast$index == index, ], score == 'FairCrps' & !ccr & lsm > 0.5 & lat > -60)
  ski2$improvement <- 1 - ski2[[methods[2]]] / ski2$none
  
  ski3 <- merge(ski[,-which(names(ski) == 'score')], ski2[,-which(names(ski2) == 'score')])
  ski3$value[ski3$method == methods[2]] <- ski3$improvement[ski3$method == methods[2]]
  
  
  png(paste0('/users/bhendj/doc/figures/recalibration/recalibration-', indname, '_CRPSS_with_and_without_calibration.png'), width=12, height=5, units='in', res=200)
  print(ggplot(ski, aes(x=lon, y=lat)) + 
          facet_grid(method ~ initmon, labeller = relabeller) + 
          geom_tile(aes(fill=-tanh(log(1 - value)))) + 
          coord_equal() + 
          scale_fill_gradientn(name="CRPSS'", 
                               colours=rev(brewer.pal(11, 'RdBu')), 
                               limits=c(-1,1)) + 
          geom_path(data=mm, aes(x=x, y=y), colour=grey(0.4), size=0.5) + 
          theme_bw()
  ); dev.off()
  
  png(paste0('/users/bhendj/doc/figures/recalibration/recalibration-', indname,'_CRPSS_and_improvement_with_and_without_calibration.png'), width=12, height=5, units='in', res=200)
  print(ggplot(ski3, aes(x=lon, y=lat)) + 
          facet_grid(method ~ initmon, labeller = relabeller) + 
          geom_tile(aes(fill=-tanh(log(1 - value)))) + 
          coord_equal() + 
          scale_fill_gradientn(name="CRPSS'", colours=rev(brewer.pal(11, 'RdBu')), limits=c(-1,1)) + 
          geom_path(data=mm, aes(x=x, y=y), colour=grey(0.4), size=0.5) + 
          theme_bw()
  ); dev.off()
  
  
  ski <- subset(skill.long[skill.long$index == index, ], score == "EnsCorr" & method %in% c('none', methods[2]) & ! ccr  & lsm > 0.5 & lat > -60)
  
  png(paste0('/users/bhendj/doc/figures/recalibration/recalibration-', indname, '_corr_with_and_without_calibration.png'), width=12, height=5, units='in', res=200)
  print(ggplot(ski, aes(x=lon, y=lat)) + 
          facet_grid(method ~ initmon, labeller = relabeller) + 
          geom_tile(aes(fill=value)) + 
          coord_equal() + 
          scale_fill_gradientn(name="correlation", colours=rev(brewer.pal(11, 'RdBu')), limits=c(-1,1)) + 
          geom_path(data=mm, aes(x=x, y=y), colour=grey(0.4), size=0.5) + 
          theme_bw()
  ); dev.off()
  
  ski <- subset(skill.long[skill.long$index == index, ], score == "Ens2AFC" & method %in% c('none', methods[2]) & ! ccr  & lsm > 0.5 & lat > -60)
  
  png(paste0('/users/bhendj/doc/figures/recalibration/recalibration-', indname, '_2AFC_with_and_without_calibration.png'), width=12, height=5, units='in', res=200)
  print(ggplot(ski, aes(x=lon, y=lat)) + 
          facet_grid(method ~ initmon, labeller = relabeller) + 
          geom_tile(aes(fill=value)) + 
          coord_equal() + 
          scale_fill_gradientn(name="2AFC", colours=rev(brewer.pal(11, 'RdBu')), limits=c(0,1)) + 
          geom_path(data=mm, aes(x=x, y=y), colour=grey(0.4), size=0.5) + 
          theme_bw()
  ); dev.off()
  
}  

## compare index to underlying variable
ski <- subset(skill.long, score == "FairCrpss" & index %in% c('tas', 'DD10') & method %in% c(methods[2]) & ! ccr  & lsm > 0.5 & lat > -60)

png('/users/bhendj/doc/figures/recalibration/recalibration-DD10_vs_tas_with_calibration.png', width=12, height=5, units='in', res=200)
print(ggplot(ski, aes(x=lon, y=lat)) + 
        facet_grid(index ~ initmon, labeller = relabeller) + 
        geom_tile(aes(fill=-tanh(log(1 - value)))) + 
        coord_equal() + 
        scale_fill_gradientn(name="CRPSS'", colours=rev(brewer.pal(11, 'RdBu')), limits=c(-1,1)) + 
        geom_path(data=mm, aes(x=x, y=y), colour=grey(0.4), size=0.5) + 
        theme_bw()
); dev.off()

ski <- subset(skill.long, score == "FairCrpss" & index %in% c('tas', 'HDD/CDD') & method %in% c(methods[2]) & ! ccr  & lsm > 0.5 & lat > -60)

png('/users/bhendj/doc/figures/recalibration/recalibration-HDD_vs_tas_with_calibration.png', width=12, height=5, units='in', res=200)
print(ggplot(ski, aes(x=lon, y=lat)) + 
        facet_grid(index ~ initmon, labeller = relabeller) + 
        geom_tile(aes(fill=-tanh(log(1 - value)))) + 
        coord_equal() + 
        scale_fill_gradientn(name="CRPSS'", colours=rev(brewer.pal(11, 'RdBu')), limits=c(-1,1)) + 
        geom_path(data=mm, aes(x=x, y=y), colour=grey(0.4), size=0.5) + 
        theme_bw()
); dev.off()

ski.cast <- dcast(ski, ... ~ index, value.var='value')
ski.cast$tas[is.na(ski.cast[['HDD/CDD']])] <- NA
ski <- melt(ski.cast, id.vars=1:9, variable.name='index')

png('/users/bhendj/doc/figures/recalibration/recalibration-HDD_vs_tas_with_calibration_masked.png', width=12, height=5, units='in', res=200)
print(ggplot(ski, aes(x=lon, y=lat)) + 
        facet_grid(index ~ initmon, labeller = relabeller) + 
        geom_tile(aes(fill=-tanh(log(1 - value)))) + 
        coord_equal() + 
        scale_fill_gradientn(name="CRPSS'", colours=rev(brewer.pal(11, 'RdBu')), limits=c(-1,1)) + 
        geom_path(data=mm, aes(x=x, y=y), colour=grey(0.4), size=0.5) + 
        theme_bw()
); dev.off()
# 
# for (index in c('tas', 'DD10', 'DD02', 'HDD/CDD')){
#   
#   print(index)
#   indname <- gsub('HDD/CDD', 'HDDCDD', index)
#   ## effect of recalibration
#   ski <- subset(skill.long[skill.long$index == index,], score == "FairCrpss" & method %in% c(methods[2]) & lsm > 0.5 & lat > -60)
#   
#   png(paste0('/users/bhendj/doc/figures/recalibration/recalibration-', indname, '_calibration_vs_recalibration.png'), width=12, height=5, units='in', res=200)
#   print(ggplot(ski, aes(x=lon, y=lat)) + 
#           facet_grid(ccr ~ initmon, labeller = relabeller) + 
#           geom_tile(aes(fill=-tanh(log(1 - value)))) + 
#           coord_equal() + 
#           scale_fill_gradientn(name="CRPSS'", colours=rev(brewer.pal(11, 'RdBu')), limits=c(-1,1)) + 
#           geom_path(data=mm, aes(x=x, y=y), colour=grey(0.4), size=0.5) + 
#           theme_bw()
#   ); dev.off()
#   
#   ski <- subset(skill.long[skill.long$index == index,], score == "FairCrpss" & method %in% c(methods[2]) & lsm > 0.5 & lat > -60)
#   ski2 <- dcast(subset(skill.long[skill.long$index == index, ], score == 'FairCrps' &  lsm > 0.5 & lat > -60), ... ~ ccr, value.var = 'value')
#   ski2$improvement <- 1 - ski2[['TRUE']] / ski2[['FALSE']]
#   
#   ski3 <- merge(ski[,-which(names(ski) == 'score')], ski2[,-which(names(ski2) == 'score')])
#   ski3$value[ski3$ccr] <- ski3$improvement[ski3$ccr]
#   
#   png(paste0('/users/bhendj/doc/figures/recalibration/recalibration-', indname, '_calibration_vs_recalibration_improvement.png'), width=12, height=5, units='in', res=200)
#   print(ggplot(ski3, aes(x=lon, y=lat)) + 
#           facet_grid(ccr ~ initmon, labeller = relabeller) + 
#           geom_tile(aes(fill=-tanh(log(1 - value)))) + 
#           coord_equal() + 
#           scale_fill_gradientn(name="CRPSS'", colours=rev(brewer.pal(11, 'RdBu')), limits=c(-1,1)) + 
#           geom_path(data=mm, aes(x=x, y=y), colour=grey(0.4), size=0.5) + 
#           theme_bw()
#   ); dev.off()
#   
#   if (index == 'tas'){
#     
#     ski <- subset(skill.long[skill.long$index == index, ], score == "FairCrpss" & method %in% c('smoothccr') & lsm > 0.5 & lat > -60)
#     
#     png(paste0('/users/bhendj/doc/figures/recalibration/recalibration-', indname, '_smoothccr_vs_recalibration.png'), width=12, height=5, units='in', res=200)
#     print(ggplot(ski, aes(x=lon, y=lat)) + 
#             facet_grid(ccr ~ initmon, labeller = relabeller) + 
#             geom_tile(aes(fill=-tanh(log(1 - value)))) + 
#             coord_equal() + 
#             scale_fill_gradientn(name="CRPSS'", colours=rev(brewer.pal(11, 'RdBu')), limits=c(-1,1)) + 
#             geom_path(data=mm, aes(x=x, y=y), colour=grey(0.4), size=0.5) + 
#             theme_bw()
#     ); dev.off()
#     
#     ski <- subset(skill.long[skill.long$index == index,], score == "FairSprErr" & method %in% c('smoothccr') & lsm > 0.5 & lat > -60)
#     
#     png(paste0('/users/bhendj/doc/figures/recalibration/recalibration-SprErr_', indname, '_smoothccr_vs_recalibration.png'), width=12, height=5, units='in', res=200)
#     print(ggplot(ski, aes(x=lon, y=lat)) + 
#             facet_grid(ccr ~ initmon, labeller = relabeller) + 
#             geom_tile(aes(fill=log(value))) + 
#             coord_equal() + 
#             scale_fill_gradientn(name="log(SprErr)", colours=rev(brewer.pal(11, 'RdBu')), limits=c(-1.5, 1.5)) + 
#             geom_path(data=mm, aes(x=x, y=y), colour=grey(0.4), size=0.5) + 
#             theme_bw()
#     ); dev.off()
#     
#   }  
# }
# 
# 
# ## Distribution of skill scores for various indices and calibration methods
# 
# 
# png('/users/bhendj/doc/figures/recalibration/recalibration-boxplot_CRPS.png', width=8, height=5, units='in', res=200)
# print({
#   ggplot(subset(skill.long, score %in% c('FairCrps') & lsm > 0.5 & lat > -60 & !ccr), 
#          aes(x=index, y=value)) + 
#     geom_boxplot(aes(fill=method), outlier.size=0.6) + 
#     facet_grid(. ~ initmon, labeller = relabeller) + 
#     ggtitle("Mean CRPS (land only)") + 
#     scale_y_continuous(trans=log2_trans(), 
#                        breaks=c(0.001, 0.01, 0.1, 1, 10), 
#                        limits=c(0.001, 10)) + 
#     ylab("") + 
#     scale_fill_manual(name="Calibration methods",
#                       breaks=methods, values=methcols,
#                       labels=relabeller('method', methods)) + 
#     theme_bw()  
# })
# dev.off()
# 
# 
# 
# 
# 
# png('/users/bhendj/doc/figures/recalibration/recalibration-boxplot_CRPSS.png', width=8, height=5, units='in', res=200)
# print({  
#   ggplot(subset(skill.long, score %in% c('FairCrpss') & lsm > 0.5 & lat > -60 & !ccr), 
#          aes(x=index, y=-tanh(log(1 - value)))) + 
#     geom_hline(yintercept=0) +
#     geom_boxplot(aes(fill=method), outlier.size=0.6) + 
#     facet_grid(. ~ initmon, labeller = relabeller) + 
#     ggtitle("CRPSS' (land only)") + 
#     scale_y_continuous(limits=c(-1, 1)) + 
#     ylab("") + 
#     scale_fill_manual(name="Calibration methods",
#                       breaks=methods, values=methcols,
#                       labels=relabeller('method', methods)) +
#     theme_bw()
# })
# dev.off()
# 
# 
# 
# 
# 
# png('/users/bhendj/doc/figures/recalibration/recalibration-boxplot_CRPSS_improvement.png', width=8, height=5, units='in', res=200)
# print({  
#   ggplot(subset(none.diff.long, score %in% c('FairCrps') & lsm > 0.5 & lat > -60 & !ccr), 
#          aes(x=index, y= -tanh(log(value)))) + 
#     geom_hline(yintercept=0) + 
#     geom_boxplot(aes(fill=method), outlier.size=0.6) + 
#     facet_grid(. ~ initmon, labeller = relabeller) + 
#     ggtitle("CRPSS' wrt to no calibration (land only)") +
#     scale_y_continuous(limits=c(-1, 1)) + 
#     ylab("") + 
#     scale_fill_manual(name="Calibration methods",
#                       breaks=methods, values=methcols,
#                       labels=relabeller('method', methods)) +
#     theme_bw()
# })
# dev.off()
# 
# 
# 
# 
# 
# png('/users/bhendj/doc/figures/recalibration/recalibration-boxplot_correlation.png', width=8, height=5, units='in', res=200)
# print({  
#   ggplot(subset(skill.long, score %in% c('EnsCorr') & lsm > 0.5 & lat > -60 & ! ccr), 
#          aes(x=index, y=value)) + 
#     geom_boxplot(aes(fill=method), outlier.size=0.6) + 
#     facet_grid(. ~ initmon, labeller = relabeller) + 
#     ggtitle("Correlation with the ensemble mean (land only)") + 
#     ylim(c(-1,1)) + 
#     ylab("") + 
#     scale_fill_manual(name="Calibration methods",
#                       breaks=methods, values=methcols,
#                       labels=relabeller('method', methods)) +
#     theme_bw()
# })
# dev.off()
# 
# 
# 
# 
# 
# 
# png('/users/bhendj/doc/figures/recalibration/recalibration-boxplot_spread_to_error.png', width=8, height=5, units='in', res=200)
# print({  
#   ggplot(subset(skill.long, score %in% c('FairSprErr') & lsm > 0.5 & lat > -60 & !ccr), 
#          aes(x=index, y=value)) + 
#     geom_hline(aes(yintercept=1)) + 
#     geom_boxplot(aes(fill=method), outlier.size=0.6) + 
#     facet_grid(. ~ initmon, labeller = relabeller) + 
#     ggtitle("Spread to error ratio (land only)") + 
#     scale_y_continuous(trans=log2_trans(), 
#                        breaks=c(0.1, 0.33, 1, 3), 
#                        limits=c(0.1, 5)) + 
#     ylab("") + 
#     scale_fill_manual(name="Calibration methods",
#                       breaks=methods, values=methcols,
#                       labels=relabeller('method', methods)) +
#     theme_bw()
# })
# dev.off()
# 
# 
# 
# 
# png('/users/bhendj/doc/figures/recalibration/recalibration-boxplot_RMSE.png', width=8, height=5, units='in', res=200)
# print({  
#   ggplot(subset(skill.long, score %in% c('EnsRmse') & lsm > 0.5 & lat > -60 & !ccr), 
#          aes(x=index, y=value)) + 
#     geom_boxplot(aes(fill=method), outlier.size=0.6) + 
#     facet_grid(. ~ initmon, labeller = relabeller) + 
#     ggtitle("Root mean square error of ensemble mean (land only)") + 
#     scale_y_continuous(trans=log2_trans(), 
#                        breaks=c(0.033, 0.1, 0.33, 1, 3, 10), 
#                        limits=c(0.033, 10)) + 
#     ylab("") + 
#     scale_fill_manual(name="Calibration methods",
#                       breaks=methods, values=methcols,
#                       labels=relabeller('method', methods)) +
#     theme_bw()
# })
# dev.off()
# 
# 
# 
# 
# 
# png('/users/bhendj/doc/figures/recalibration/recalibration-boxplot_RMSESS_improvement.png', width=8, height=5, units='in', res=200)
# print({  
#   ggplot(subset(none.diff.long, score %in% c('EnsRmse') & lsm > 0.5 & lat > -60 & !ccr), 
#          aes(x=index, y=-tanh(log(value)))) + 
#     geom_boxplot(aes(fill=method), outlier.size=0.6) + 
#     facet_grid(. ~ initmon, labeller = relabeller) + 
#     ggtitle("RMSESS' wrt no calibration (land only)") +
#     scale_y_continuous(limits=c(-1, 1)) + 
#     ylab("") + 
#     scale_fill_manual(name="Calibration methods",
#                       breaks=methods, values=methcols,
#                       labels=relabeller('method', methods)) +
#     theme_bw()
# })
# dev.off()
# 
# 
# 
# 
# 
# png('/users/bhendj/doc/figures/recalibration/recalibration-boxplot_bias.png', width=8, height=5, units='in', res=200)
# print({  
#   ggplot(subset(skill.long, score %in% c('EnsMe') & lsm > 0.5 & lat > -60 & !ccr), 
#          aes(x=index, y=value)) + 
#     geom_boxplot(aes(fill=method), outlier.size=0.6) + 
#     facet_grid(. ~ initmon, labeller = relabeller) + 
#     ggtitle("Mean error (land only)") +
#     ylim(c(-2,2)) + 
#     ylab("") + 
#     scale_fill_manual(name="Calibration methods",
#                       breaks=methods, values=methcols,
#                       labels=relabeller('method', methods))  +
#     theme_bw()
# })
# dev.off()
# 

## plot for klimaseminar
scores <- c('FairSprErr', 'EnsCorr', 'FairCrpss')
submethods <- c("smooth-forward", "smoothRecal-forward")

skill.long  <- read_scores(model='ecmwf-system4',
                           index=indices,
                           score=scores,
                           method=submethods,
                           init=initmonths,
                           ccr=c(FALSE,TRUE),
                           grid='global2',
                           lead=2,
                           detrend=FALSE,
                           cleanup=TRUE)

skill <- skill.long %>% filter(lsm > 0.5 & lat > -60) %>%
  mutate(lon = ifelse(lon > 188, lon - 360, lon),
         index = gsub("NDD", "DD", as.character(index)),
         index = gsub("CDD", "HDD/CDD", as.character(index)),
         index = gsub("HDDch", "HDD/CDD", as.character(index)),
         index = ifelse(substr(index, 1, 3) == "PDD", paste0('DD', formatC(100 - as.numeric(substr(index, 4,5)), width=2, flag=0)), index),
         gridarea = 4*cos(lat/180*pi)) %>%
  mutate(method2 = factor(paste0(c('debias', 'debiasRecal')[as.numeric(method)], ifelse(ccr, '_CCR', '')), c('debias', 'debias_CCR', 'debiasRecal', 'debiasRecal_CCR'))) %>%
  filter(method2 %in% c("debias_CCR", "debiasRecal"))
skill$index <- factor(as.character(skill$index), unique(as.character(skill$index)))


for (score in scores){
  
png(paste0('/users/bhendj/doc/figures/klimaseminar-recalibration-',score,'.png'), width=8, height=5, units='in', res=200)
print({
  ggplot(skill[skill$score == score,], 
         if (score == 'FairSprErr') { 
           aes(x=index, y=log(value)) 
         } else {
           aes(x=index, y=value)
         }) + 
    facet_grid(. ~ initmon, labeller=relabeller) + 
    geom_hline(yintercept=0) +
    geom_boxplot(aes(fill=method2), outlier.size=0.6, labeller=relabeller) + 
    coord_cartesian(ylim=if (score == 'FairSprErr'){
      c(-2,1)
    } else if (score == 'EnsCorr'){
      c(-1,1)
    } else {
      c(-2,1)
    }) + 
    ggtitle("") + 
    ylab(relabeller('score', score))+ 
    theme_bw()  
})
dev.off()
}


