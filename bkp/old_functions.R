plot_boxplot <- function(data, group1="index", 
                         group2="method", colours=NULL, 
                         weight=FALSE, ...){
  scores <- unique(data$score)
  inits <- unique(data$initmon)
  if (is.null(colours)){
    methods <- unique(data$method)
    colours <- hcl(seq(methods) / length(methods)*360, l=50, c=70)
    names(colours) <- methods
    if (any(methods == 'none')) colours['none'] <- gray(0.6)    
  }
  
  ylimtmp <- data %>% group_by(score) %>% summarise(ymin=min(value, na.rm=T), ymax=max(value, na.rm = T))
  ylims <- as.matrix(ylimtmp[,-1])
  rownames(ylims) <- as.character(ylimtmp$score)
  if ("EnsCorr" %in% scores) ylims['EnsCorr',] <- c(-1,1)
  if ("FairCrpss" %in% scores) ylims['FairCrpss',1] <- -2
  
  par(mfrow=c(length(scores), length(inits)), mar=c(5,3,3,1))
  j <- 0
  for (tscore in scores){
    for (tinit in inits){
      j <- j + 1
      
      dsub <- data %>% filter(initmon == tinit, score == tscore) %>%
        droplevels
      group1s <- unique(dsub[[group1]])
      group2s <- unique(dsub[[group2]])
      xat <- outer(seq(along=group2s), (seq(along=group1s) - 1)*(length(group2s) + 2), '+')
      
      
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
      boxplot(as.formula(paste("value ~ ", group2, "*", group1)),
              axes=FALSE, data=dsub, at=xat, add=TRUE,
              las=1, lty=1, pch=16, cex=0.5, 
              col=colours, 
              weight=if (weight) cos(dsub$lat / 180*pi) else NULL)
      box()
      if (tscore == scores[length(scores)]){
        axis(1, at=colMeans(xat), relabeller(group1, group1s), tick=F, cex.axis=par("cex.lab"))
      }
      axis(3, at=par("usr")[1], paste0(letters[j], ') ', relabeller("score", tscore), ' of forecasts for ', c(`11`='DJF (Nov. init.)', `05`="JJA (May init.)")[tinit]), hadj=0, line=0, cex.axis=par('cex.lab'), tick=F)
      if (tinit == inits[1] & tscore == scores[1]){
        legend('bottom', relabeller(group2, group2s), 
               fill=colours, 
               ncol=2, cex=par("cex.axis"), inset=0.02, bg='white')
      }
    }
  }
}
