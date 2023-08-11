
MainDir <- 'D:/Work/mtDNA'
library(glue)
source(glue('{MainDir}/SimPed.R'))
source(glue('{MainDir}/infer_mt_haplo.R'))

# test simped & sim_mtSame functions
Ped_IN <- SimPed(N_dam=50, EPP=e, N_offspring = c('1'=0.2, '2'=0.3, '3'=0.3, '4'=0.1, '5'=0.1))
mtSame <- sim_mtSame(Ped_IN)
Geno <- SimGeno(Ped_IN, nSnp=400, ParMis=1.0, SnpError=1e-4)
SeqOUT <- sequoia(Geno, Err=1e-4, quiet=FALSE, mtSame=mtSame, LifeHistData=Ped_IN[,c(1,4,5)])



#===============================================================================
# ==== AR & ER ====
#===============================================================================
# estimate assignment & error rate
# iterate over parameter combinations, and store each result in a separate data file
# filename contains parameter combo, for quick checking if this combo has been run already,
#  so that you can easily interrupt and run the remainder later.

# function EstConf()  edited to be able to simulate pedigrees & mtSame
source(glue('{MainDir}/ConfProb_simped.R'))

# folder to store output in
# dir.create(glue('{MainDir}/ARER_EPP_2023-08-09'))
setwd(glue('{MainDir}/ARER_2023-08-09'))

EPP_vals <- seq(0, 0.5, 0.1)
pkg_v <- packageVersion('sequoia')

# Note: no automatic elegant solution for with/without BY + Complx='simp';
# edited args.seq + withLH + FileName by hand

for (m in 1:2) {  # without/with mtM
  if (m==2 & pkg_v == '2.5.6')  next
  for (e in EPP_vals) {
    FileName <- glue::glue('EC_{pkg_v}_EPP{e}')  # _wBY_simp
    FileName <- paste0(FileName, ifelse(m==1, '.RDS', '_mt.RDS'))
    if (file.exists(FileName))  next
    cat('\n\n', format(Sys.time(), '%H:%M:%S'), '\t', e, '\t', m, '\n\n')

    EC <- EstConf_simped(simped_fun = SimPed,
                         args.simped = list(N_dam=200, EPP=e,
                                            N_offspring = c('1'=0.2, '2'=0.3, '3'=0.3, '4'=0.1, '5'=0.1)),
                         withLH = FALSE,
                         mt = m==2,
                         args.sim = list(nSnp=400, SnpError=1e-4, ParMis=1.0),
                         args.seq = list(Module='ped', Err=1e-4, CalcLLR=FALSE),  # , Complex='simp'
                         nSim=10,
                         nCores=5)
    EC$pairwise <- list()
    for (i in 1:10) {
      EC$pairwise[[i]] <- ComparePairs(EC$Pedigree.reference[[i]], EC$Pedigree.inferred[[i]],
                                       GenBack=1, patmat=TRUE)
    }

    saveRDS(EC[c('ConfProb', 'PedComp.fwd', 'RunParams', 'RunTime', 'pairwise')], file = FileName)
  }
}




#===============================================================================
# ==== plot ====
#===============================================================================

setwd('D:/Work/mtDNA/ARER_2023-08-09')

EPP_vals <- seq(0, 0.5, 0.1)
PV <- c('2.5.6', '2.7.0', '2.7.0.2', '2.7.0.2')

# create array to store the data in to be plotted
ARERT <- array(dim = c(length(PV), 2, length(EPP_vals), 3),
               dimnames = list(c('2.5.6', '2.7.0.0','2.7.0.2', '2.7.0.2 +mt'),
                               c('noBY', 'withBY'), EPP_vals,
                               c('AR', 'ER', 'Time')))
for (v in seq_along(PV)) {
  for (b in 1:2) {
    for (e in seq_along(EPP_vals)) {
      FileName <- glue::glue('EC_{PV[v]}_EPP{EPP_vals[e]}')
      if (b==2)  FileName <- paste0(FileName, '_wBY_simp')
      FileName <- paste0(FileName, ifelse(v!=4, '.RDS', '_mt.RDS'))
      if (!file.exists(FileName))  next
      EC <- readRDS(FileName)
      ARERT[v,b,e,'AR'] <- 1 - mean(apply(EC$PedComp.fwd[,'GD','Match',], 1, sum) / apply(EC$PedComp.fwd[,'GD','Total',], 1, sum))
      ARERT[v,b,e,'ER'] <- mean(apply(EC$PedComp.fwd[,'TT',c('Mismatch','P2only'),], 1, sum) / apply(EC$PedComp.fwd[,'TT','Total',], 1, sum))
      ARERT[v,b,e,'Time'] <- mean(EC$RunTime)/60   # time in minutes
    }
  }
}


YLAB <- c("AR" = "False negative rate",
          "ER" = "Error rate",
          "Time" = "Runtime (min)")
COL <- setNames(hcl.colors(n=length(PV), palette="Zissou 1"), dimnames(ARERT)[[1]])

par(mfcol=c(3, 1), mai=c(.3, .6, 0, 0), omi=c(.5,.4,.9,.1), xpd=F)
for (q in 1:2) {
  for (AE in c("AR", "ER", "Time")) {
    if (AE =="Time") Ylim <- c(0, 10)
    if (AE == "AR")  Ylim <- c(0, .6)
    if (AE == "ER")  Ylim <- c(0, .12)
    plot(1,1, type="n", xlim=c(0, 0.5), ylim=Ylim, ylab=YLAB[[AE]], xlab="", las=1,
         cex.lab=1.2, cex.axis=1.2)
    abline(h=axTicks(2), col="lightgrey")
    if (AE == 'Time')  mtext('EPP', side=1, line=4, cex=1.1)

    for (v in seq_along(PV)) {
      for (b in 1:2) {
        lines(EPP_vals, ARERT[v,b,,AE], type="b", col=COL[v], lwd=2, lty=b)
      }
    }
    if (AE == "AR") {
      legend("topleft", names(COL), lty=1, lwd=2, col=COL, inset=.02, title="Version")
    }
  }
}



#===============================================================================
# ==== pairwise comparisons ====
#===============================================================================

setwd('D:/Work/mtDNA/ARER_2023-08-09')
library(sequoia)
library(dplyr)

EPP_vals <- seq(0, 0.5, 0.1)
PV <- c( '2.7.0.2', '2.7.0.2', '2.5.6')
NSIM <- 10

# create array to store data in
# first run sham ComparePairs() to get dims & dimnames of its results
PC_tmp <- sequoia::ComparePairs(sequoia::Ped_HSg5, sequoia::Ped_HSg5, GenBack=1, patmat=TRUE)
PCC <- array(dim = c(length(PV), 2, length(EPP_vals), NSIM, dim(PC_tmp)),
               dimnames = c(list(c('2.7.0.2', '2.7.0.2 +mt','2.5.6'),
                               c('noBY', 'withBY'), EPP_vals, 1:NSIM),
                               dimnames(PC_tmp)))
for (v in seq_along(PV)) {
  for (b in 1:2) {
    for (e in seq_along(EPP_vals)) {
      FileName <- glue::glue('EC_{PV[v]}_EPP{EPP_vals[e]}')
      if (b==2)  FileName <- paste0(FileName, '_wBY_simp')
      FileName <- paste0(FileName, ifelse(v!=2, '.RDS', '_mt.RDS'))
      if (!file.exists(FileName))  next
      EC <- readRDS(FileName)
      PCC[v,b,e,,,] <- plyr::laply(EC$pairwise, function(x) x)
    }
  }
}


# counts as proportion of Ped1
Ped1_Totals <- apply(PCC[,,,,,-8],MARGIN=c(1:5), sum, na.rm=TRUE)  # excl 'X'=not genotyped
PCP <- sweep(PCC, MARGIN=c(1:5), STATS=Ped1_Totals, FUN='/')
# mean across iterations
PCM <- apply(PCP, c(1:3,5,6), mean, na.rm=TRUE)

RR <- c('FS', 'MHS', 'PHS', 'U')
COL <- c('FS' = 'forestgreen', 'MHS'='firebrick2', 'PHS'='dodgerblue2', 'U'='darkgrey')
LTY <- c(2,1,3)

pdf('Plot_EPP_pairwise.pdf', width=8, height=9)
par(mfcol=c(4,2), mai=c(.4,.5,.2,.3), omi=c(.3,0,.5,1.5))
for (b in 1:2) {
  for (Rx in RR) {
    plot(1,1, type="n", xlim=c(0, 0.5), ylim=c(0,1), ylab='', xlab="", las=1,
         cex.lab=1.2, cex.axis=1.2, main=paste('Actual', Rx))
    abline(h=axTicks(2), col="lightgrey")
    if (Rx=='FS')  mtext(dimnames(PCC)[[2]][b], side=3, line=3, xpd=NA)
    if (Rx=='U')  mtext('EPP', side=1, line=3, xpd=NA)

    for (v in seq_along(PV)) {
      for (Ry in RR) {
        lines(EPP_vals, PCM[v,b,,Rx,Ry], type="b", lwd=2, lty=LTY[v],
              col=ifelse(v==3, adjustcolor(COL[Ry],alpha.f=0.5), COL[Ry]))
      }
    }
  }
}
legend(0.55,5.5, legend=names(COL), col=COL, lwd=2, title='Assigned as', xpd=NA, cex=1.5)
legend(0.55,4.2, legend=dimnames(PCC)[[1]], lty=LTY, lwd=2, title='Sequoia version', xpd=NA, cex=1.5)
dev.off()
