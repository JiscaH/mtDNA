# AR, ER & runtime versus proportion of parents genotyped
# nSnp=200, Err=5e-3, Ped_HSg5 + deer ped
# sequoia versions [as in 2017 MS], 2.5.6, 2.7, 2.7+mtDNA



#===============================================================================
# ==== sim & run ====
#===============================================================================

setwd('D:/Work/mtDNA/ARER_2023-08-09')

PM <- rev(1 - c(seq(0, 0.25, 0.05), seq(0.4, 1, 0.2)))
LL <- c(400, 200)
EE <- c(1e-4, 5e-3)
PEDZ <- c('HSg5', 'swan')
pkg_v <- packageVersion('sequoia')

NSIM <- 10

z <- 1
#for (z in 2:1) {
  if (PEDZ[z]=='HSg5') {
    Ped <- sequoia::Ped_HSg5
    LH <- sequoia::LH_HSg5
    mtM <- as.matrix(read.table('../PedHSg5_mtM.txt', header=TRUE, row.names=1))
  } else if (PEDZ[z]=='swan') {
    Ped <- read.table('../test_swan/Ped_swan.txt', header=TRUE)
    LH <- read.table('../test_swan/LH_swan.txt', header=TRUE)
    mtM <- as.matrix(read.table('../test_swan/swan_mtM.txt', header=TRUE, row.names=1))
  }

  for (y in 1:2) {  # high quality / medium quality SNP data
    for (m in 2:1) {  # without/with mtM
      if (m==2 & pkg_v == '2.5.6')  next
      for (x in seq_along(PM)) {
        if (PEDZ[z]=='HSg5' & PM[x] >0.95 & m==1)  next   # TODO: later
        FileName <- glue::glue('EC_{pkg_v}_{PEDZ[z]}_L{LL[y]}_E{EE[y]}_PM{PM[x]*10}')
        if (m==1) {
          ARGS_SEQ <- list(Module='ped', Err=EE[y], CalcLLR=FALSE)
          FileName <- paste0(FileName, '.RDS')
        } else {
          ARGS_SEQ <- list(Module='ped', Err=EE[y], CalcLLR=FALSE, mtSame=mtM)
          FileName <- paste0(FileName, '_mt.RDS')
        }
        if (file.exists(FileName))  next
        cat('\n\n', format(Sys.time(), '%H:%M:%S'), '\t', PEDZ[z], '\t', LL[y], '\t', PM[x], '\t', m, '\n\n')

        EC <- sequoia::EstConf(Pedigree = Ped,
                               LifeHistData = LH,
                               args.sim = list(nSnp=LL[y], SnpError=EE[y], ParMis=PM[x]),
                               args.seq = ARGS_SEQ,
                               nSim=NSIM,
                               nCores=5)

        EC$pairwise <- list()
        for (i in 1:NSIM) {
          EC$pairwise[[i]] <- ComparePairs(Ped_HSg5, EC$Pedigree.inferred[[i]],
                                           GenBack=2, patmat=TRUE)
        }
        saveRDS(EC[c('ConfProb', 'PedComp.fwd', 'RunParams', 'RunTime', 'pairwise')], file = FileName)


        saveRDS(EC[c('ConfProb', 'PedComp.fwd', 'RunParams', 'RunTime')], file = FileName)
      }
    }
  }
#}



#===============================================================================
# ==== plot ====
#===============================================================================

setwd('D:/Work/mtDNA/ARER_2023-08-09')

# max no. correct assigned / max no. errors
Nped <- cbind(AR = c('FS'=2*1014, 'HSg5'=2*960, 'deer'=1642+1202,
                     'griffin'=167+163, 'swan'=2*399),
              ER = 2*c(1157, 1000, 1998, 200, 500))


PM <- rev(1 - c(seq(0, 0.25, 0.05), seq(0.4, 1, 0.2)))
LL <- c(400, 200)
EE <- c(1e-4, 5e-3)
PEDZ <- c('HSg5', 'swan')
PV <- c('2.5.6', '2.7.0.2', '2.7.0.2')

ARERT <- array(dim = c(length(PEDZ),2,length(PV),length(PM),3),
               dimnames = list(PEDZ, c('goodQ', 'mediumQ'),
                               c('2.5.6', '2.7.0.2', '2.7.0.2 +mt'),
                               PM, c('AR', 'ER', 'Time')))
for (v in seq_along(PV)) {
  for (z in seq_along(PEDZ)) {
    for (q in 1:2) {
      for (x in seq_along(PM)) {
        FileName <- glue::glue('EC_{PV[v]}_{PEDZ[z]}_L{LL[q]}_E{EE[q]}_PM{PM[x]*10}')
        FileName <- paste0(FileName, ifelse(v==3 | v==4, '_mt.RDS', '.RDS'))
        #       cat(v, z, q, x, FileName, file.exists(FileName), '\n')
        if (!file.exists(FileName))  next
        EC <- readRDS(FileName)
        ARERT[z,q,v,x,'AR'] <- 1 - mean(apply(EC$PedComp.fwd[,'TT','Match',], 1, sum) / Nped[PEDZ[z], "AR"], na.rm=TRUE)
        ARERT[z,q,v,x,'ER'] <- mean(apply(EC$PedComp.fwd[,'TT',c('Mismatch','P2only'),], 1, sum) / Nped[PEDZ[z], "ER"], na.rm=TRUE)
        ARERT[z,q,v,x,'Time'] <- mean(EC$RunTime)/60   # time in minutes
      }
    }
  }
}

# table(!is.na(ARERT['swan',,,,'AR']))
# apply(ARERT, c(1,5), range, na.rm=TRUE)


YLAB <- c("AR" = "False negative rate",
          "ER" = "Error rate",
          "Time" = "Runtime (min)")
Yzero <- 1e-5
COL <- setNames(hcl.colors(n=length(PV), palette="Zissou 1"), dimnames(ARERT)[[3]])
MaxTime <- c('swan' = 10, 'HSg5' = 60)


PZ <-  'HSg5' #  'swan'  #

#pdf('plot_ARER_parmis_swan_prelim2.pdf', width=7, height=8)
par(mfcol=c(3, 2), mai=c(.3, .6, 0, 0), omi=c(.5,.4,.7,.1), xpd=F)
for (q in 1:2) {
  for (AE in c("AR", "ER", "Time")) {
    if (AE == "AR")  Ylim <- c(Yzero, 0.6)
    if (AE == "ER")  Ylim <- c(Yzero, 0.035)
    if (AE =="Time") Ylim <- c(0, MaxTime[[PZ]])
    plot(1,1, type="n", xlim=c(0, 1), ylim=Ylim, ylab="", xlab="", las=1,
         cex.lab=1.2, cex.axis=1.2)  # , log = ifelse(AE=="Time", "", "y"))
    abline(h=axTicks(2), col="lightgrey")
    if (AE == "AR")  mtext(c('400 SNPs, Err=1e-4', '200 SNPs, Err=5e-3')[q],
                           side=3, line=1, cex=1.2)
    if (AE == 'Time')  mtext('Prop. non-genotyped parents', side=1, line=4, cex=1.1)
    if (q==1)  mtext(YLAB[[AE]], side=2, line=5, cex=1.1)

    for (v in seq_along(PV)) {
      Y <- ARERT[PZ,q,v,,AE]
      #      if (AE != "Time")  Y[Y==0] <- Yzero   # else not plotted on log scale
      lines(PM, Y, type="b", col=COL[v], lwd=2)
    }
    if (AE == "AR" & q==1) {
      legend("topleft", names(COL), lty=1, lwd=2, col=COL, inset=.02, title="Version")
    }
  }
}
#dev.off()



# ARERT['HSg5',,'2.5.6','1',]
#                AR      ER      Time
# goodQ   0.5109375 0.38405  57.03325
# mediumQ 0.4636979 0.41005 101.98348


#===============================================================================
# ==== GPs ====
#===============================================================================

# How do new edits affect number of assigned grandparents?

setwd('D:/Work/mtDNA/ARER_2023-08-01')

# max no. correct assigned / max no. errors
Nped <- cbind(AR = c('FS'=2*1014, 'HSg5'=2*960, 'deer'=1642+1202,
                     'griffin'=167+163, 'swan'=2*399),
              ER = 2*c(1157, 1000, 1998, 200, 500))


PM <- rev(1 - c(seq(0, 0.25, 0.05), seq(0.4, 1, 0.2)))
LL <- c(400, 200)
EE <- c(1e-4, 5e-3)
PEDZ <- c('HSg5', 'swan')
PV <- c('2.5.6', '2.7.0.2', '2.7.0.2')

ARER_D <- array(dim = c(length(PEDZ),2,length(PV),length(PM), 3, 2),
               dimnames = list(PEDZ, c('goodQ', 'mediumQ'),
                               c('2.5.6', '2.7.0.2', '2.7.0.2 +mt'),
                               PM, c('DG', 'DD', 'DT'), c('AR', 'ER')))
for (v in seq_along(PV)) {
  for (z in seq_along(PEDZ)) {
    for (q in 1:2) {
      for (x in seq_along(PM)) {
        FileName <- glue::glue('EC_{PV[v]}_{PEDZ[z]}_L{LL[q]}_E{EE[q]}_PM{PM[x]*10}')
        FileName <- paste0(FileName, ifelse(v==3 | v==4, '_mt.RDS', '.RDS'))
        #       cat(v, z, q, x, FileName, file.exists(FileName), '\n')
        if (!file.exists(FileName))  next
        EC <- readRDS(FileName)
        PC_sum <- apply(EC$PedComp.fwd, 2:3, sum)[c('DG', 'DD', 'DT'), ]
        ARER_D[z,q,v,x,,'AR'] <- 1 - PC_sum[,'Match'] / PC_sum[,'Total']
        ARER_D[z,q,v,x,,'ER'] <- apply(PC_sum[,c('Mismatch','P2only')],1,sum) / PC_sum[,'Total']
      }
    }
  }
}


YLAB <- c("AR" = "False negative rate",
          "ER" = "Error rate",
          "Time" = "Runtime (min)")
Yzero <- 1e-5
COL <- c('DG' = '#018571', 'DD'='#a6611a', 'DT'='black')  # green - brown - black
LTY <- setNames(c(3,2,1), dimnames(ARER_D)[[3]])


PZ <- 'swan'

pdf('plot_ARER_parmis_swan_D.pdf', width=7, height=8)
par(mfcol=c(3, 2), mai=c(.3, .6, 0, 0), omi=c(.5,.7,.7,.1), xpd=F)
for (q in 1:2) {
  for (y in c('DG', 'DD', 'DT')) {
    if (y %in% c('DG', 'DD')) {
      AE <- 'AR'
    } else {
      AE <- 'ER'
    }
    if (AE == "AR")  Ylim <- c(Yzero, 1)
    if (AE == "ER")  Ylim <- c(Yzero, 0.05)
    plot(1,1, type="n", xlim=c(0, 1), ylim=Ylim, ylab="", xlab="", las=1,
         cex.lab=1.2, cex.axis=1.2)  # , log = ifelse(AE=="Time", "", "y"))
    abline(h=axTicks(2), col="lightgrey")
    if (y == "DG")  mtext(c('400 SNPs, Err=1e-4', '200 SNPs, Err=5e-3')[q],
                           side=3, line=1, cex=1.2)
    if (q==1)  mtext(YLAB[[AE]], side=2, line=4, cex=1.1)
    if (q==1)  mtext(y, side=2, line=6, cex=1.1, xpd=NA)

    for (v in seq_along(PV)) {
      lines(PM, ARER_D[PZ,q,v,,y,AE], type="b", col=COL[y], lty=LTY[v], lwd=2)
    }
    if (y == "DG" & q==1) {
      legend("topleft", names(LTY), lty=LTY, lwd=2, inset=.02, title="Version")
    }
  }
}
dev.off()

