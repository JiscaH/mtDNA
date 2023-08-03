
MainDir <- 'D:/Work/mtDNA'
library(glue)

source(glue('{MainDir}/SimPed.R'))

source(glue('{MainDir}/infer_mt_haplo.R'))


# test simped function
for (e in seq(0, 0.5, 0.1)) {
  cat('\n\n', e, '\n')
  Pedigree <- SimPed(N_dam=50, EPP=e, N_offspring = c('1'=0.2, '2'=0.3, '3'=0.3, '4'=0.1, '5'=0.1))
  mtSame <- sim_mtSame(Pedigree)
  Geno <- SimGeno(Pedigree, nSnp=400, ParMis=1.0, SnpError=1e-4)
  SeqOUT <- sequoia(Geno, Err=1e-4, quiet=TRUE, mtSame=mtSame)
  PC <- PedCompare(Pedigree, SeqOUT$Pedigree, Plot=FALSE)
  print(PC$Counts['TT',,])
}



# estimate assignment & error rate
source(glue('{MainDir}/ConfProb_simped.R'))

# folder to store output in
dir.create(glue('{MainDir}/ARER_EPP_2023-08-03'))
setwd(glue('{MainDir}/ARER_EPP_2023-08-03'))

EPP_vals <- seq(0, 0.5, 0.1)
pkg_v <- packageVersion('sequoia')

for (m in 1:2) {  # without/with mtM
  if (m==2 & pkg_v != '2.7.0')  next
  for (e in EPP_vals) {
    FileName <- glue::glue('EC_{pkg_v}_EPP{e}')
    FileName <- paste0(FileName, ifelse(m==1, '.RDS', '_mt.RDS'))
    if (file.exists(FileName))  next
    cat('\n\n', format(Sys.time(), '%H:%M:%S'), '\t', e, '\t', m, '\n\n')

    EC <- EstConf_simped(args.simped = list(N_dam=200, EPP=e,
                                            N_offspring = c('1'=0.2, '2'=0.3, '3'=0.3, '4'=0.1, '5'=0.1)),
                         mt = m==2,
                         args.sim = list(nSnp=400, SnpError=1e-4, ParMis=1.0),
                         args.seq = list(Module='ped', Err=1e-4, CalcLLR=FALSE),
                         nSim=10,
                         nCores=5)
    saveRDS(EC[c('ConfProb', 'PedComp.fwd', 'RunParams', 'RunTime')], file = FileName)
  }
}




# plot
# TODO
