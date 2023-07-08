# Count number of SNPs that are different between all pairs of individuals
# 
# Assume genotype data is formatted as:
# a matrix with 1 row per individual, 1 column per SNP, 
# 0/1/(2) copies of reference allele, rownames = individual IDs

# # Example:
# # simulate (diploid) genotype data with some duplicates
# Geno_act <- sequoia::SimGeno(sequoia::Ped_griffin, nSnp=100, ParMis=0, CallRate=1, SnpError=0)
# Geno_sim <- sequoia::MkGenoErrors(Geno_act, SnpError=0.05, CallRate=0.9)
# Geno_dups <- sequoia::MkGenoErrors(Geno_act[sample.int(nrow(Geno_act), 50), ], SnpError=0.05, CallRate=0.9)
# rownames(Geno_dups) <- paste0(rownames(Geno_dups), '_dup')
# Geno_sim <- rbind(Geno_sim, Geno_dups)
# 
# # count differences & plot
# difs <- Count_difs(Geno_sim)
# hist(difs['Different',,], breaks=c(0:100))
# plot(jitter(difs['NonMissing',,]), difs['Different',,], pch=16, col=adjustcolor(1, alpha.f=0.5))


Count_difs <- function(Geno,
                       Missing_codes = c(NA, -9, -1),  # allow different codings for missing
                       fill_matrix = FALSE)   # fill lower triangle & diagonal
{
  nInd <- nrow(Geno)
  IDs <- rownames(Geno)
  
  # array to store results in 
  OUT <- array(NA, dim=c(2,nInd,nInd),
               dimnames = list(c('NonMissing', 'Different'), IDs, IDs))
  # NonMissing: number of SNPs at which both individuals are successfully genotyped
  
  # only loop over upper triangle of nInd x nInd matrix
  for (i in 1:(nInd-1)) {
    for (j in (i+1):nInd) {
      
      # logical vector: both i + j succesfully genotyped
      snpd_both <- !Geno[i,] %in% Missing_codes & !Geno[j,] %in% Missing_codes
      OUT['NonMissing',i,j] <- sum(snpd_both)
      OUT['Different',i,j] <- sum(Geno[i,snpd_both] != Geno[j,snpd_both])
      
      # copy values to lower triangle
      if (fill_matrix) {
        OUT['NonMissing',j,i] <- OUT['NonMissing',i,j]
        OUT['Different',j,i] <- OUT['Different',i,j]
      }
      
    }
  }
  
  # diagonal
  if (fill_matrix) {
    for (i in 1:nInd) {
      OUT['NonMissing',i,i] <- sum(!Geno[i,] %in% Missing_codes)
      OUT['Different',i,i] <- 0
    }
  }
  
  return(OUT)
}

