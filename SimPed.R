# simulate pedigree
# single cohort, socially monogamous, EPP between 0 and 0.5
# litter size as input vector with proportions
# single extra-pair male per female; no restriction on # females per male


SimPed <- function(N_dam = 100,
                   EPP = 0.1,
                   N_offspring = c('1'=0.2, '2'=0.3, '3'=0.3, '4'=0.1, '5'=0.1))
{
  # check input ====
  if (EPP < 0 | EPP > 0.5)  stop('EPP must be between 0 and 0.5')
  if (!sequoia:::is.wholenumber(N_dam))  stop('N_dam must be a whole number')
  if (sum(N_offspring) != 1) {
    N_offspring <- N_offspring / sum(N_offpsring)
    warning("scaled N_offpsring to sum to 1, new values: ", round(N_offspring))
  }
  if (!is.null(names(N_offspring))) {
    nums <- suppressWarnings(as.numeric(names(N_offspring)))
    if (any(is.na(nums)) || any(!sequoia:::is.wholenumber(nums))) {
      stop('names of N_offspring must be (quoted) whole numbers')
    }
  } else {
    names(N_offspring) <- seq_len(length(N_offspring))
  }

  # sample number of offspring by each female
  litter_size <- sample(x = as.numeric(names(N_offspring)),
                      size = N_dam,
                      replace = TRUE,
                      prob = N_offspring)

  # sample which offspring are sired by extra-pair male
  dad_EP <- rbinom(sum(litter_size), size=1, prob=EPP)

  # sample extra-pair male for each female
  dam_id <- paste0('f', formatC(1:N_dam, width=nchar(N_dam), flag='0'))
  sire_id <- paste0('m', formatC(1:N_dam, width=nchar(N_dam), flag='0'))
  EP_mate <- setNames(rep(NA, N_dam), dam_id)
  for (i in 1:N_dam) {
    EP_mate[i] <- sample(sire_id[-i], size=1)   # exclude social male from options
  }

  # create pedigree
  Pedigree <- data.frame(id = paste0("i", formatC(1:sum(litter_size),
                                                  width=nchar(sum(litter_size)),
                                                  flag="0")),
                         dam = rep(dam_id, times = litter_size),
                         sire = NA)

  for (i in 1:nrow(Pedigree)) {
    if (dad_EP[i]==0) {
      Pedigree$sire[i] <- gsub('f', 'm', Pedigree$dam[i])
    } else {
      Pedigree$sire[i] <- EP_mate[Pedigree$dam[i]]
    }
  }

  return( Pedigree )
}
