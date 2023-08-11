#' @title Simulate a pedigree
#'
#' @description Create a pedigree for a single cohort of a socially monogamous
#'    population, with a user-setable rate of extra-pair paternities. 
#'
#' @param N_dam  number of dams (= number of nests)
#' @param EPP  extra-pair paternity (EPP) rate, between 0 and 0.5.
#' @param N_offspring vector with proportions giving the distribution of number 
#'   of offspring per female (=per nest) 
#' @param withLH logical, also create columns with BirthYear and Sex?
#'
#' @return dataframe with pedigree with columns id-dam-sire. If withLH=TRUE, 
#'  also columns BirthYear + Sex.
#'
#' @details  Can only simulate a single cohort with unrelated parents yet.
#'
#' @examples
#'  Ped <- SimPed(N_dam=200, EPP=0.25, N_offspring=c('1'=0.3, '2'=0.69, '3'=0.01))
#'

# litter size as input vector with proportions
# single extra-pair male per female; no restriction on # females per male
# withLH: also simulate birthyear & sex column


SimPed <- function(N_dam = 100,
                   EPP = 0.1,
                   N_offspring = c('1'=0.2, '2'=0.3, '3'=0.3, '4'=0.1, '5'=0.1),
                   withLH = TRUE)
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
  dam_id <- paste0('i_f', formatC(1:N_dam, width=nchar(N_dam), flag='0'))
  sire_id <- paste0('i_m', formatC(1:N_dam, width=nchar(N_dam), flag='0'))
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
  
  if (withLH) {
    Pedigree$BirthYear <- with(Pedigree, ifelse(id %in% c(dam, sire), 1, 2))
    Pedigree$Sex <- with(Pedigree, ifelse(id %in% dam, 1, 
                                          ifelse(id %in% sire, 2, 3)))
    Pedigree$Sex[Pedigree$Sex==3] <- rbinom(n=sum(Pedigree$Sex==3), size=1, prob=0.5) +1
  }

  return( Pedigree )
}
