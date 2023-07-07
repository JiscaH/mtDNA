# Infer mitochondrial haplotype for as many individuals as possible based on
# a pedigree + known mt haplotype for some individuals. 
# Infers for both ancestors, descendants, and other matrilineal relatives. 

# IN:
# Pedigree: data.frame with columns id - dam - sire
# mtHaps: named vector

# OUT:
# mtHaps extended with haplotypes for additional individuals, if any were found.


infer_mt_haplotype <- function(Pedigree, mtHaps) {
  Ped <- sequoia::PedPolish(Pedigree)
  Ped$G <- sequoia::getGenerations(Ped)
  nG <- max(Ped$G)
  
  # go bottom -> top to derive haplotypes of ancestors, 
  # then top -> bottom to get haplotypes of all their descendants
  gg <- c(nG:1, 1:nG)  
  for (x in seq_along(gg)) {
    these <- with(Ped, G==gg[x] & !is.na(dam))
    if (x <= nG) {  # way up
      mtHaps_new <- with(Ped, setNames( mtHaps[id[these]], dam[these] ) )
      mtHaps_new <- mtHaps_new[!duplicated(names(mtHaps_new))]
    } else {  # way down
      mtHaps_new <- with(Ped, setNames( mtHaps[dam[these]], id[these] ) )
    }
    if (all(is.na(mtHaps_new)))  next   # especially possible on way up
  
    # check for conflicts
    z <- intersect(names(mtHaps), names(mtHaps_new))
    if (any(mtHaps[z] != mtHaps_new[z])) {
      stop('Conflicting mt haplotype for individuals: ', names(which(mtHaps[z] != mtHaps_new[z])),
           ' Please check maternal pedigree links and/or mt haplotypes!')  
    }
    
    mtHaps <- c(mtHaps, mtHaps_new)
    mtHaps <- mtHaps[!duplicated(names(mtHaps))]
  }
  
  return( mtHaps )
}
