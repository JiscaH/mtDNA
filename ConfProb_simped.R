#============================================================================
#============================================================================
#' @title Confidence Probabilities
#'
#' @description Estimate confidence probabilities ('backward') and assignment
#'   error rates ('forward') per category (genotyped/dummy) by repeatedly
#'   simulating genotype data from a reference pedigree using
#'   \code{\link{SimGeno}}, reconstruction a pedigree from this using
#'   \code{\link{sequoia}}, and counting the number of mismatches using
#'   \code{\link{PedCompare}}.
#'
#' @details The confidence probability is taken as the number of correct
#'   (matching) assignments, divided by all assignments made in the
#'   \emph{observed} (inferred-from-simulated) pedigree. In contrast, the false
#'   negative & false positive assignment rates are proportions of the number of
#'   parents in the \emph{true} (reference) pedigree. Each rate is calculated
#'   separatedly for dams & sires, and separately for each category
#'   (\strong{G}enotyped/\strong{D}ummy(fiable)/\strong{X} (none)) of
#'   individual, parent and co-parent.
#'
#'  This function does not know which individuals in the actual \code{Pedigree}
#'  are genotyped, so the confidence probabilities need to be added to the
#'  \code{Pedigree} as shown in the example at the bottom.
#'
#'  A confidence of \eqn{1} means all assignments on simulated data were correct for
#'  that category-combination. It should be interpreted as (and perhaps modified
#'  to) \eqn{> 1 - 1/N}, where sample size \code{N} is given in the last column
#'  of the \code{ConfProb} and \code{PedErrors} dataframes in the output. The
#'  same applies for a false negative/positive rate of \eqn{0} (i.e. to be
#'  interpreted as \eqn{< 1/N}).
#'
#' @section Assumptions:
#'   Because the actual true pedigree is (typically) unknown, the provided
#'   reference pedigree is used as a stand-in and assumed to be the true
#'   pedigree, with unrelated founders. It is also assumed that the probability
#'   to be genotyped is equal for all parents; in each iteration, a new random
#'   set of parents (proportion set by \code{ParMis}) is mimicked to be
#'   non-genotyped. In addition, SNPs are assumed to segregate independently.
#'
#'   An experimental version offering more fine-grained control is available at
#'   https://github.com/JiscaH/sequoiaExtra .
#'
#' @section Object size:
#'   The size in Kb of the returned list can become pretty big, as each of the
#'   inferred pedigrees is included. When running \code{EstConf} many times for
#'   a range of parameter values, it may be prudent to save the required summary
#'   statistics for each run rather than the full output.
#'
#' @param simped_fun  function returning a pedigree with columns id-dam-sire-birthyear-sex.
#' @param args.simped list of arguments to \code{SimPed}. \code{source()} this
#'   function first, available from https://github.com/JiscaH/sequoiaExtra .
#' @param  withLH  logical, pass in each run the BirthYear + Sex simulated by \code{SimPed}
#'   to \code{sequoia}
#' @param mt  logical, with/without \code{mtSame}. Simulated with every female
#'   founder having a unique mitochondrial haplotype. requires
#'   \code{source('infer_mt_haplotype.R')}, available from same github repo.
#' @param args.sim  list of arguments to pass to \code{\link{SimGeno}}, such as
#'   \code{nSnp} (number of SNPs), \code{SnpError} (genotyping error rate) and
#'   \code{ParMis} (proportion of non-genotyped parents). Set to \code{NULL} to
#'   use all default values.
#' @param args.seq  list of arguments to pass to \code{\link{sequoia}}, such as
#'   \code{Module} ('par' or 'ped'), \code{Err} (assumed genotyping error rate),
#'   and \code{Complex}. May include (part of) \code{SeqList}, a list of sequoia
#'   output (i.e. as a list-within-a-list). Set to \code{NULL} to use all
#'   default values.
#' @param nSim number of iterations of simulate - reconstruct - compare to
#'   perform, i.e. number of simulated datasets.
#' @param nCores number of computer cores to use. If \code{>1}, package
#'   \pkg{parallel} is used. Set to NULL to use all but one of the available
#'   cores, as detected by \code{parallel::detectCores()} (using all cores tends
#'   to freeze up your computer).
#' @param quiet suppress messages. \code{TRUE} runs \code{SimGeno} and
#'   \code{sequoia} quietly, \code{'very'} also suppresses other messages and
#'   the iteration counter when \code{nCores=1} (there is no iteration counter
#'   when \code{nCores>1}).
#'
#' @return A list, with elements:
#'   \item{ConfProb}{See below}
#'   \item{PedErrors}{See below}
#'   \item{Pedigree.reference}{the pedigree from which data was simulated}
#'   \item{LifeHistData}{}
#'   \item{Pedigree.inferred}{a list with for each iteration the inferred
#'     pedigree based on the simulated data}
#'   \item{SimSNPd}{a list with for each iteration the IDs of the individuals
#'     simulated to have been genotyped}
#'   \item{PedComp.fwd}{array with \code{Counts} from the 'forward'
#'     \code{PedCompare}, from which \code{PedErrors} is calculated}
#'   \item{RunParams}{a list with the call to \code{EstConf} as a semi-nested
#'   list (args.sim, args.seq, nSim, nCores), as well as the default parameter
#'   values for \code{SimGeno} and \code{sequoia}.}
#'   \item{RunTime}{\code{sequoia} runtime per simulation in seconds, as
#'     measured by \code{\link{system.time}()['elapsed']}.}
#'
#' Dataframe \code{ConfProb} has 7 columns:
#' \item{id.cat, dam.cat, sire.cat}{Category of the focal individual, dam, and
#'   sire, in the pedigree inferred based on the simulated data. Coded as
#'   G=genotyped, D=dummy, X=none}
#' \item{dam.conf}{Probability that the dam is correct, given the categories of
#'   the assigned dam and sire (ignoring whether or not the sire is correct)}
#' \item{sire.conf}{as \code{dam.conf}, for the sire}
#' \item{pair.conf}{Probability that both dam and sire are correct, given their
#'   categories}
#' \item{N}{Number of individuals per category-combination, across all
#'   \code{nSim} iterations}
#'
#' Array \code{PedErrors} has three dimensions:
#' \item{class}{\itemize{
#'   \item \code{FalseNeg}(atives): could have been assigned but was not
#' (individual + parent both genotyped or dummyfiable; P1only in
#' \code{PedCompare}).
#'   \item \code{FalsePos}(itives): no parent in reference pedigree, but
#' one was assigned based on the simulated data (P2only)
#'   \item \code{Mismatch}: different parents between the pedigrees
#'   }}
#' \item{cat}{Category of individual + parent, as a two-letter code where the
#'   first letter indicates the focal individual and the second the parent;
#'   G=Genotyped, D=Dummy, T=Total}
#' \item{parent}{dam or sire}
#'
#'
#' @seealso \code{\link{SimGeno}, \link{sequoia}, \link{PedCompare}}.
#'
#' @importFrom plyr adply
#'
#' @examples
#' # estimate proportion of parents that are genotyped (= 1 - ParMis)
#' sumry_grif <- SummarySeq(SeqOUT_griffin, Plot=FALSE)
#' tmp <- apply(sumry_grif$ParentCount['Genotyped',,,],
#'              MARGIN = c('parentSex', 'parentCat'), FUN = sum)
#' props <- sweep(tmp, MARGIN='parentCat', STATS = rowSums(tmp), FUN = '/')
#' 1 - props[,'Genotyped'] / (props[,'Genotyped'] + props[,'Dummy'])
#'
#' # Example for parentage assignment only
#' conf_grif <- EstConf(Pedigree = SeqOUT_griffin$Pedigree,
#'                LifeHistData = SeqOUT_griffin$LifeHist,
#'                args.sim = list(nSnp = 200,   # no. in actual data, or what-if
#'                                SnpError = 5e-3,  # best estimate, or what-if
#'                                CallRate=0.8,     # from SnpStats()
#'                                ParMis=c(0.39, 0.20)),  # calc'd above
#'                args.seq = list(Err=5e-3, Module="par"),  # as in real run
#'                nSim = 2,   # try-out, proper run >=20 (10 if huge pedigree)
#'                nCores=1)
#'
#' # parent-pair confidence, per category (Genotyped/Dummy/None)
#' conf_grif$ConfProb
#'
#' # Proportion of true parents that was correctly assigned
#' 1 - apply(conf_grif$PedErrors, MARGIN=c('cat','parent'), FUN=sum, na.rm=TRUE)
#'
#' # add columns with confidence probabilities to pedigree
#' # first add columns with category (G/D/X)
#' Ped.withConf <- getAssignCat(Pedigree = SeqOUT_griffin$Pedigree,
#'                              SNPd = SeqOUT_griffin$PedigreePar$id)
#' Ped.withConf <- merge(Ped.withConf, conf_grif$ConfProb, all.x=TRUE,
#'                       sort=FALSE)  # (note: merge() messes up column order)
#' head(Ped.withConf[Ped.withConf$dam.cat=="G", ])
#'
#' # save output summary
#' \dontrun{
#' conf_griff[['Note']] <- 'You could add a note'
#' saveRDS(conf_grif[c('ConfProb','PedComp.fwd','RunParams','RunTime','Note')],
#'    file = 'conf_200SNPs_Err005_Callrate80.RDS')
#' }
#'
#' ## P(actual FS | inferred as FS) etc.
#' PairL <- list()
#' for (i in 1:length(conf_grif$Pedigree.inferred)) {  # nSim
#'   cat(i, "\t")
#'   PairL[[i]] <- ComparePairs(conf_grif$Pedigree.reference,
#'                              conf_grif$Pedigree.inferred[[i]],
#'                              GenBack=1, patmat=TRUE, ExcludeDummies = TRUE,
#'                              Return="Counts")
#' }
#' # P(actual relationship (Ped1) | inferred relationship (Ped2))
#' PairA <- plyr::laply(PairL, function(M)
#'                      sweep(M, MARGIN='Ped2', STATS=colSums(M), FUN="/"))
#' PairRel.prop <- apply(PairA, 2:3, mean, na.rm=TRUE)  # mean across simulations
#' round(PairRel.prop, 3)
#' # or: P(inferred relationship | actual relationship)
#' PairA2 <- plyr::laply(PairL, function(M)
#'                       sweep(M, MARGIN='Ped1', STATS=rowSums(M), FUN="/"))
#'
#' @export


EstConf_simped <- function(simped_fun = SimPed,
                           args.simped = list(N_dam=100, EPP=0.1, N_offspring = c(0.2,0.3,0.3,0.2)),
                           withLH = TRUE,
                           mt = FALSE,
                           args.sim = list(nSnp = 400, SnpError = 1e-3, ParMis=c(0.4, 0.4)),
                           args.seq = list(Module="ped", Err=1e-3, Tassign=0.5, CalcLLR = FALSE),
                           nSim = 10,
                           nCores = 1,
                           quiet=TRUE)
{
  
  # check input ----
  #  if (is.null(Pedigree))  stop("Please provide Pedigree")
  #  if (is.null(LifeHistData))  stop("Please provide LifeHistData")
  if (!is.null(args.sim) & !is.list(args.sim))  stop("args.sim should be a list or NULL")
  if (!is.null(args.seq) & !is.list(args.seq))  stop("args.seq should be a list or NULL")
  if (!sequoia:::is.wholenumber(nSim) || nSim<1 || length(nSim)>1)
    stop("nSim must be a single positive number")
  if (packageVersion('sequoia') < '2.7')  mt <- FALSE
  
  if (!quiet %in% c(TRUE, FALSE, "very"))  stop("'quiet' must be TRUE, FALSE, or 'very'")
  quiet.EC <- ifelse(quiet == "very", TRUE, FALSE)
  quiet <- ifelse(quiet %in% c("very", TRUE), TRUE, FALSE)
  if (!"quiet" %in% names(args.sim))  args.sim <- c(args.sim, list(quiet = quiet))
  if (!"quiet" %in% names(args.seq))  args.seq <- c(args.seq, list(quiet = quiet))
  
  if ("Err" %in% names(args.sim)) {
    args.sim[["SnpError"]] <- args.sim[["Err"]]
    args.sim[["Err"]] <- NULL    # common confusion, otherwise fuzy matching with 'ErrorFM'.
  }
  
  Ped.ref_L <- list()
  LH_L <- vector('list', nSim)
  mtSame_L <- vector('list', nSim)
  for (i in 1:nSim) {
    Pedigree_i <- do.call(simped_fun, args.simped)
    Ped.ref_L[[i]] <- sequoia::PedPolish(Pedigree_i, KeepAllColumns=FALSE)
    if (withLH)  LH_L[[i]] <- Pedigree_i[,c(1,4,5)]   # skip dam & sire columns
    if (mt) {
      mtSame_L[[i]] <- sim_mtSame(Ped.ref_L[[i]])
    }
  }
  
  if ("Module" %in% names(args.seq)) {
    ParSib <- ifelse(args.seq$Module == "ped", "sib", "par")
  } else if ("MaxSibIter" %in% names(args.seq)) {
    ParSib <- ifelse(args.seq$MaxSibIter > 0, "sib", "par")
  } else {
    ParSib <- "sib"   # default Module = "ped"
  }
  
  if (!quiet.EC) {
    if (ParSib == "par") {
      message("Simulating parentage assignment only ...")
    } else {
      message("Simulating full pedigree reconstruction ...")
    }
  }
  
  # no. cores ----
  if (!is.null(nCores) && (!sequoia:::is.wholenumber(nCores) || nCores<1))
    stop("nCores must be a positive number, or NULL")
  
  if (is.null(nCores) || nCores>1) {
    if (!requireNamespace("parallel", quietly = TRUE)) {
      if (interactive() & !quiet.EC) {
        message("Installing pkg 'parallel' to speed things up... ")
      }
      utils::install.packages("parallel")
    }
    maxCores <- parallel::detectCores()
    if (is.null(nCores)) {
      nCores <- maxCores -1
    } else if (nCores > maxCores) {
      nCores <- maxCores
      warning("Reducing 'nCores' to ", maxCores, ", as that's all you have",
              immediate.=TRUE)
    }
    if (nCores > nSim) {
      nCores <- nSim
    }
    if (!quiet.EC)  message("Using ", nCores, " out of ", maxCores, " cores")
  }
  
  
  utils::flush.console()
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # function to simulate genotypes & infer pedigree ----
  
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  SimInfer <- function(i, RefPedigree, args.sim,
                       LHs, args.seq, mtSames,
                       quiet.EC, ParSib)
  {
    if (!quiet.EC)  cat("i=", i, "\t", format(Sys.time(), "%H:%M:%S"), "\n")
    # nothing printed by parallel::parSapply (?)
    OUT.i <- list()
    GM <- do.call(sequoia::SimGeno, c(list(Pedigree=RefPedigree[[i]]), args.sim))
    OUT.i$SimSNPd <- rownames(GM)
    args.seq.all <-  c(list(GenoM = GM,
                            LifeHistData = LHs[[i]],
                            DummyPrefix = c("sim_F", "sim_M"),
                            Plot = FALSE),
                       args.seq)
    # add mtSame separately, so that script also runs with sequoia < 2.7
    if (!is.null(mtSames[[i]]))   args.seq.all <- c(args.seq.all, mtSame = mtSames[[i]])  
    OUT.i$RunTime <- system.time(Seq.i <- do.call(sequoia::sequoia,
                                                  args.seq.all))["elapsed"]
    rm(GM)
    gc()
    if (ParSib == "par") {
      OUT.i$Pedigree.inferred <- Seq.i[["PedigreePar"]]
    } else {
      OUT.i$Pedigree.inferred <- Seq.i[["Pedigree"]]
    }
    return( OUT.i )
  }
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  # call above function, for 1 or >1 cores ----
  if (nCores>1) {
    cl <- parallel::makeCluster(nCores)
    AllOUT <- parallel::parLapply(cl, X=seq.int(nSim), fun=SimInfer,
                                  RefPedigree = Ped.ref_L,
                                  mtSames = mtSame_L,
                                  args.sim, LHs=LH_L, args.seq,
                                  quiet.EC=TRUE, ParSib)
    #                                  chunk.size=1)
    parallel::stopCluster(cl)
  } else {
    AllOUT <- plyr::llply(seq.int(nSim), .fun=SimInfer,
                          RefPedigree = Ped.ref_L,
                          mtSames = mtSame_L,
                          args.sim, LHs=LH_L, args.seq,
                          quiet.EC, ParSib)
  }
  RunTime <- sapply(AllOUT, "[[", "RunTime")
  Pedigree.inferred <- sapply(AllOUT, "[", "Pedigree.inferred")
  SimSNPd <- sapply(AllOUT, "[", "SimSNPd")
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # confidence probabilities ----
  nSimz <- ifelse(nSim>1, nSim,2)  # else problems w R auto-dropping dimension
  CatNames <- c("G", "D", "X")
  ClassNames <- c("Match", "Mismatch", "P1only", "P2only", "_")
  
  PC.rev.cd <- array(0, dim = c(nSimz, 3,3,3, 5,5),
                     dimnames = list(iter = seq_len(nSimz),
                                     id.cat = CatNames, dam.cat = CatNames, sire.cat = CatNames,
                                     dam.class = ClassNames, sire.class = ClassNames))
  for (i in 1:nSim) {
    PC.rev.cd[i,,,,,] <- sequoia::PedCompare(Ped1 = Pedigree.inferred[[i]],
                                             Ped2 = Ped.ref_L[[i]],
                                             SNPd = SimSNPd[[i]],
                                             Symmetrical=FALSE, Plot=FALSE)$Counts.detail
  }
  
  Ntot <- apply(PC.rev.cd, c('id.cat', 'dam.cat', 'sire.cat'), sum)
  OK <- list('G' = 'Match',
             'D' = 'Match',
             'X' = c('P2only', '_'))
  confA <- array(dim = c(3,3,3,3),
                 dimnames = c(list(paste0(c('dam', 'sire', 'pair'), '.conf')),
                              dimnames(PC.rev.cd)[2:4]))
  for (i in c('G','D','X')) {
    confA['dam.conf' ,,i,] <- apply(PC.rev.cd[,,i,,OK[[i]],], c('id.cat', 'sire.cat'), sum) / Ntot[,i,]
    confA['sire.conf',,,i] <- apply(PC.rev.cd[,,,i,,OK[[i]]], c('id.cat', 'dam.cat'), sum) / Ntot[,,i]
    for (j in c('G','D','X')) {
      confA['pair.conf',,i,j] <- apply(PC.rev.cd[,,i,j,OK[[i]],OK[[j]]], 'id.cat', sum) / Ntot[,i,j]
    }
  }
  
  confA[c("dam.conf" , "pair.conf"),,"X",] <- NA  # no dam
  confA[c("sire.conf", "pair.conf"),,,"X"] <- NA  # no sire
  
  Conf.df <- plyr::adply(confA, .margins=2:4)
  Conf.df <- merge(Conf.df,
                   plyr::adply(Ntot, .margins=3:1, function(x) data.frame(N=x)))
  Conf.df <- Conf.df[Conf.df$id.cat != 'X',]
  if (ParSib == "par") {
    Conf.df <- Conf.df[Conf.df$id.cat == 'G' & Conf.df$dam.cat %in% c('G','X') &
                         Conf.df$sire.cat %in% c('G','X'), ]
  }
  Conf.df <- Conf.df[order(Conf.df$id.cat, Conf.df$dam.cat, Conf.df$sire.cat), ]
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # assignment errors (ignores co-parent) ----
  PedComp.fwd <-  array(0, dim=c(nSimz, 7,5,2),
                        dimnames = list(iter = seq_len(nSimz),
                                        cat = c("GG", "GD", "GT", "DG", "DD", "DT", "TT"),
                                        class = c("Total", "Match", "Mismatch", "P1only", "P2only"),
                                        parent = c("dam", "sire")))
  for (i in 1:nSim) {
    PedComp.fwd[i,,,] <- sequoia::PedCompare(Ped1 = Ped.ref_L[[i]],
                                             Ped2 = Pedigree.inferred[[i]],
                                             SNPd = SimSNPd[[i]],
                                             Symmetrical=FALSE, Plot=FALSE)$Counts
  }
  
  PedComp.tmp <- apply(PedComp.fwd, 2:4, sum)
  PedErrors <- sweep(PedComp.tmp[,c("P1only", "P2only","Mismatch"),], c(1,3),
                     PedComp.tmp[,"Total",], "/")
  PedErrors[c("GG", "GD", "DG", "DD"), "P2only", ] <- NA  # if parent in ref. pedigree, by def not P2only
  dimnames(PedErrors)[['class']] <- c("FalseNeg", "FalsePos", "Mismatch")
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # out ----
  RunParams <- list(EstConf = sequoia:::namedlist(args.sim, args.seq, nSim, nCores),
                    SimGeno_default = formals(SimGeno),
                    sequoia_default = formals(sequoia),
                    sequoia_version = as.character(utils::packageVersion("sequoia")))
  
  return( list(ConfProb = Conf.df,
               PedErrors = PedErrors,
               Pedigree.reference = Ped.ref_L,
               LifeHistData = LH_L,
               Pedigree.inferred = Pedigree.inferred,
               SimSNPd = SimSNPd,
               PedComp.fwd = PedComp.fwd,
               RunParams = RunParams,
               RunTime = RunTime) )
}
