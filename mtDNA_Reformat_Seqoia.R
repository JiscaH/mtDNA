#load libraries needed
library(vcfR)
library(adegenet)
library(sequoia)
library(graph4lg)

############################################################################
######## Anya's code to convert vcf file to dataframe of haplotypes ########
############################################################################

#Load vcf file using vcfR
Vcfbbear_mtDNA <- read.vcfR("D:/Hervey PhD Materials/Black Bear GTseq/bbear_mtDNA_calls.vcf")

#Convert vcf to genind object using adegent
##NOTE** added ploidy = 1 to read data as haploid not diploid
Gen_bbear_mtDNA <- vcfR2genind(Vcfbbear_mtDNA, ploidy = 1) # took some time but worked!
Gen_bbear_mtDNA <- vcfR2genind(Vcfbbear_mtDNA) # took some time but worked!

#Convert genind to dataframe in terms of genotype info
Gen_bbear_mtDNA_df <- genind2df(Gen_bbear_mtDNA)

#write dataset to csv for reference
write.csv(Gen_bbear_mtDNA_df,"D:/Hervey PhD Materials/Black Bear GTseq/BBear_Haplotypes.csv")

###################################################################################################
######## Code to convert dataframe of haplotypes into pairwise matrix of matches/mismatches #######
###################################################################################################

#we want to keep the IDs stored as a vector to add to a final dataframe later
IDs<-rownames(Gen_bbear_mtDNA_df)

#To create a final pairwise matrix between individuals summing the number of mismatches 
#we will first create an empty list to store the pairwise matrices for each locus
pairwise_matrices <- list()

#Next we will use a for loop to genearate a pairwise matrix
#for each locus where the value is 0 for matching, 1 for mismatching
#or NA is at least one of the two individuals compared are not genotyped
#at the locus of interest
for (col in colnames(Gen_bbear_mtDNA_df)) {
  #Get the locus column values
  values <- Gen_bbear_mtDNA_df[[col]]
  
  #Initialize an empty matrix
  matrix_size <- length(values)
  pairwise_matrix <- matrix(NA, nrow = matrix_size, ncol = matrix_size)
  
  #Generate pairwise matrix where 0 is a match, 1 is a mismatch, and NA if 
  #at least one individual is not genotyped at the locus of interest
  for (i in 1:matrix_size) {
    for (j in 1:matrix_size) {
      if (is.na(values[i]) || is.na(values[j])) {
        pairwise_matrix[i, j] <- NA
      } else if (values[i] == values[j]) {
        pairwise_matrix[i, j] <- 0
      } else {
        pairwise_matrix[i, j] <- 1
      }
    }
  }
  #Change diagnol value between the same sample as NA
  diag(pairwise_matrix) <- NA
  
  #Add pairwise matrix to the list
  pairwise_matrices[[col]] <- pairwise_matrix
}

#Now that we have a list of matrices where each matrix contains the pairwise comparisons between all samples
#and 1 represents mismatches genotype, 0 represents matching genotype and NA represents at least one sample
#for the comparison was not gentoyped at the locus
#now we want to sum the number of mismatches between each sample but to do this we need to convert NAs to 0s
#using a custom function and then include it in the reduce function which will be used to add across the list 
#of matrices
modifiedSum <- function(x, y) {
  replace(x, is.na(x), 0) + replace(y, is.na(y), 0)
}
pairwise_matrix_mismatch<-Reduce(modifiedSum, pairwise_matrices)

#convert back to dataframe and add IDs back
pairwise_matrix_mismatch_df<-as.data.frame(pairwise_matrix_mismatch)
colnames(pairwise_matrix_mismatch_df)<-IDs
rownames(pairwise_matrix_mismatch_df)<-IDs

#print pairwise matrix of mismatches to csv file
write.csv(pairwise_matrix_mismatch_df,"D:/Hervey PhD Materials/Black Bear GTseq/BBear_pairwise_mismatches.csv")

#now create pairwise matrix for Sequoia format where mismatching haplotypes are assigned as 0 and 
#matching are assigned as 1
pairwise_matrix_mismatch_df[pairwise_matrix_mismatch_df==0]<- -9
pairwise_matrix_mismatch_df[pairwise_matrix_mismatch_df>0]<- 0
pairwise_matrix_mismatch_df[pairwise_matrix_mismatch_df==-9]<- 1

#Output as pairwise matrix for sequoia
write.csv(pairwise_matrix_mismatch_df,"D:/Hervey PhD Materials/Black Bear GTseq/BBear_pairwise_mismatches_sequoia.csv")

###################################################################
######## Anya's code for reformatting to be read by Sequoia########
###################################################################

##NOTE** had to slightly adjust code after changing how the vcf file was read into adegenet as a haploid dataset

#first change NAs to -9 for Sequoia format
Gen_bbear_mtDNA_df[is.na(Gen_bbear_mtDNA_df)] <- -9

#convert dataframe to numeric using a function to apply as.numeric across all haplotypes
i<-c(1:ncol(Gen_bbear_mtDNA_df))
Gen_bbear_mtDNA_df[,i]<-apply(Gen_bbear_mtDNA_df[,i],2,
                              function(x)as.numeric(as.character(x)))
#convert to matrix
Gen_bbear_mtDNA_df<-as.matrix(Gen_bbear_mtDNA_df)

# Check format for Sequoia
GenoM_bbear <- Gen_bbear_mtDNA_df
CheckGeno(GenoM_bbear)
