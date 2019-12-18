library(tidyverse)

popProp <- read_csv("PWS472-QuaggaCol.csv")%>%
  filter(str_detect(Allele, "Dbug"))
geneCounts <- read_csv("PWS472-QuaggaCol.csv")%>%
  filter(!str_detect(Allele, "Dbug"))

# Takes a vector of frequencies and multiples each frequency
# by every other frequency except itself. Returns the sum 
# of these caluculations.
calculateExpectedHeterozygosity = function(vec){
  numAlleles = length(vec)
  h = 0
  for(i in 1:numAlleles){
    for(j in 1: numAlleles){
      if(i != j){
        h= h + (vec[i] * vec[j])
      }
    }
  }
  return(h)
}
# Create an empty matrix for the expected heterozygousity of each population for
# each of the six loci
expectedHeterozygosity = matrix(data = NA, nrow = 6, ncol = 9)
rNames = paste("Dbug", 1:6, sep = "")
rownames(expectedHeterozygosity) = rNames
cNames = c("CE-54", "EE-54", "LPI-54", "LPII-54", "WO-54", "RC-54", "EO-54", "SF-54", "Total")
colnames(expectedHeterozygosity) = cNames

#Fill in the matrix 
for(j in 1: nrow(expectedHeterozygosity)){
  locus = paste("Dbug", j, sep = "")
  totalInstance = c()
  totalCount = 0
  for(i in 1:(ncol(expectedHeterozygosity)-1)){
    pop = colnames(expectedHeterozygosity)[i]
    freq = pull(filter(popProp, str_detect(Allele, locus)), pop) #The allele frequency
    geneCount = pull(geneCounts, pop)[j] #The total number of alleles in the population
    totalInstance = c(totalInstance, freq*geneCount) 
    totalCount = totalCount + geneCount
    expectedHeterozygosity[j, i] = calculateExpectedHeterozygosity(freq) 
  }
  expectedHeterozygosity[j, "Total"] = calculateExpectedHeterozygosity(totalInstance/totalCount)
}

generateFSTmatrix = function(matrix, dbug){
  #Creates empty average heterozygousity matrix with the populations along each axis
  Hs.pairs = matrix(data = NA, nrow = 8, ncol = 8)
  rownames(Hs.pairs) = cNames[1:nrow(Hs.pairs)]
  colnames(Hs.pairs) = cNames[1:nrow(Hs.pairs)]
  
  # Fills matrix with expected subpopulation heterozygousity
  for(i in 1:nrow(Hs.pairs)){
    for(j in 1:ncol(Hs.pairs)){
      popC = colnames(Hs.pairs)[j]
      popR = rownames(Hs.pairs)[i]
      geneCountC = pull(geneCounts, popC)[dbug]
      geneCountR = pull(geneCounts, popR)[dbug]
      averageH = ((expectedHeterozygosity[dbug, popR]*geneCountR)+
                    (expectedHeterozygosity[dbug, popC]*geneCountC))/(geneCountR+geneCountC)
      Hs.pairs[popR, popC] = averageH
    }
  }
  FST.pairs = Hs.pairs
  for(i in (1:nrow(FST.pairs))){
    FST.pairs[i,] = 1 - (FST.pairs[i,]/expectedHeterozygosity[dbug, "Total"])
  }
  return(FST.pairs)
}
print("Dbug1 Pairwise Fsts")
print(generateFSTmatrix(expectedHeterozygosity, 1))
print("Dbug2 Pairwise Fsts")
print(generateFSTmatrix(expectedHeterozygosity, 2))
print("Dbug3 Pairwise Fsts")
print(generateFSTmatrix(expectedHeterozygosity, 3))
print("Dbug4 Pairwise Fsts")
print(generateFSTmatrix(expectedHeterozygosity, 4))
print("Dbug5 Pairwise Fsts")
print(generateFSTmatrix(expectedHeterozygosity, 5))
print("Dbug6 Pairwise Fsts")
print(generateFSTmatrix(expectedHeterozygosity, 6))
