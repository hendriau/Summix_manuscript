## Simulation code

library(gtools)
library(dplyr)
library(tidyr)

# initialize data frame with probability values
probvalues = data.frame(
  startp = c(0,0.01,0.05,0.1,0.25),
  endp = c(0.015,0.055,0.105,0.255,0.505)
)

# Ancestries used

Ancestries = c('AFR', 'EAS', 'EUR', 'IAM', 'SAS')

# Number of Ancestries simulated
# IMPORTANT: Generates values by doing 1 - (# of ancestries)
# IMPORTANT: Sum of first n-1 values is less then 1, then uses values and last value is 1-(first values)

AncestryNum = 5

# generate data frame of all possible probability values, 2 to 5 ancestries
if (AncestryNum == 2){
  tf = data.frame('X1' = probvalues$endp)
}
if (AncestryNum > 2){
  tf = data.frame(combinations(5, (AncestryNum - 1), probvalues$endp, repeats.allowed = TRUE))
}

# clean data frame with probability values to only sum < 1, ending prob values
tfc = tf %>% 
  mutate(sumc = rowSums(.)) %>% 
  filter(sumc < 1) %>% 
  select(starts_with('X'))

# create starting prob values data frame
tfs = data.frame(matrix(0, nrow = dim(tfc)[1], ncol = dim(tfc)[2]))

# match starting values to their ending values
for (n in (1:dim(tfc)[2])){
  for (i in (1:dim(tfc)[1])){
    tfs[i,n] = probvalues$startp[which(probvalues$endp == tfc[i,n])]
  }
}

# create ending probability frame for final ancestry
tfsf = (1 - rowSums(tfc)); tfcf = (1 - rowSums(tfs))
tfs$f = tfsf; tfc$f = tfcf

# name columns, e.g. S1 = starting ancestry 1, E1 = ending ancestry 1
names(tfs) = paste('S', (1:dim(tfs)[2]), sep = '')
names(tfc) = paste('E', (1:dim(tfs)[2]), sep = '')

# merge starting/ending probabilies data frame in correct order, S1 E1 S2 E2...
neworder <- order(c(2*(seq_along(tfs) - 1) + 1,
                    2*seq_along(tfc)))
tfm = cbind(tfs, tfc)[,neworder]

# create data frame of all ancestry permutations
ancpf = data.frame(permutations(5, AncestryNum, Ancestries))
names(ancpf) = paste('Anc', (1:dim(ancpf)[2]), sep = '')

# create data frame for ancestry and probability combinations merge
tframe = data.frame()

#merge all possible ancestry and probability combinations
for (j in 1:dim(ancpf)[1]){
  
  # reusable frame
  sframe = data.frame(matrix(nrow = dim(tfm)[1]))
  
  # merge 3 data frames together in correct order
  for (i in 1:dim(ancpf)[2]){
    ancvec = rep(ancpf[j,i], dim(tfm)[1])
    cframe = data.frame(ancvec, get('tfm')[[paste('S', i, sep = '')]], get('tfm')[[paste('E', i, sep = '')]])
    sframe = cbind(sframe, cframe)
  }
  
  # remove artifact vector
  sframe[,1] = NULL
  
  # bind to final frame
  tframe = rbind(tframe, sframe)
}

# name final frame correctly
# not generalized, comment out for different ancestries
names(tframe) = c('A1', 'S1', 'E1'
                  ,'A2', 'S2', 'E2'
                  ,'A3', 'S3', 'E3'
                  ,'A4', 'S4', 'E4'
                  ,'A5', 'S5', 'E5'
)

# bind together combinations for duplicate testing
# not generalized, comment out for different ancestries
r1 = tframe %>% 
  unite(t1, 1:3) %>% 
  unite(t2, 2:4) %>% 
  unite(t3, 3:5) %>% 
  unite(t4, 4:6) %>% 
  unite(t5, 5:7)

# final frame of all combinations
# removed all duplicate combinations
result = tframe[!duplicated(t(apply(r1, 1, sort))),]

# write result
write.table(result,
            file = '/parameters/2_anc_parameters.txt',
            row.names = FALSE)

# parameter file for 1 ancestry simulaitons

ancone = data.frame(
  A1 = c('AFR', 'EAS', 'EUR', 'IAM', 'SAS'),
  S1 = rep(1, 5),
  E1 = rep(1, 5)
)

write.table(ancone,
            file = '/parameters/1_anc_parameters.txt',
            row.names = FALSE)
