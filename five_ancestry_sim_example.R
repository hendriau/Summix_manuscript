## Simulation Script

# Load in Libraries
library(dplyr)
library(nloptr)

# Number which calls specific paramaters for simulation, increasing
ivalnum = 1

# Number of simulations to run
testnum = 1000
# Number of people to simulate
ntot = 10000

# Load in reference data and paramaters
load("/referencedata/reference.Rdata")
parameters = read.csv("/parameters/5_anc_parameters.txt", sep="")

# Set limits for each cluster core
nodecontrol = ceiling(dim(parameters)[1]/20)
tivec = c(1, round(seq(1, dim(parameters)[1], by = dim(parameters)[1]/nodecontrol))[-1], dim(parameters)[1])

# Ancestry estimation function
# Loaded in manually so each core on cluster would have access to it
ancestr = function(refmatrix, obsvector){
  
  testmatrix = cbind(refmatrix, obsvector)
  
  starting = numeric(ncol(refmatrix))
  for (i in 1:(ncol(refmatrix))){
    starting[i] = 1/ncol(refmatrix)
  }
  
  fn.ancmix = function(x){
    minfunc = 0
    for (i in 1:ncol(refmatrix)){
      minfunc = minfunc + x[i]*testmatrix[i]
    }
    minfunc = minfunc - testmatrix[ncol(refmatrix) + 1]
    minfunc = sum((minfunc)**2)
    return(minfunc)
  }
  
  gr.ancmix <- function(x){
    gradvec = matrix(0,ncol(refmatrix),1)
    gradfunc = 0
    for (i in 1:ncol(refmatrix)){
      gradfunc = gradfunc + x[i]*testmatrix[i]
    }
    gradfunc = gradfunc - testmatrix[ncol(refmatrix) + 1]
    for (i in 1:ncol(refmatrix)){
      gradvec[i] = sum(2 * testmatrix[i] * gradfunc)
    }
    return(gradvec)
  }
  
  heq.ancmix = function(x){
    equality = 0
    for (i in 1:ncol(refmatrix)){
      equality = equality + x[i]
    }
    return(equality - 1)
  }
  
  hin.ancmix <- function(x){
    h = numeric(ncol(refmatrix))
    for (i in 1:ncol(refmatrix)){
      h[i] = x[i]
    }
    return(h)
  }
  
  start_time = Sys.time()
  suppressMessages({S = slsqp(starting,
            fn = fn.ancmix,
            gr = gr.ancmix,
            hin = hin.ancmix,
            heq = heq.ancmix)})
  end_time = Sys.time()
  ttime = end_time - start_time
  
  val = c( S$par,
           S$value,
           S$iter,
           ttime
  )
  
  return(val)
}

# Determines proportions of population
popvecprop = function(propvec){
  prop = numeric(length = (dim(propvec)[2]/2))
  
  starts = propvec %>% 
    select(starts_with('S'))
  ends = propvec %>% 
    select(starts_with('E'))
  
  for (i in (1:length(prop))){
    prop[i] = runif(1, min = starts[1,i], max = ends[1,i])
  }
  
  prop[length(prop)] = 1 - sum(prop[-length(prop)])
  
  return(prop)
}

# Pulls reference data from parameters
AncFrame = data.frame(parameters[,grep('A', names(parameters))])
# Pulls Observed data from parameters
SEframe = data.frame(parameters[,grep('S|E', names(parameters))])
# Calculates observed intervals
testint = c(ifelse(ivalnum == 1, tivec[ivalnum], tivec[ivalnum] + 1), tivec[ivalnum+1])
# Initialize output data frame
finalframe = data.frame(matrix(vector(), 0, (21 + dim(parameters)[2]),
                               dimnames=list(c(), c('P_Num', 'T_Num', 'Seed', names(parameters),
                                                    'AFR_prop', 'EAS_prop', 'EUR_prop', 'IAM_prop', 'SAS_prop',
                                                    'R_AFR', 'R_EAS', 'R_EUR', 'R_IAM', 'R_SAS', 'R_val', 'R_iter', 'R_time',
                                                    'R_AFR_acc', 'R_EAS_acc', 'R_EUR_acc', 'R_IAM_acc', 'R_SAS_acc'))),
                        stringsAsFactors=F)


# Simulations
for (m in testint[1]:testint[2]){
  
  # Pull data from reference
  propvec = SEframe[m,]
  ancvec = AncFrame[m,]
  ancnum = dim(AncFrame)[2]
  outframe = finalframe[FALSE,]
  
  for (j in 1:testnum){
    
    # Set Seed
    seed = as.integer(Sys.time())
    set.seed(seed)
    
    # Select 100K SNPS from reference data
    refdat = reference1000GNAM %>% 
      sample_n(100000) %>% 
      select(SNP, CHR, A1, A2, AF_afr, AF_eas, AF_eur, namAFflip, AF_sas) %>% 
      rename(rsID = SNP, AFR = AF_afr, EAS = AF_eas, EUR = AF_eur, IAM = namAFflip, SAS = AF_sas)
    
    # Generate vector of proportions of ancestry
    prop = popvecprop(propvec)
    
    # Number of people in simulation population of each ancestry
    popnum = numeric(ancnum)
    for (i in 1:ancnum){
      if (i == ancnum){
        popnum[i] = (ntot - sum(popnum))
      }
      else {
        popnum[i] = floor((ntot * prop[i]))
      }
    }
    
    # Initialize matrix for simulated population
    pop_matrix = matrix(0, nrow = 100000, ncol = 3)
    
    # Generate allele counts for population using multinom 
    for (i in 1:ancnum){
      popmatrixadd = t(sapply(refdat[[as.character(ancvec[1,i])]], function(x){x2<-as.numeric(x); rmultinom(1, popnum[i], prob=(c(x2**2, 2*x2*(1-x2), (1-x2)**2)))}))
      pop_matrix = pop_matrix + popmatrixadd
    }
    
    # MAF threshold to filter by
    MAF_thresh = 0.01
    
    # Calcluate allele frequencies
    master_frame_gen1 <- data.frame(refdat[,c(1:4)], pop_matrix)
    master_frame_gen1$AF <- (2 * master_frame_gen1[,5] + master_frame_gen1[,6]) / (2 * ntot)
    # Filter by MAF
    master_frame_gen2 <- master_frame_gen1[master_frame_gen1$AF > MAF_thresh & master_frame_gen1$AF < (1-MAF_thresh),]
    
    # Pull reference data
    refdatm = refdat %>% 
      select(rsID, AFR, EAS, EUR, IAM, SAS)
    # Set Observed data
    obsvecm = master_frame_gen2 %>% 
      select(rsID, AF)
    
    # Merge reference and observed
    mergeframe = merge(refdatm, obsvecm, by = 'rsID') %>% 
      select(AFR, EAS, EUR, IAM, SAS, AF)
    
    # Run SQP algorithm
    Rslsqp = ancestr(mergeframe[,1:5], mergeframe[,6])
    
    # Save test info for each simulation
    testinfo = data.frame(P_Num = m, T_Num = j, Seed = seed)
    # Save parameters
    parinfo = parameters[m,]
    # Convert to character
    ancvec[] <- lapply(ancvec, as.character)
    
    # Save data into final data frame and output
    propinfo = data.frame(AFR_prop = ifelse('AFR' %in% ancvec, prop[which(ancvec == 'AFR')], 0),
                          EAS_prop = ifelse('EAS' %in% ancvec, prop[which(ancvec == 'EAS')], 0),
                          EUR_prop = ifelse('EUR' %in% ancvec, prop[which(ancvec == 'EUR')], 0),
                          IAM_prop = ifelse('IAM' %in% ancvec, prop[which(ancvec == 'IAM')], 0),
                          SAS_prop = ifelse('SAS' %in% ancvec, prop[which(ancvec == 'SAS')], 0))
    Rinfo = data.frame(R_AFR = Rslsqp[1],
                       R_EAS = Rslsqp[2],
                       R_EUR = Rslsqp[3],
                       R_IAM = Rslsqp[4],
                       R_SAS = Rslsqp[5],
                       R_val = Rslsqp[6],
                       R_iter = Rslsqp[7],
                       R_time = Rslsqp[8])
    # Accuracy, "True proportions" - Estimated
    Rest = data.frame(R_AFR_acc = (Rinfo[,'R_AFR'] - propinfo[,'AFR_prop']),
                      R_EAS_acc = (Rinfo[,'R_EAS'] - propinfo[,'EAS_prop']),
                      R_EUR_acc = (Rinfo[,'R_EUR'] - propinfo[,'EUR_prop']),
                      R_IAM_acc = (Rinfo[,'R_IAM'] - propinfo[,'IAM_prop']),
                      R_SAS_acc = (Rinfo[,'R_SAS'] - propinfo[,'SAS_prop']))
    
    # Bind Data
    outline = cbind(testinfo, parinfo, propinfo, Rinfo, Rest)
    outframe = rbind(outframe, outline)
  }
  
  # Bind Data
  finalframe = rbind(finalframe, outframe)
  
}

# Write output to text file
write.table(finalframe,
            file = paste("/5ancestry/anc",dim(AncFrame)[2], "_testnum_", ivalnum, ".txt", sep = ''),
            row.names = FALSE)
