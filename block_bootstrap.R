################
##  Function that performs block bootstrap of genome/exome data.
##  An example call is located at the bottom of this file.  
##  
##  Inputs:
##    - data::DataFrame
##        Dataframe containing named columns from standardized dataset of 
##        reference, gnomAD, and the corresponding genetic block number. 
##        Must contain columns CHR, RSID, POS, cM_Block, ALT, REF, and any numeric
##        columns referenced in the `ancestries` and `target_anc` arguments.
##
##    - ancestries::Vector::String
##        List of column names of the reference ancestries in the `data` argument
##
##    - target_anc::String
##        The target ancestry to generate confidence intervals for
##  
##    - lower_upper_quantiles::Vector::Numeric
##        The probabilities for the lower and upper quantiles of the confidence
##        interval, in order of lower, then upper. 
##    
##    - NRESAMPLES::Integer
##        The number of resamples to perform for the bootstrap process
##
##    - NWORKERS::Integer
##        The number of processes/cores/sessions to spawn for calculating resamples.
##        To fully leverage the `future` package for multicore processing, run 
##        a script which "sources" this file and call from the command line
##        with Rscript.  Generally, each worker will consume a fixed amount of 
##        memory that will be between 600 MiB and 4 GiB.  Adjust the number of 
##        workers accordingly for your workstation.    
##
##  RetVals:
##    - Dataframe containing the upper and lower quantiles from the bootstrap process. 
##      For example: 
##  +----------------+----------+--------------------+--------------------+
##  | ref_ancestries | estimate | lower_0025         | upper_0975         |
##  +----------------+----------+--------------------+--------------------+
##  | ref_AF_eur     | 0.1234   | 0.78999650329E097  | 0.795453823816841  |
##  +----------------+----------+--------------------+--------------------+
##  | ref_AF_afr     | 0.1234   | 0.0456665228368799 | 0.0477531318024074 |
##  +----------------+----------+--------------------+--------------------+
##  | ref_AF_sas     | 0.1234   | 0.0799985730403849 | 0.0873194853156217 |
##  +----------------+----------+--------------------+--------------------+
##  | ref_AF_eas     | 0.1234   | 0.0315489944964883 | 0.0356779163897794 |
##  +----------------+----------+--------------------+--------------------+
##  | ref_AF_iam     | 0.1234   | 0.0419186747641208 | 0.0445834200423259 |
##  +----------------+----------+--------------------+--------------------+
################            # Stacked Allele Frequency Table, {RSID, CHR, POS, REF, ALT, AF_1, AF_2, ...}
block_bootstrap <- function(genome_data = NULL,      
                            # Proportion/Reference Ancestries - Allele Freq. Columns
                            ancestries = c("ref_AF_afr_1000G",
                                           "ref_AF_iam_1000G",
                                           "ref_AF_eur_1000G", 
                                           "ref_AF_eas_1000G",
                                           "ref_AF_sas_1000G"),
                            # (to-)Estimate Ancestry - Allele Freq. Column
                            target_anc = "gnomad_AF_oth",
                            # Probabilities, boundaries of Confidence Interval - resampling-replicate distribution of est.
                            lower_upper_quantiles = c(0.025, 0.975),
                            # Quantity of replicates with resampling to compute
                            NRESAMPLES = 1000,
                            # Qty. processing units allocated - advanced, refer to documentation
                            NWORKERS = 1,
                            # Suppress output, warnings, messages if FALSE
                            debug = FALSE)
{
  ################################
  ## Library Init
  ################################
  library(summix, quietly = !debug)
  library(progress, quietly = !debug)
  library(tidyverse, quietly = !debug)
  library(future, quietly = !debug)
  ## These should be loaded by tidyverse
  #library(tidyr) 
  #library(tibble)
  
  start_time = Sys.time()

  ################################
  ## Split by cM_Block
  ################################
  print("Generating cM Blocks")
  
  # Only hold onto select columns
  genome_data = genome_data[,c("CHR", "RSID", "POS", "cM_Block", "ALT", "REF", ancestries, target_anc)]
  
  ## Save the row-index of each SNP for bootstrap, later
  genome_data$row_idx = seq.int(nrow(genome_data)) 
  
  ## create a unique block identifier for each chromosome and block pair
  genome_data$cM_Block2<-paste(genome_data$CHR, genome_data$cM_Block, sep=".")
  ## Generate each 'block'
  blocks <- split(genome_data, genome_data$cM_Block2)
  
  ################################
  ## Do the Bootstrap
  ################################
  print("Beginning Resampling dispatch")
  
  plan(multiprocess, workers = NWORKERS)
  print(paste0("Available cores: ", availableCores()))
  print(paste0("Using ", NWORKERS, " workers."))
  options(future.globals.maxSize= 0.3*1024^3) ## Set max copyable size for each future to ~300 MiB
  
  allfutures = NULL;
  allfutures = list(1:NRESAMPLES)
  
  # Bootstrap Loop
  pb <- progress_bar$new(total = NRESAMPLES)
  for (i in 1:NRESAMPLES) {
    pb$tick()
    
    allfutures[[i]] <- future({
      ## Sample block numbers, with replacement
      cm = c( sample.int(length(table(genome_data$cM_Block2)),
                         size = length(table(genome_data$cM_Block2)),
                         replace = TRUE) )
      
      ## From resampled w/ replacement blocks, get
      ## indicies of SNPs to fill dataset
      indices = unlist(sapply(blocks[cm], function(x){ x$row_idx }))
      
      sample <- genome_data[indices,]  ## dataset for this sampling
      sample = sample[,c("CHR", "RSID", "POS", "ALT", "REF", ancestries, target_anc)]
      # sample = tidyr::drop_na(sample, any_of(target_anc))
      sample = sample[!is.na(sample[,target_anc]),]
      
      ## Initial Values
      guess <- c(.2, .2, .2, .2, .2)
      
      out = summix(data = sample,
                    reference = ancestries,
                    observed = target_anc, 
                    pi.start = guess)
      
      cm = NULL;
      indices = NULL;
      sample = NULL;
      gc()
      
      out
    })
  }
  print("Done!")
  
  # Do Garbage Collection on data after use
  rm(blocks)
  gc()
  
  print("Beginning DMS Evaluation of Resample")
  
  ## Now pull values from the evaluated futures
  res = NULL;
  pb <- progress_bar$new(total = NRESAMPLES)
  for (i in 1:NRESAMPLES) {
    pb$tick()
    
    while(!resolved(allfutures[[i]])) {
      # Dead Man's Switch
    }
    if (is.null(res)) {
      res = value(allfutures[[i]])
    }
    else {
      res = rbind(res, value(allfutures[[i]]))
    }
  }
  print("Done!")
  
  ################################
  ## Get the Quantiles 
  ################################
  print("Beginning quantiles")
  
  ## Identify and create `str(...)` column names for user-supplied 
  ## quantile/Conf.Int endpoints.  
  lower_name = paste0("lower_", gsub('\\.', '', toString(lower_upper_quantiles[[1]])))
  upper_name = paste0("upper_", gsub('\\.', '', toString(lower_upper_quantiles[[2]])))
  
  ## Setup output dataframe
  .nrow = length(ancestries)
  estimate_info = tibble(ref_ancestries = ancestries,
                         estimate = numeric(.nrow),
                         mean = numeric(.nrow),
                         median = numeric(.nrow),
                         "{lower_name}" := numeric(.nrow),
                         "{upper_name}" := numeric(.nrow))
  
  ## Perform proportion estimation across entirety of input RSIDs and CHRs
  guess <- c(.2, .2, .2, .2, .2)
  out = summix(data = genome_data,
                reference = ancestries,
                observed = target_anc, 
                pi.start = guess)
  
  ## Construct values for the estimate info outputs.
  ## Performs maps across the values of "ref_ancestries" to form quantities.  
  anc = "ref_ancestries"  # need to use indirection within `map` functions
  estimate_info <- estimate_info %>%
    mutate(estimate = map_dbl(.data[[anc]], ~{ out[[.]] })) %>%
    mutate(mean = map_dbl(.data[[anc]], ~{ mean(res[[.]]) })) %>%
    mutate(median = map_dbl(.data[[anc]], ~{ median(res[[.]]) })) %>%
    mutate("{lower_name}" := map_dbl(.data[[anc]], ~{ 
      quantile(res[[.]], probs=lower_upper_quantiles[[1]]) 
    })) %>%
    mutate("{upper_name}" := map_dbl(.data[[anc]], ~{ 
      quantile(res[[.]], probs=lower_upper_quantiles[[2]]) 
    }))
  
  print(paste0("Total time for ", target_anc))
  print(Sys.time() - start_time)
  
  return(list(estimate_info=estimate_info, replicates=res))
}  


################
##  An example function which calls the bootstrap function above
##
##  Inputs: 
##    - est_ancestries::Vector::String
##        List of all ancestries to generate confidence intervals for
##    - datadir::String
##        Filepath of where origin dataset exists
##    - name_template::String
##        The naming convention of all the origin dataset files.  Place a % character
##        where the chromosome number is. 
##    - outdir::String
##        Filepath of where to write the output dataframes
##    - nworkers::Integer
##        Number of workers to use, passed to the bootstrap function argument.  
##        Refer to that function documentation above for an explanation. 
##
##  Outputs:
##    - Writes out dataframe of the bootstrap confidence intervals at locations
##      specified in arguments.  
##
##  Notes:
##    - To run on all ancestries, try loading Chrom 21, and get all the names for ``est_ancestries``
##      argument, as in the following example:
##      ```r
##          library(R.utils)
##          library(data.table)
##          
##          datadir = "/ref_gnomad_mergedata/genome_gnomad-ref-bherer_all_v2/"
##          name_template = "genome_cmInfo_chrom%_output.txt.gz"
##          
##          chrom_num = 21
##          split_name = unlist(strsplit(name_template, split="%"))
##          filestring = paste0(datadir, "/", split_name[1], chrom_num, split_name[2])
##            
##          d = fread(filestring,
##                    na.strings = ('NA'),
##                    colClasses = c(RSID = 'character',
##                                   CHR = 'integer',
##                                   POS = 'integer',
##                                   cM_Block = 'integer'),
##                    data.table=FALSE,
##                    showProgress=TRUE)
##          
##          tmp = names(d)
##          tmp2 = tmp[grep("gnomad_*", tmp)]
##
##          excall_block_bootstrap(est_ancestries = tmp2)
##      ```
##  
excall_block_bootstrap <- function(est_ancestries = c("gnomad_AF_afr", "gnomad_AF_amr", "gnomad_AF_asj",
                                                      "gnomad_AF_eas", "gnomad_AF_fin", "gnomad_AF_nfe", 
                                                      "gnomad_AF_sas", "gnomad_AF_oth"),
                                   datadir = "/referencemerge_with_pel/exomes_with_bherer/",
                                   name_template = "exome_cmInfo_chrom%_output.txt.gz",
                                   outdir = "/exomes_all/",
                                   nworkers = 4)
{
  library(R.utils)
  library(data.table)
  library(tidyverse)
  library(progress)
  
  ################################
  ## Build rbind of all chromosomes
  ################################
  print("Building stacked dataset")
  
  bherer = NULL;
  
  pb <- progress_bar$new(total = 22)
  for (chrom_num in c(1:22)) {
    pb$tick()
    
    split_name = unlist(strsplit(name_template, split="%"))
    filestring = paste0(datadir, "/", split_name[1], chrom_num, split_name[2])
    
    d = fread(filestring,
              na.strings = ('NA'),
              colClasses = c(RSID = 'character',
                             CHR = 'integer',
                             POS = 'integer',
                             cM_Block = 'character'),
              data.table=FALSE,
              showProgress=TRUE)
    
    if (chrom_num == 1) {
      bherer = d
    }
    else {
      bherer = rbind(bherer, d)
    }
  }
  # Drop unused variable
  rm(d)
  
  ################################
  ## Call the Bootstrap 
  ################################
  sliced_frame = NULL
  for (anc in est_ancestries) {
    column_anc = bherer[,anc]
    if (! all(is.na(column_anc))) {
      sliced_frame = bherer[which(! is.na(column_anc)),
                            c("RSID", "cM_Block", "CHR", "POS", "ALT", "REF",
                              "ref_AF_eur_1000G", "ref_AF_afr_1000G", "ref_AF_sas_1000G",
                              "ref_AF_eas_1000G", "ref_AF_iam_1000G", anc)]
      
      retval = block_bootstrap(genome_data = sliced_frame,
                               ancestries = c("ref_AF_afr_1000G",
                                              "ref_AF_iam_1000G",
                                              "ref_AF_eur_1000G", 
                                              "ref_AF_eas_1000G",
                                              "ref_AF_sas_1000G"),
                               target_anc = anc,
                               lower_upper_quantiles = c(0.025, 0.975),
                               NRESAMPLES = 1000,
                               NWORKERS = nworkers) 
      sliced_frame = NULL
      
      out_filestring = paste0(outdir, "/bootstrap_CI-", anc, ".txt.gz")
      fwrite(retval[['estimate_info']], file=out_filestring, 
             na='NA', sep='\t',
             row.names=FALSE, compress='gzip', quote=FALSE)
      
      out_filestring = paste0(outdir, "/bootstrap_CI-", anc, "-replicate_distrib.txt.gz")
      fwrite(retval[['replicates']], file=out_filestring, 
             na='NA', sep='\t',
             row.names=FALSE, compress='gzip', quote=FALSE)
    }
  }
  
  return(TRUE)
}
