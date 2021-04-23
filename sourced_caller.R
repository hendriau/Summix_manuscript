library(R.utils)
library(data.table)

source("/cM_block-estimates/block_bootstrap.R")

## For exome
ddir = "/referencemerge_with_pel/exomes_with_bherer/"
template = "exome_cmInfo_chrom%_output.txt.gz"
odir = "/exomes_all/"

# Loads in correct header from chromosome 21, full bootstrap run on all autosomes
chrom_num = 21
split_name = unlist(strsplit(template, split="%"))
filestring = paste0(ddir, "/", split_name[1], chrom_num, split_name[2])
  
d = fread(filestring,
          na.strings = ('NA'),
          colClasses = c(RSID = 'character',
                         CHR = 'integer',
                         POS = 'integer',
                         cM_Block = 'integer'),
          data.table=FALSE,
          showProgress=TRUE)

tmp = names(d)
torun = tmp[grep("_(afr|amr|nfe|eas|sas|oth)$", tmp)]

d = NULL;

ran = FALSE
ran = excall_block_bootstrap(est_ancestries = 'gnomad_AF_afr',
                             datadir = ddir,
                             name_template = template,
                             outdir = odir,
                             nworkers = 16)

print(ran)
