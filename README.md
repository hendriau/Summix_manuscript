## Summix
R scripts, data, and results used for the Summix manuscript.

#### ancestry_adjusted_data
Ancestry adjusted AF results for gnomAD AFR/AMR AFs matched to 1000G AFR/PEL populations.

#### summix_data
GnomAD exome and genome data merged with 1000G reference data and filtered to MAF>1% for at least one reference group and gnomAD group. This data was used to estimate ancestry proportions and ancestry adjusted AF for the 2021 manuscript (https://www.cell.com/ajhg/fulltext/S0002-9297(21)00221-4). 

#### R Scripts
*block_bootstrap.R*: Functions which perform the block bootstrap analysis of exome/genome data.

*five_ancestry_sim_example.R*: Example script used for simulations, for the five ancestry simulations.

*simulation_parameters_generate.R*: Generates simulation parameters, removes duplicate replicates.

*source_called.R*: An example script which executes the block bootstrap function.
