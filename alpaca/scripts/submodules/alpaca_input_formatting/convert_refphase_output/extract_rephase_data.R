library(optparse)

# Create option list
option_list <- list(
  make_option(c("--refphase_rData"), type="character", 
              help="path to refphase .RData file"),
  make_option(c("--output_dir"), type="character", 
              help="output directory")
)

# Create the parser
parser <- OptionParser(
  description = "Extract tsv files from refphase output",
  option_list = option_list
)

# Parse the arguments
args <- parse_args(parser)

# Extract values
results_object <- args$refphase_rData
output_dir <- args$output_dir
print('Loading refphase results from .RData file')
print(results_object)
tmp = load(results_object)
refphase_results = get(tmp[1])
phased_segments = as.data.frame(refphase_results$phased_segs)
snps = as.data.frame(refphase_results$phased_snps)
purity_ploidy = refphase_results$sample_data
print('Writing phased segments, snps and purity ploidy to tsv files')
print(output_dir)
write.table(phased_segments, file = paste0(output_dir, "/phased_segs.tsv"), sep = "\t", row.names = FALSE)
write.table(snps, file = paste0(output_dir, "/phased_snps.tsv"), sep = "\t", row.names = FALSE)
write.table(purity_ploidy, file = paste0(output_dir, "/purity_ploidy.tsv"), sep = "\t", row.names = FALSE)