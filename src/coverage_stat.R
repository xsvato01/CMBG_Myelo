library(data.table)

# get arguments
args = commandArgs(trailingOnly=TRUE)

# load file
#file_name <- basename(args[2])
file <- data.table(read.table(args[1], sep="\t"))
#file <- file[,-c("V3","V6","V7","V8")]
colnames(file) <- c("Chrom", "Start", "End",  "Exon", "Base_number", "Coverage")

#file_name <- basename(args[1])
#path <- dirname(args[1])

# create table with coverage per exon statistics
#file_ID <- "CXCR4"
Exon_stat <- list()

Exon_stat <- file[, list(Coverage_max = max(Coverage), Coverage_min = min(Coverage), Coverage_mean = round(mean(Coverage), digit = 2),
                         Coverage_median = median(as.double(Coverage))), by=.(Chrom, Exon, Start, End)]
write.table(Exon_stat, file = args[2], quote = F, row.names = F, sep = "\t")
#print(Exon_stat)