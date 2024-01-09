rm(list = ls())
require(readr)
require(data.table)
require(stringr)
require(edgeR)
require(limma)
require(SmartSVA)
library(optparse)

options(stringsAsFactors = F)
option_list = list(
  make_option("--data", action = "store", default = "/z/Comp/lu_group/Members/yixuan/AGEde/zenodo_AGEde_pipeline/data", type = "character"),
  make_option("--ge", action = "store", default = NA, type = "character"),
  make_option("--outdir", action = "store", default = "out", type = "character"),
  make_option("--outname", action = "store", default = NA, type = "character")
)
opt = parse_args(OptionParser(option_list=option_list))


# Read gene expression matrix for all genes and all samples (revised HYFA output format)
cat("Reading GTEx gene expressions\n")
exps <- as.data.frame(fread(opt$ge))
t <- unique(exps$Tissue)
# Replace colnames (gene symbol after transpose) with ENSEMBL id
ref = fread(paste0(opt$data, '/GTEx_genereads/ref.txt.gz'), header=F)
colnames(exps)[2:(ncol(exps)-1)] = ref$V1[match(colnames(exps)[2:(ncol(exps)-1)], ref$V2)]
rm(ref)
gc()
# Get samples and exp_genes for each tissue, and save them
cat("Generating data for", t, "\n")
exp_genes <- read.delim(paste0(opt$data, "/GTEx_genefiltering/expressed_genes/", t, ".genes.tsv"), stringsAsFactors = FALSE, header = FALSE)[,1]
chrs <- read.delim(paste0(opt$data, "/GTEx_genefiltering/chrs/", t, ".chrs.tsv"), stringsAsFactors = FALSE, header = FALSE)[,1]
exp_genes <- exp_genes[!chrs %in% 'chrY']
rownames(exps) <- exps[, 1]
exps <- exps[, -c(1, 2)]
ovp <- intersect(exp_genes, colnames(exps))
exps <- as.data.frame(t(exps[, match(ovp, colnames(exps))]))
cat(paste0("Size: ", nrow(exps), " genes and ", ncol(exps), " samples\n"))
# Write matrix with exps for the tissue
fwrite(exps, paste0(opt$outdir, "/aux/", t, ".", opt$outname, ".v8.gene_exps.txt.gz"), sep='\t', quote=F, row.names=T, col.names=T, na=NA)


# Calculate surrogate variables
cat("Calculating SVs for", t, "\n")
covs_imp <- read.table(paste0(opt$data, "/sample_covariates_updated/covs.basalcovs.imp_all.v8.txt"), row.names=1, check.names=F, header=T)
covs <- read.table(paste0(opt$data, "/sample_covariates_updated/covs.basalcovs.", t, ".v8.txt"), check.names=F, row.names=1, header=T)
covs <- rbind(covs, covs_imp[!(rownames(covs_imp) %in% rownames(covs)),])
covs <- covs[match(colnames(exps), rownames(covs)),]
# Prepare metadata. Missing values are imputed as the median.
covs$AGE <- as.numeric(covs$AGE %in% 'F')+1
covs$SMTSISCH[is.na(covs$SMTSISCH)] <- median(covs$SMTSISCH, na.rm=T)
covs$SMRIN[is.na(covs$SMRIN)] <- median(covs$SMRIN, na.rm=T)
if (opt$outname != 'lab_y') {
  fwrite(covs, paste0(opt$outdir, "/", t, ".", ifelse(opt$outname=='lab_yhat', 'lab_x', 'unlab_x'), ".covs.txt"), sep='\t', quote=F, row.names=T, col.names=T, na=NA)
}
# Models to evaluate (form is the full model, form0 is the null one)
form <- "~ AGE+SMTSISCH+SMRIN+AGE"
form0 <- "~ SMTSISCH+SMRIN+AGE"
if (sum(as.numeric(!is.na(covs$SMTSISCH))) < 1){
  covs$SMTSISCH <- NULL
  form <- "~ AGE+SMRIN+AGE"
  form0 <- "~ SMRIN+AGE"
  cat("Ischemic time not available. Drop variable. \n")
}

allSv <- fread(paste0(opt$outdir, "/aux/covs.svs.", t, ".", opt$outname, ".v8.FINAL.txt"))