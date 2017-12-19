# Prepare example data for analysis
#
# Description:
#   learn mediating effect of methylation in the relationship
# between tissue and expression


# import sample labels
samples <- read.delim("data/Sample_info_RNAseq_RIN.txt")
# expression
exprs <- read.table("data/human_chimp_orth_exp_methyl_7725_hum.txt", sep="")
# methylation data
methyl <- read.csv("data/chimp_human_orth_7725_avg_methyl_per_ts_gene.txt", sep="")
# import previous computed results
#df_tmp <- get(load("Limma_output_fit1.rda"))

# <---- expression data matching/cleaning
exprs_tmp <- exprs[,2:48]
#exprs_tmp <- exprs_tmp[rownames(exprs_tmp) %in% rownames(df_tmp),]  
exprs_subset <- exprs_tmp[,1:31]
colnames(exprs_subset) <- gsub(".x", "", colnames(exprs_subset))

# <---- methylation data matching/cleaning
methyl_tmp <- methyl[,-17]
methyl_tmp[,32] <- as.character(methyl_tmp[,32])

methyl_tmp <-  methyl_tmp[methyl_tmp[,32] %in% rownames(exprs_subset),]
methyl_subset <- methyl_tmp[,1:31]
rownames(methyl_subset) <- methyl_tmp[,32]

# check agreement between methylation data and expression data
all.equal(rownames(methyl_subset), rownames(exprs_subset))
all.equal(colnames(methyl_subset), colnames(exprs_subset))

# <---- sample data matching
samples_sub <- samples[samples$Sample_ID %in% colnames(exprs_subset),]


# Subsetting data to human samples heart vs. kidney pair. Running linear models on exprs ~ tissue + RIN and also on exprs|methylation ~ tissue + RIN.

# <---- subset data to the tissue pair
# compare human heart and kidney
index_sample <- which(samples_sub$Species == "human"& (samples_sub$Tissue == "heart"|samples_sub$Tissue == "kidney"))
exprs_pair <- exprs_subset[,index_sample]
methyl_pair <- methyl_subset[,index_sample]
species <- droplevels.factor(samples$Species[index_sample])
tissue <- droplevels.factor(samples$Tissue[index_sample])
RIN <- samples$RIN[index_sample]

dim(exprs_pair)
dim(methyl_pair)

df <- list(exprs_pair=exprs_pair, 
           methyl_pair=methyl_pair, 
           tissue=tissue, 
           RIN=RIN)

save(df, file = "data/example-kidney-v-heart-human.rda")
