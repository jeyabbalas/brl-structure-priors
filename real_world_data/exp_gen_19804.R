# Installing core bioconductor packages
source("http://bioconductor.org/biocLite.R")
biocLite()

biocLite("GEOquery")
library(GEOquery)
# GEO accession ID
geo_id = "GSE19804"

# Download dataset
filePaths = getGEOSuppFiles(geo_id)
geo_dir = paste("./", paste(geo_id,"/",sep=""), sep="")
geo_tarfile = paste(geo_dir, paste(geo_id,"_RAW.tar",sep=""), sep="")
celdir = paste(geo_dir, paste(geo_id,"_RAW/",sep=""), sep="")
out_file = paste(celdir, paste(geo_id,"_exp.txt",sep=""), sep="")
matrix_file = paste(celdir, paste(geo_id,"_series_matrix.txt",sep=""), sep="")


# Untar the archived files
untar(geo_tarfile, exdir=celdir)



################
# Read CEL Files
################
library(affy)
setwd(celdir)

data.raw = ReadAffy()


# Additional information
slotNames(data.raw)
annotation_type <- annotation(data.raw)



####################
# RMA normalization 
####################
data.eset = rma(data.raw)
data.exp = exprs(data.eset)

# Remove suffix ".CEL.gz" from sample IDs.
sample_ids = colnames(data.exp)
sample_ids = sapply(strsplit(sample_ids, "\\."), "[", 1)
colnames(data.exp) <- sample_ids

# Format to 5 decimals
data.exp=format(data.exp, digits=5)



#############
# Annotation
#############
# Install the relevant .db annotation package eg. hgu133plus2.db
#Check
biocLite("hgu133plus2.db")
library(hgu133plus2.db)

probe_ids = rownames(data.exp)
gene_symbols = unlist(mget(probe_ids, hgu133plus2SYMBOL, ifnotfound=NA))

annotated = as.data.frame(cbind(probe_ids, gene_symbols))



################
# IQR filtering
################
# Merging annotation with expression data
data.exp.df <- as.data.frame(data.exp)
data.exp.df$probe_ids <- rownames(data.exp.df)
data.annotated = merge(data.exp.df, annotated, by.x="probe_ids", by.y="probe_ids")
write.table(data.annotated, paste(celdir,"exp_annotated.txt"),sep="\t",row.names=F)
data.annotated = read.delim(paste(celdir,"exp_annotated.txt"),sep="\t",check.names=F)


# Sorting by gene symbols
data.annotated.sorted = data.annotated[order(data.annotated$gene_symbols),]
logdata = data.annotated.sorted[,!(colnames(data.annotated.sorted) %in% c("probe_ids", "gene_symbols"))]
unlogdata = 2^logdata


# Calculating IQR for all probes using unlog data
iqr <- apply(unlogdata,1,IQR)
data.iqr = cbind(data.annotated.sorted[,(colnames(data.annotated.sorted) %in% c("probe_ids", "gene_symbols"))], iqr, unlogdata)
write.table((data.iqr), paste(celdir,"rma_unlog.IQR.txt"), sep="\t",row.names=F)


# Keep probe with highest iqr in case of multiple probes
names(iqr) = data.annotated.sorted$probe_ids
iqrs = split.default(iqr, data.annotated.sorted$gene_symbols)
maxes = sapply(iqrs, function(x) names(which.max(x)))
singleprobe = data.iqr[data.iqr$probe_ids %in% maxes, !(colnames(data.iqr) == "probe_ids")]


## remove row with gene symbol NA
newdata = singleprobe
write.table(newdata, paste(celdir,"singleprobe_unlogged.txt"),sep="\t", row.names=FALSE,quote=FALSE)

d = newdata[,!(colnames(newdata) %in% c("gene_symbols", "iqr"))]
gene_symbols <- newdata[,(colnames(newdata) %in% c("gene_symbols"))]
logd = cbind(gene_symbols, log2(d))
write.table(logd, paste(celdir,"singleprobe_logged.txt"),sep="\t", row.names=FALSE,quote=FALSE)




#############
# Phenotype
#############
getPhenotypeMap = function (matrix_file){
        library("GEOquery")
        gse = getGEO(filename=matrix_file)
        geo_acc_id = gse@phenoData@data[["geo_accession"]]
        #Check
        phenotype = gse@phenoData@data[["characteristics_ch1"]]
        names(phenotype) = geo_acc_id
        return(phenotype)
}

logd = read.delim(paste(celdir,"singleprobe_logged.txt"),sep="\t",check.names=F, stringsAsFactors=FALSE)

phenoMap <- getPhenotypeMap(matrix_file)

sample_ids = colnames(data.exp)
id_variable = character(length=length(sample_ids)+1)
target_variable = character(length=length(sample_ids)+1)
#BRL formatted headers
id_variable[1] = "#Sample"
target_variable[1] = "@Class"
count = 2
for(id in sample_ids){
        #id = unlist(strsplit(id, "_"))[1]
        id_variable[count] = id
        target_variable[count] = phenoMap[[id]]
        count = count+1
}

names(target_variable) = names(logd)
names(id_variable) = names(logd)
logd.target <- rbind(target_variable, logd)
logd.target <- rbind(id_variable, logd.target)

# Transpose
logd.target <- t(logd.target)

###################
# Print final data
###################
write.table(logd.target, file = out_file, quote = FALSE, sep = "\t", row.names=FALSE, col.names=FALSE)
