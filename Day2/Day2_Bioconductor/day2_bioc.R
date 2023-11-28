## ----setup, message=FALSE, warning=FALSE, echo=FALSE, eval=TRUE, include=FALSE--------
knitr::opts_chunk$set(eval=TRUE, echo=TRUE, fig.align="center", fig.show='hold', size='tiny', fig.path='figures/beamer-')

options(digits=4)


## ----eval=FALSE, size='tiny'----------------------------------------------------------
## if (!requireNamespace("BiocManager", quietly = TRUE))
##     install.packages("BiocManager")
## BiocManager::install(version = "3.18")


## ----eval=FALSE, size='tiny'----------------------------------------------------------
## BiocManager::install("minfi")


## ----eval=FALSE, size='tiny'----------------------------------------------------------
## #For large files, increase timeout from default e.g. 120sec
## options(timeout=120)
## BiocManager::install()


## ----eval=FALSE, size='tiny'----------------------------------------------------------
## #RNA-Seq workflow
## BiocManager::install("rnaseqGene")
## browseVignettes("rnaseqGene")


## ----eval=FALSE, size='tiny'----------------------------------------------------------
## browseVignettes("DESeq2")


## ----eval=TRUE, size='tiny'-----------------------------------------------------------
df <- data.frame(x = 1:10, y = letters[1:10])
class(df)
t(df) # t is the generic function
t.data.frame(df) # S3 method for data.frame objects


## ----eval=TRUE, size='tiny'-----------------------------------------------------------
methods(t) # what other methods for generic function t() ? 
methods(t.test) # bad name for a generic function


## ----eval=TRUE, size='tiny'-----------------------------------------------------------
methods(class="lm") # what methods to apply to lm objects?


## ----eval=FALSE, size='tiny'----------------------------------------------------------
## setClass("myNewClass",
## representation(x = "numeric", y = "numeric"))
## n1 <- new("myNewClass", x = rnorm(10), y = 1:10)
## n1[1]
## n1$x # Doesn't work!
## n1
## n1@x # @-operator is used to access slots


## ----eval=TRUE,message=FALSE,warning=FALSE, size='tiny'-------------------------------
# Base functions / Infrastructure for Bioconductor:
library(Biobase)
?ExpressionSet

#simulate 6 observations of 100 genes 
#by sampling from a uniform distribution
eset <- ExpressionSet(assayData=matrix(runif(600), 
                                       nrow=100, ncol=6))
eset


## ----eval=TRUE, size='tiny'-----------------------------------------------------------
exprs(eset)[1:5,]
#each slot can be accessed using the @ operator, 
#like the $ for S3 classes
eset@assayData


## ----eval=TRUE, size='tiny'-----------------------------------------------------------
#generate gender
gender <- rep(c("Male", "Female"), 3) 
#generate some categorical
type <- as.factor(rep(c("A", "B", "C"), 2)) 
pdata <- data.frame(gender, type)
#annotated data frame
pdata <- new("AnnotatedDataFrame", data = pdata) 
#add pdata to eset
phenoData(eset) <- pdata 
eset


## ----eval=TRUE, size='tiny'-----------------------------------------------------------
#each slot can be accessed using the @ operator, 
#like the $ for S3 classes
str(eset)


## ----eval=TRUE, size='tiny'-----------------------------------------------------------
# use @ access annotated data frame and $ for gender
eset@phenoData$gender 
#same but different; method to access phenoData
phenoData(eset)$gender 
# directly access gender; works only for phenoData
eset$gender 


## ----eval=TRUE, size='tiny'-----------------------------------------------------------
# assayData can be accessed/extracted by exprs()
head(exprs(eset)) 
# via @ only locked environment
eset@assayData 
# does not allow users
# to directly manipulate expression data


## ----eval=TRUE, size='tiny'-----------------------------------------------------------
eset[1:10, ] # 10 genes, phenoData stays the same


## ----eval=TRUE, size='tiny'-----------------------------------------------------------
eset <- eset[, order(eset$gender)] #reordering by gender
exprs(eset)[1:2,]  # both assayData 
eset$gender #and phenoData are reordered
eset$type


## ----eval=FALSE, size='tiny'----------------------------------------------------------
##   library(org.Hs.eg.db) #Genome wide annotation for Human


## ----eval=FALSE, size='tiny'----------------------------------------------------------
##   library(hgu133plus2.db) #Affymetrix Human Genome U133 Plus 2


## ----eval=FALSE, size='tiny'----------------------------------------------------------
##   library(GO.db) #Annotation describing the Gene Ontology


## ----eval=FALSE, size='tiny'----------------------------------------------------------
##   #TxDb objects map UTRs, CDS, exons and transcripts
##   library(TxDb.Hsapiens.UCSC.hg38.knownGene)


## ----eval=FALSE, size='tiny'----------------------------------------------------------
##   # Full genome sequences for UCSC version hg38 stored
##   # as Biostrings object
##   library(BSgenome.Hsapiens.UCSC.hg38)


## ----eval=FALSE, size='tiny'----------------------------------------------------------
##   library(biomaRt)


## ----eval=TRUE,warning=FALSE,message=FALSE, size='tiny'-------------------------------
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
ensids <- c("ENSG00000130720", "ENSG00000103257", "ENSG00000156414")
select(org.Hs.eg.db, keys=ensids, 
       columns=c('ENSEMBL','SYMBOL'), keytype="ENSEMBL")


## ----eval=TRUE,warning=FALSE,message=FALSE, size='tiny'-------------------------------
#retrieve annotation via http://www.biomart.org/
library(biomaRt) 
ensembl <- useDataset("hsapiens_gene_ensembl",
                      mart = useMart("ensembl")) 
affyids <- c("202763_at", "209310_s_at")
getBM(attributes = c("affy_hg_u133_plus_2", "hgnc_symbol",
                     "chromosome_name","start_position",
                     "end_position", "band"),
    filters = "affy_hg_u133_plus_2",
    values = affyids, mart = ensembl)


## ----eval=TRUE,warning=FALSE,message=FALSE, size='tiny'-------------------------------
#Representation and manipulation of genomic intervals
library(GenomicRanges)
library(Gviz)
data(cpgIslands)
class(cpgIslands)
cpgIslands[1:3,]


## ----eval=TRUE,warning=FALSE,message=FALSE, size='tiny'-------------------------------
#Representation and manipulation of genomic intervals
start(cpgIslands)
end(cpgIslands)
seqnames(cpgIslands) # the chromosomes for each ranges
seqlevels(cpgIslands) # the possible chromosomes
seqlengths(cpgIslands) # the lengths for each chromosome


## ----eval=TRUE,warning=FALSE,message=FALSE, size='tiny'-------------------------------
# simulate some methylation values
elementMetadata(cpgIslands) <- matrix(runif(50), nrow = 10) 
cpgIslands[1:3,]


## ----eval=TRUE,warning=FALSE,message=FALSE, size='tiny'-------------------------------
# again, plot elements are organized as tracks
chr <- seqlevels(cpgIslands) # only chromosome 7
gen <- genome(cpgIslands) # human genome build 19
gen
atrack <- AnnotationTrack(cpgIslands, name = "CpG") #annotation
dtrack <- DataTrack(cpgIslands, type = "p", jitter.x = T,
                    name = "beta",  ylim = c(0, 1.05),
                    legend = TRUE, cex = 1.2) 
dtrack # data track that shows the methylation values as points


## ----eval=TRUE,warning=FALSE,message=FALSE,out.width='70%', size='tiny'--------------
# axis that shows basepair position:
gtrack <- GenomeAxisTrack()
# ideogram :
itrack <- IdeogramTrack(genome = gen, chromosome = chr)
plotTracks(list(itrack, gtrack, dtrack, atrack))

