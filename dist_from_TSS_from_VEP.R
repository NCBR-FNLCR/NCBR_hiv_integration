library(GenomicRanges)

dist2TSS <- function(bedlikeFile) {
    gtf <- rtracklayer::import("/data/CCBR_Pipeliner/db/PipeDB/Indices/hg19_basic/genes_only.gtf")
    gtf$gene_id <- matrix(unlist(strsplit(gtf$gene_id, split="\\.")),ncol=2,byrow=T)[,1]
    gtf <- sort(gtf,ignore.strand=T)

    bedlikeIn <- read.table(bedlikeFile,header=T,comment.char="")
    names(bedlikeIn)[which(names(bedlikeIn) == "STRAND")] <- "GENE_STRAND"
    bedlike <- makeGRangesFromDataFrame(bedlikeIn, keep.extra.columns=T)
    bedlike <- sort(bedlike,ignore.strand=T)

    bedlikeInGene <- bedlike[which(bedlike$Gene != "-")]
    bedlikeNoGene <- bedlike[which(bedlike$Gene == "-")]

    IDs <- numeric(length(bedlikeInGene))
    for (i in 1:length(bedlikeInGene)) {
    	IDs[i] <- which(gtf$gene_id == bedlikeInGene$Gene[i])
    }

    gtf1 <- gtf[IDs]
    bedlikeInGene$Dist_to_TSS <- distance(flank(bedlikeInGene,0), flank(gtf1,0), ignore.strand=T)

    gtf2 <- gtf[nearest(flank(bedlikeNoGene,0),flank(gtf,0),ignore.strand=T)]
    bedlikeNoGene$Dist_to_TSS <- distance(flank(bedlikeNoGene,0), flank(gtf2,0), ignore.strand=T)
    bedlikeNoGene$SYMBOL <- gtf2$gene_name
    bedlikeNoGene$Gene <- gtf2$gene_id
    bedlikeNoGene$BIOTYPE <- gtf2$gene_type

    allbedlike <- c(bedlikeInGene, bedlikeNoGene)
    allbedlike <- sort(allbedlike, ignore.strand=T)

    return(data.frame(allbedlike))
}


bedlikeFile <- "merged_sorted_insertions_vep_ATAC.bed"
outFile <- "merged_sorted_insertions_vep_TSS.bed"

allbedlike <- dist2TSS(bedlikeFile)

write.table(allbedlike, outFile, quote=F, sep="\t", row.names=F) 