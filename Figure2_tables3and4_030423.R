#February 2023
#script to compare CODIS with non-CODIS markers, given that features have
#already been extracted by extract_proximity_info.R. In particular, this script
#compares the means for extracted variables for CODIS as a set compared with the
#means of random sets of comparison markers.

if(!("data.table" %in% installed.packages())){install.packages("data.table")}
library(data.table)

strs <- fread("STR_proximity.csv")

#Define the subset of the 1.6 million STRs for comparison
#in the main text, this is non-codis autosomal tetramers
#To get table S1 / figure S1, use
#compare <- strs[strs[[1]] != "chrX" & strs[[1]] != "chrY" & strs$CODIS == 0 & strs$motif.len %in% 1:6,]
#To get the main text, use
compare <- strs[strs[[1]] != "chrX" & strs[[1]] != "chrY" & strs$CODIS == 0 & strs$motif.len == 4,]

#Just the CODIS markers
codis <- strs[strs$CODIS == 1,]

codis.means <- colMeans(codis[,10:36, with = FALSE])

#Select 10,000 random sets of 20 comparison loci and compute their means.
set.seed(8675309)
n.comp <- 10000
compare.means.dist <- matrix(nrow = n.comp, ncol = 27)
for(i in 1:n.comp){
  compare.means.dist[i,] <- colMeans(compare[sample(1:nrow(compare))[1:nrow(codis)], 10:36, with = FALSE])
}

#Compute the percentiles at which the CODIS means fall compared with
#randomly selected sets of STRs
percentiles <- numeric(27)
for(i in 1:27){
  percentiles[i] <- mean(compare.means.dist[,i] <= codis.means[i])
}
as.data.frame(percentiles, row.names = names(codis)[10:36])


xlabs <- c("distance to nearest gene (Mb)", "genes w/in 1kb", "genes w/in 10kb", "genes w/in 100kb",
           "distance to ClinVar path (Mb)","ClinVar path w/in 1kb","ClinVar path w/in 10kb","ClinVar path w/in 100kb",
           "distance to GWAS hit (Mb)", "GWAS hits w/in 1kb","GWAS w/in 10kb","GWAS w/in 100kb", 
           "GWAS traits w/in 1kb", "GWAS traits w/in 10kb", "GWAS traits w/in 100kb",
           "distance to DNase site (Mb)", "DNAse sites w/in 1kb", "DNase sites w/in 10kb", "DNase sites w/in 100kb")
inds <-  c(5:8, 13:27)

pal <- c('#1b9e77','#d95f02','#7570b3')


pdf("Figure2.pdf", width = 8, height = 10)
par(mfrow = c(5,4), mar = c(3.5,3.5,0.5,0.5), las = 1, mgp = c(2.2,0.5,0))
for(i in 1:length(inds)){
  toplot <- compare.means.dist[,inds[i]]
  if(i %in% c(1, 5, 9, 16)){toplot <- toplot / 1e6}
  myhist <- hist(toplot, main = "", xlab = xlabs[i], breaks = 50, ylab = "", yaxt = "n", col = "darkgray", lty = "blank")
  axis(2, at = c(0, max(myhist$counts)), labels = c("0", as.character(round(max(myhist$counts) ))))
  if(i %in% c(1, 5, 9, 13, 16)){title(ylab = "Frequency")}
  cod.line <- codis.means[inds[i]]
  if(i %in% c(1, 5, 9, 16)){cod.line <- cod.line / 1e6}
  lines(rep(cod.line,2), c(0, 10^6), lty = 2, lwd = 2, col = pal[2])
  if(i == 12){plot(compare.means.dist[,inds[i]], type = "n", bty = "n", xaxt = "n", yaxt = "n", ylab = "", xlab = "")}
}
dev.off()


#Define a new variable encoding intragenic status vs. not
mean(compare$refseq.sel.gene.nearest == 0) #39% intragenic for tetramers
ig <- compare$refseq.sel.gene.nearest == 0

#table 4
cor(cbind(ig, compare[,c(20, 12, 16, 24, 28, 31, 35)]), method = "spearman")

#table S2 analogous to Table 4 but at 100kb
cor(cbind(ig, compare[,c(21, 13, 17, 25, 29, 32, 36)]), method = "spearman")


###########################################################################
###########################################################################
#Repeat the random sets, matching on intragenic fraction (1/2) for supplement
compare.ig <- compare[ig==1,]
compare.notig <- compare[ig == 0,]
set.seed(8675309)
n.comp <- 10000
compare.means.dist.igmatched <- matrix(nrow = n.comp, ncol = 27)
for(i in 1:n.comp){
  c.ig <- compare.ig[sample(1:nrow(compare.ig))[1:10],]
  c.notig <- compare.notig[sample(1:nrow(compare.notig))[1:10],]
  compare.means.dist.igmatched[i,] <- colMeans(rbind(c.ig, c.notig)[, 10:36, with = FALSE])
}

percentiles.igmatched <- numeric(27)
for(i in 1:27){
  percentiles.igmatched [i] <- mean(compare.means.dist.igmatched[,i] <= codis.means[i])
}
#Table S3
as.data.frame(percentiles.igmatched , row.names = names(codis)[10:36])


