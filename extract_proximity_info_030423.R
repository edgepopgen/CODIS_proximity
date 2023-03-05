#February 2023
#code to extract positions of features near STR markers from the
#hipstr reference.

if(!("data.table" %in% installed.packages())){install.packages("data.table")}
library(data.table)


#General function to get distances from each str in a table to the nearest
#feature noted in another table (feat.tab), plus the number of features within a
#cutoff distance for a vector of such cutoffs. 
#str.tab is a table of strs, needs to have a column called mid with the midpoint coord,
#and chromosome needs to be in a variable called chr. chromosome labeling scheme must match
#the feature table exactly.
#feature tab is a data frame of features. chromosome info is in the (numeric) column feat.chrcol
#has either one position marker for a feature (feat.poscol)
#or, if the features are ranges, has two position coordinates.
#dist.cuts is a vector of distances (in basepairs) within which to count observed features
#(if an str midpoint is inside a feature, i.e. between the two position coords,
#then the distance to that feature is 0, otherwise it is the distance to the closest end.
#returns a matrix with 1 row per str. first column is the distance to the nearest
#feature; remaining columns are #s of features within the distances specified in dist.cuts.
nearby <- function(str.tab, feat.tab, feat.chrcol, feat.poscol, feat.poscol2 = NULL, 
                   dist.cuts = c(1000,10000,100000), feat.aggcol = NULL, aggfun = sum){
  nstr <- dim(str.tab)[1]
  result <- matrix(nrow = nstr, ncol = 1 + length(dist.cuts) + length(dist.cuts) * !is.null(feat.aggcol) )
  chr = "nonchr"
  for(i in 1:nstr){
    if(str.tab$chr[i] != chr){
      chr <- str.tab$chr[i]
      feat.chr <- feat.tab[feat.tab[[feat.chrcol]] == str.tab$chr[i],]
    }
    if(is.null(feat.poscol2)){
      dists <- abs(str.tab$mid[i] - feat.chr[[feat.poscol]])
    }
    if(!is.null(feat.poscol2)){
      dists.1 <- str.tab$mid[i] - feat.chr[[feat.poscol]]
      dists.2 <- str.tab$mid[i] - feat.chr[[feat.poscol2]]
      inside <- dists.1 * dists.2 <= 0
      dists <- apply(cbind(abs(dists.1), abs(dists.2)), 1, min)
      dists[inside == 1] <- 0
    }
    result[i,1] <- min(dists)
    for(j in 1:length(dist.cuts)){
      result[i,1+j] <- sum(dists <= dist.cuts[j])
      if(!is.null(feat.aggcol)){
        result[i, 1 + length(dist.cuts) + j] <- aggfun(feat.chr[dists <= dist.cuts[j],][order(dists[dists <= dist.cuts[j]]),][[feat.aggcol]]) 
      }
    }
    
    if(i %% 10000 == 0){print(paste(i,"STRs completed"))}
  }
  result
}

#Read in STRs

str.loc <- fread("hg19.hipstr_reference.bed")
names(str.loc) <- c("chr", "start", "end", "motif.len", "V5", "name", "motif")

#List of CODIS STR names
STR_NAMES=data.frame("D22S1045","TPOX","D2S441","D2S1338","vWA","D12S391",
                     "D5S818","CSF1PO","D1S1656","D10S1248","TH01","D13S317",
                     "D16S539","D18S51","D19S433","D21S11","D3S1358","FGA",
                     "D7S820","D8S1179")

str.loc$CODIS <- as.numeric(str.loc$name %in% STR_NAMES)
str.loc$mid <- (str.loc$start + str.loc$end)/2



###################################################################
###################################################################
#RefSeq Select
refseq.sel <- fread("20230120_data_integrator_RefSeqSelect")

refseq.sel$TSS <- refseq.sel$ncbiRefSeqSelect.txStart
refseq.sel$TSS[refseq.sel$ncbiRefSeqSelect.strand == "-"] <- refseq.sel$ncbiRefSeqSelect.txEnd[refseq.sel$ncbiRefSeqSelect.strand == "-"]
dim(refseq.sel) #21,432 protein coding genes

refseq.sel.TSS.dists <- nearby(str.loc, refseq.sel, 2, 16)
refseq.sel.TSS.dists.dt <- as.data.table(refseq.sel.TSS.dists) 
setnames(refseq.sel.TSS.dists.dt, new = c("refseq.sel.TSS.nearest", "refseq.sel.TSS.win1kb", 
                                             "refseq.sel.TSS.win10kb", "refseq.sel.TSS.win100kb"))

refseq.sel.gene.dists <- nearby(str.loc, refseq.sel, 2, 4, 5)
refseq.sel.gene.dists.dt <- as.data.table(refseq.sel.gene.dists) 
setnames(refseq.sel.gene.dists.dt, new = c("refseq.sel.gene.nearest", "refseq.sel.gene.win1kb", 
                                             "refseq.sel.gene.win10kb", "refseq.sel.gene.win100kb"))

rm(refseq.sel.TSS.dists)
rm(refseq.sel.gene.dists)
rm(refseq.sel)


###################################################################
###################################################################
#Hapmap SNPs at high freq in CEU/CEPH

hapmap <- fread("20230120_data_integrator_HapMapCEU") # 4029798 SNPs

hapmap$n <- hapmap$hapmapSnpsCEU.homoCount1 + hapmap$hapmapSnpsCEU.heteroCount + hapmap$hapmapSnpsCEU.homoCount2
hapmap$af1 <- (hapmap$hapmapSnpsCEU.homoCount1 + hapmap$hapmapSnpsCEU.heteroCount/2) / hapmap$n
#summary(hapmap$af1)
#mean(hapmap.af1 >= 0.01 & hapmap.af1 <= 0.99)
hapmap <- hapmap[hapmap$af1 >= 0.01 & hapmap$af1 <= 0.99,] #reduce to 2705918 sites

hapmap.dists <- nearby(str.loc, hapmap, 1, 2)
hapmap.dists.dt <- as.data.table(hapmap.dists)
setnames(hapmap.dists.dt, new = c("hapmapCEU.nearest",
                  "hapmapCEU.win1kb","hapmapCEU.win10kb","hapmapCEU.win100kb"))

rm(hapmap.dists)
rm(hapmap)

###################################################################
###################################################################
#Clinvar pathogenic variants

clinvar <- fread("20230120_data_integrator_ClinVar") 
clinvar <- clinvar[clinvar$clinvarMain.clinSign == "Pathogenic"]

clinvar.dists <- nearby(str.loc, clinvar, 1, 2)
clinvar.dists.dt <- as.data.table(clinvar.dists)
setnames(clinvar.dists.dt, new = c("clinvarpath.nearest", "clinvarpath.win1kb",
                                   "clinvarpath.win10kb", "clinvarpath.win100kb")) 
rm(clinvar.dists)
rm(clinvar)

###################################################################
###################################################################
#GWAS catalog - filtered to only include unique rsids
#distance to nearest gwas hit
#no. gwas w/in 1kb, 10kb, 100kb
#no. traits w/in 1kb, 10kb, 100kb

gwas.cat <- fread("20230120_data_integrator_GWAScat") #392271 rows

#183014 unique rsids
gwas.agg <- aggregate(gwas.cat, list(gwas.cat$gwasCatalog.name), function(x){x[1]})
#aggregation decreases to 183014 rows

rm(gwas.cat)

gwas.allhits.dists <- nearby(str.loc,gwas.agg,2,3)
gwas.allhits.dists.dt <- as.data.table(gwas.allhits.dists)
setnames(gwas.allhits.dists.dt, new = c("gwas.allhits.nearest", "gwas.allhits.win1kb", 
                                        "gwas.allhits.win10kb", "gwas.allhits.win100kb"))

rm(gwas.allhits.dists)
rm(gwas.agg)


###################################################################
###################################################################
#GWAS catalog - filtered to commonly studied traits
#no. traits w/in 1kb, 10kb, 100kb

gwas.cat <- fread("20230120_data_integrator_GWAScat") #392271 rows


studies.per.trait <- aggregate(gwas.cat$gwasCatalog.pubMedID, by = list(gwas.cat$gwasCatalog.trait), function(x){length(unique(x))})
traits.threeplus <- studies.per.trait$Group.1[studies.per.trait$x >= 3]
dim(gwas.cat[gwas.cat$gwasCatalog.trait %in% traits.threeplus,]) # 206640 rows remain, 493 traits
#traits.fourplus <- studies.per.trait$Group.1[studies.per.trait$x >= 4]
#dim(gwas.cat[gwas.cat$gwasCatalog.trait %in% traits.fourplus,])

gwas.cat <- gwas.cat[gwas.cat$gwasCatalog.trait %in% traits.threeplus,]

gwas.agg <- aggregate(gwas.cat, list(gwas.cat$gwasCatalog.name, gwas.cat$gwasCatalog.trait), function(x){x[1]})
#aggregation decreases to 146039 rows

rm(studies.per.trait)
rm(gwas.cat)

gwas.dists <- nearby(str.loc,gwas.agg,3,4, feat.aggcol = 12, aggfun = function(x){length(unique(x))})
gwas.trait.dists.dt <- as.data.table(gwas.dists[,5:7])
setnames(gwas.trait.dists.dt, new = c("gwastraits.1kb", "gwastraits.10kb", "gwastraits.100kb"))
rm(traits.threeplus)
rm(gwas.dists)
rm(gwas.agg)



###################################################################
###################################################################
#DNase clusters w/in 1, 10, 100

DNase.loc <- fread("20230120_data_integrator_ENCODEregulation_DNaseClustersV3")
#1949038 rows

DNase.hs <- DNase.loc[DNase.loc$wgEncodeRegDnaseClustered.score == 1000,]
rm(DNase.loc)

DNase.dists <- nearby(str.loc,DNase.hs,1,2)
DNase.dists.dt <- as.data.table(DNase.dists)
setnames(DNase.dists.dt, new = c("DNase.nearest", 
          "DNase.win1kb", "DNase.win10kb", "DNase.win100kb"))

rm(DNase.dists)
rm(DNase.hs)

                                   
combined.dt <- cbind(str.loc, refseq.sel.TSS.dists.dt, refseq.sel.gene.dists.dt, 
                                 hapmap.dists.dt, clinvar.dists.dt, 
                     gwas.allhits.dists.dt, gwas.trait.dists.dt, DNase.dists.dt)
fwrite(combined.dt, file = "STR_proximity.csv")

################################################################
#CODIS ONLY

###################################################################
###################################################################
#dbSNP common SNPs, only for CODIS because it's slow to do on everything

snp.loc <- fread("20230120_data_integrator_dbSNP153Common")
dim(snp.loc) #14908376 common variants
dbsnp.codis.dists <- nearby(str.loc[str.loc$CODIS == 1,], snp.loc, 1, 2)
rm(snp.loc)
dbsnp.codis.dists.dt <- as.data.table(dbsnp.codis.dists)
setnames(dbsnp.codis.dists.dt, new  = c("dbsnp153.common.nearest","dbsnp153.common.win1kb",
                                                  "dbsnp153.common.win10kb","dbsnp153.common.win100kb"))
rm(dbsnp.codis.dists)


###################################################################
###################################################################
#revisit ClinVar, pull out specific names of phenotypes

clinvar <- fread("20230120_data_integrator_ClinVar") 
clinvar <- clinvar[clinvar$clinvarMain.clinSign == "Pathogenic"]

codis.loc <- str.loc[str.loc$CODIS == 1,]
clinvar[clinvar[[1]] == "chr1" & clinvar[[2]] >= 230905396 - 10000 & clinvar[[2]] <= 230905396 + 10000 ]

clinvar.dists.codis <- nearby(codis.loc, clinvar, 1, 2, feat.aggcol = 18, aggfun = function(x){paste(x, sep = "", collapse = "|")})

clinvar.phenmat <- matrix(nrow = dim(clinvar.dists.codis)[1], ncol = 3)
clinvar.nphenmat <- matrix(nrow = dim(clinvar.dists.codis)[1], ncol = 3)
for(i in 1:3){
  for(j in 1:dim(clinvar.dists.codis)[1]){
    phenvec <- strsplit(clinvar.dists.codis[j, 4+i], "|", fixed = TRUE)[[1]]
    phenvec <- phenvec[phenvec != "not provided"]
    phen.uni <- unique(phenvec)
    n.uni <- length(phen.uni)
    phen.uni.paste <- paste(phen.uni, sep = "", collapse = "|")
    clinvar.phenmat[j,i] <- phen.uni.paste
    clinvar.nphenmat[j,i] <- n.uni
  }
}

clinphen <- as.data.table(clinvar.phenmat)
setnames(clinphen, new = c("clinvarpath.phen1kb", "clinvarpath.phen10kb", "clinvarpath.phen100kb"))
clin.nphen <- as.data.table(clinvar.nphenmat)
setnames(clin.nphen, new = c("clinvarpath.nphen1kb", "clinvarpath.nphen10kb", "clinvarpath.nphen100kb"))

rm(clinvar)
rm(clinvar.dists.codis)
rm(clinvar.phenmat)
rm(clinvar.nphenmat)
rm(phenvec)
rm(n.uni)
rm(phen.uni.paste)


###################################################################
###################################################################
#revisit gwas catalog, pull out commonly-studied-phenotype names

gwas.cat <- fread("20230120_data_integrator_GWAScat") #392271 rows

studies.per.trait <- aggregate(gwas.cat$gwasCatalog.pubMedID, by = list(gwas.cat$gwasCatalog.trait), function(x){length(unique(x))})
traits.threeplus <- studies.per.trait$Group.1[studies.per.trait$x >= 3]
dim(gwas.cat[gwas.cat$gwasCatalog.trait %in% traits.threeplus,]) # 206640 rows remain, 493 traits

gwas.cat <- gwas.cat[gwas.cat$gwasCatalog.trait %in% traits.threeplus,]

gwas.agg <- aggregate(gwas.cat, list(gwas.cat$gwasCatalog.name, gwas.cat$gwasCatalog.trait), function(x){x[1]})
#aggregation decreases to 146039 rows

rm(studies.per.trait)
rm(gwas.cat)

gwas.dists.codis <- nearby(str.loc[str.loc$CODIS==1,],gwas.agg,3,4, feat.aggcol = 12, aggfun = function(x){paste(unique(x), sep = "", collapse = "|")})

gwasphen <- as.data.table(gwas.dists.codis[,5:7])
setnames(gwasphen, new = c("gwastraits.names.1kb", "gwastraits.names.10kb", "gwastraits.names.100kb"))
rm(gwas.dists.codis)
rm(gwas.agg)



###################################################################
###################################################################
#Revisit refseq select to get the names of nearby genes for CODIS
#RefSeq Select
refseq.sel <- fread("20230120_data_integrator_RefSeqSelect")


refseq.sel.gene.dists <- nearby(str.loc[str.loc$CODIS==1,], refseq.sel, 2, 4, 5, feat.aggcol = 12, aggfun = function(x){paste(x, sep = "", collapse = "|")})
refseq.sel.gd.dt <- as.data.table(refseq.sel.gene.dists[,5:7]) 
setnames(refseq.sel.gd.dt, new = c("genenames.1kb","genenames.10kb","genenames.100kb"))


rm(refseq.sel.gene.dists)
rm(refseq.sel)


###################################################################
###################################################################
#combine dbsnp and phenotype info with previously gathered info

all.str.prox <- fread("STR_proximity.csv")
codis.prox <- all.str.prox[all.str.prox$CODIS == 1,]
rm(all.str.prox)

codis.prox.supp <- cbind(codis.prox, dbsnp.codis.dists.dt, refseq.sel.gd.dt, clin.nphen, clinphen, gwasphen)
chr.num <- as.numeric(sub("...", "", codis.prox.supp$chr))

fwrite(codis.prox.supp[order(chr.num),], file = "CODIS_proximity.csv")

