#This script draws comparative histograms of genome-wide STR features vs CODIS 

if(!("dplyr" %in% installed.packages())){install.packages("dplyr")}
if(!("ggplot2" %in% installed.packages())){install.packages("ggplot2")}
if(!("patchwork" %in% installed.packages())){install.packages("patchwork")}
library(dplyr)
library(ggplot2)
library(patchwork)

#read in data, filter for autosomal tetramers, and add pseudocount
all_STRs = read.csv("STR_proximity.csv")
comp_STRs = filter(all_STRs, motif.len==4 & chr!='chrX' & chr!='chrY')
comp_STRs[,10:36] = comp_STRs[,10:36]+1

#pull out features we'll plot (refseq gene, clinvar path, GWAS hit, GWAS trait, DNAse sites; consider adding TSSs)
comp_STRs = select(comp_STRs, CODIS, refseq.sel.gene.nearest, refseq.sel.gene.win1kb, refseq.sel.gene.win10kb, refseq.sel.gene.win100kb, 
                   clinvarpath.nearest, clinvarpath.win1kb, clinvarpath.win10kb, clinvarpath.win100kb, 
                   gwas.allhits.nearest, gwas.allhits.win1kb, gwas.allhits.win10kb, gwas.allhits.win100kb,
                   gwastraits.1kb, gwastraits.10kb, gwastraits.100kb,
                   DNase.nearest, DNase.win1kb, DNase.win10kb, DNase.win100kb)
features = colnames(comp_STRs)[2:20]

#set up plot specific parameters and indices 
xlabs = c("distance to nearest gene", "genes w/in 1kb", "genes w/in 10kb", "genes w/in 100kb", 
          "distance to ClinVar path", "ClinVar path w/in 1kb", "ClinVar path w/in 10kb", "ClinVar path w/in 100kb",
          "distance to GWAS hit", "GWAS hits w/in 1kb", "GWAS hits w/in 10kb", "GWAS hits w/in 100kb", 
          "GWAS traits w/in 1kb", "GWAS traits w/in 10kb", "GWAS traits w/in 100kb",
          "distance to DNAse site", "DNAse sites w/in 1kb", "DNAse sites w/in 10kb", "DNAse sites w/in 100kb")
row_startis = c(1, 5, 9, 13, 16)
nearest_is = c(1, 5, 9, 16)
x100 = c(2, 3, 4, 6, 11, 14, 15, 19)
x1000 = c(7, 8, 12)
x30 = c(10, 13, 18)
x5 = c(17)
  
comp_STRs$CODIS[comp_STRs$CODIS==0] = "non-CODIS"
comp_STRs$CODIS[comp_STRs$CODIS==1] = "CODIS"

plots = vector(mode='list', length=length(features))
for(i in 1:length(features)) {
  plots[[i]] <- local({
    i <- i
    curplot <- ggplot(comp_STRs, aes_string(x=features[i], fill="CODIS")) +
      geom_histogram(bins=30, aes(y=..density..), alpha=0.6, position = 'identity') +
      scale_fill_manual(values=c("#d95f02", "dimgray")) +
      theme_classic() +
      scale_x_log10() + #breaks=c(1, 11, 101), labels=c("0", "10", "100")) + 
      theme(legend.position = "none", 
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank()) +
      ylab("") +
      xlab(xlabs[i])

    #fix x-axis to account for pseudocount
    if(sum(i==x100)>0) {
      curplot = curplot + scale_x_log10(breaks=c(1, 11, 101), labels=c("0", "10", "100"))
    } else if(sum(i==x1000)>0) {
      curplot = curplot + scale_x_log10(breaks=c(1, 11, 101, 1001), labels=c("0", "10", "100", "1000"))
    } else if(sum(i==x30)>0) {
      curplot = curplot + scale_x_log10(breaks=c(1, 4, 11, 31), labels=c("0", "3", "10", "30"))
    } else if(sum(i==x5)>0) {
      curplot = curplot + scale_x_log10(breaks=c(1, 3, 6), labels=c("0", "2", "5"))
    } else if(sum(i==nearest_is)>0) {
      curplot = curplot + scale_x_log10(breaks=c(1, 101, 10001, 1000001), labels=c("0", "1e2", "1e4", "1e6"))
    }
    
    #add in legend once  
    if(i==4) {
      curplot <- curplot + theme(legend.position = c(.8,.8))
    }
    
    #add in y axis labels to row starts
    if(sum(row_startis==i)>0) {
      curplot <- curplot + ylab("Density") 
    }
    
    print(curplot)
  })
}

combined_plot <- 
  plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + 
  plots[[5]] + plots[[6]] + plots[[7]] + plots[[8]] + 
  plots[[9]] + plots[[10]] + plots[[11]] + plots[[12]] + 
  plot_spacer() + plots[[13]] + plots[[14]] + plots[[15]] + 
  plots[[16]] + plots[[17]] + plots[[18]] + plots[[19]] + plot_layout(ncol=4)

pdf("Figure1.pdf", width=12, height=15)
combined_plot
dev.off()
