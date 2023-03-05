# CODIS_proximity

These scripts were written for data processing and analysis for the manuscript "Microsatellites used in forensics are located in regions unusually rich in trait-associated variants" by Link, Zavaleta, Reyes, Ding, Wang, Rohlfs, & Edge.

They take as input files downloaded from UCSC genome browser's Data Integrator tool, as described in the paper, along with the hipstr human hg19 reference (https://github.com/HipSTR-Tool/HipSTR-references).

The script extract_proximity_info_030423.R extracts information about the number of features of various types near the 1.6 million STRs in the hipstr reference. It produces two files as output, one csv with the results for all 1.6 million STRs, and a second with extended results just for the CODIS markers, which is posted here.

draw_Figure1_030423.R, as one might expect, draws Figure 1 from the paper.

Figure2_tables3and4_030423.R assembles 10,000 random sets of 20 STRs and compares CODIS with them. It draws Figure 2, and with a described modification to a line near the beginning, draws Figure S1. It also assembles the information appearing in Tables 3-4 and S1-S3.
