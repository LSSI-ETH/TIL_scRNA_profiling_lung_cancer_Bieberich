#Read in GEX data:
bs833_noMHC.data <- Read10X(data.dir = "/Users/fbieberich/Documents/ETH_projects/TIL_lungtumor_analysis10x/cocult_BS833_aMHC/GEX/noMHC/filtered_feature_bc_matrix/")

#Read in TCR data:
bs833_d0_TCR <- filter.clonotypes(tcr_location = "/Users/fbieberich/Documents/ETH_projects/TIL_lungtumor_analysis10x/5prime_VDJ_seq_outs/BS833/outs/filtered_contig_annotations.csv", prefix = "", VDJinfo = c("TRA", "TRB", "TRB"))
