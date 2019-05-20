##########################
### DIF Project: prioritize co-occurring actionable genes using large scale clinical‐grade cancer genomic data ###
# @copyright2019 Ziming Zhao; Briana Kubik
# Ziming Zhao <ziming.gt@gmail.com> or <ziming.zhao@jax.org>; Briana Kubik <bck3@geneseo.edu>
##########################

#####################
### Packages used ###
#####################
library(dplyr)
library(data.table)

# Data: Genie 
# Pan-cancer analyses

###################
### Import data ###
###################
setwd("/Users/zhaoz/Dropbox/DIF_R_pipeline")
mut.dat = read.csv("GENIE_data_mutations.csv", sep = "\t", header = T, stringsAsFactors = FALSE) #mutation data
CNA.dat = read.table("data_CNA_3.0.0.txt", header = T, stringsAsFactors = FALSE) #CNA data
fus.dat = read.csv("GENIE_fusion_data.csv", sep = "\t", stringsAsFactors = FALSE) #fusion data
sam.dat = read.csv("GENIE_data_sample.csv", skip = 3, sep = "\t", header = T, stringsAsFactors = FALSE) #sample data
pat.dat = read.csv("GENIE_data_patients.csv", skip = 3, header = T, stringsAsFactors = FALSE) #patient dat
CKB.dat = read.csv("CKBEfficacyForAnalysis.csv") #CKB clinical efficacy data
CKB.dat = CKB.dat[1:10] #we only want the first 10 columns

#########################
### Filter GENIE data ###
#########################
"
Remove unnecessary rows
  > Summary information (percent of each variant)
  > Remove rows whose variant column reads 'Silent', 'Intron', '3'UTR', '5'UTR', '3'Flank', or '5'Flank'
  > 8.211934% of total are 'unnecessary' for analysis and are removed
"
unique(mut.dat$Variant_Classification) #(16) unique variants before filtering

#Function to count variant types and returns variant types and the corresponding percentage given unique 
variant.percent.func = function(var.vec, data.col) { 
  hold.vec = c()
  for(i in 1:length(var.vec)) {
    hold.vec[i] = length(which(data.col == var.vec[i]))/length(data.col)*100
  }
  return(data.frame(variant = var.vec, percent = hold.vec))
}
VariantTypesPercentage=variant.percent.func(unique(mut.dat$Variant_Classification), mut.dat$Variant_Classification)

#Summarize the number and percentage of filtred variants: 21646
l1 = length(which(mut.dat$Variant_Classification == "Silent" | mut.dat$Variant_Classification == "Intron" |
                    mut.dat$Variant_Classification == "3'UTR" | mut.dat$Variant_Classification == "5'UTR" |
                    mut.dat$Variant_Classification == "3'Flank" | mut.dat$Variant_Classification == "5'Flank")) 
w1 = which(mut.dat$Variant_Classification == "Silent" | mut.dat$Variant_Classification == "Intron" |
             mut.dat$Variant_Classification == "3'UTR" | mut.dat$Variant_Classification == "5'UTR" |
             mut.dat$Variant_Classification == "3'Flank" | mut.dat$Variant_Classification == "5'Flank")

Percentage_filteredVariants = l1/length(mut.dat$Variant_Classification)*100 #Filtered variants: 8.211934% of total
Percentage_keptVariants = nrow(mut.dat[-w1,])/length(mut.dat$Variant_Classification)*100 #Kept variants: 91.78807 of total
# Pie Chart for filtered vs kept variants
par(mar = c(2,6,4,6) + 0.1)
slices2 = c(Percentage_keptVariants, Percentage_filteredVariants) 
lbls2 = c("KeptVariants","RemovedVariants")
pct2 = round(slices2/sum(slices2)*100)
lbls2 = paste(lbls2, pct2) # add percents to labels 
lbls2 = paste(lbls2, "%", sep = "") # add % to labels 
# Save the pie chart
png(file = "Percentage_Variants_PieChart.png")
pie(slices2, labels = lbls2, col = rainbow(length(lbls2)), cex = 1)
dev.off()

#PIE CHARTS for number of variants within each variant type

par(mar = c(2,6,4,6) + 0.1)
FrameIndels=length(which(mut.dat$Variant_Classification == "In_Frame_Del"| mut.dat$Variant_Classification == 'Frame_Shift_Del'|
                           mut.dat$Variant_Classification == 'Frame_Shift_Ins'| mut.dat$Variant_Classification =='In_Frame_Ins'))
Splice=length(which(mut.dat$Variant_Classification == "Splice_Site"| mut.dat$Variant_Classification == 'Splice_Region'))
TSS=length(which(mut.dat$Variant_Classification == 'Translation_Start_Site'))
NonSilent=length(which(mut.dat$Variant_Classification == "Nonsense_Mutation"| mut.dat$Variant_Classification == 'Missense_Mutation'|
                         mut.dat$Variant_Classification == 'Nonstop_Mutation'))
Silent=length(which(mut.dat$Variant_Classification == "Silent"))
Intron=length(which(mut.dat$Variant_Classification == "Intron"))

slices = c(FrameIndels,Splice, TSS, Intron,Silent, NonSilent) 
lbls = c("FrameIndels", "Splice",  "Translation_Start_Site","Intron",
         "Silent","NonSilent")
pct = round(slices/sum(slices)*100)
lbls = paste(lbls, pct) # add percents to labels 
lbls = paste(lbls, "%", sep = "") # add % to labels 
# Save the pie chart
png(file = "variantTypes_PieChart.png")
pie(slices, labels = lbls, col = rainbow(length(lbls)), cex = 0.8,main="Pie Chart of Variant Types")
dev.off()

# Filter variants
mut.dat = mut.dat %>%
  filter(Variant_Classification != "Silent" & Variant_Classification != "Intron" &
           Variant_Classification != "3'UTR" & Variant_Classification != "5'UTR" & 
           Variant_Classification != "3'Flank" & Variant_Classification != "5'Flank")


length(mut.dat$Hugo_Symbol) #241946 entries after filtering
length(unique(mut.dat$Hugo_Symbol)) #839 unique genes after filtering
length(unique(mut.dat$Tumor_Sample_Barcode)) #34853 unique patients after filtering

###################################
### Summary Stats of GENIE data ###
###################################
mut.uni.genes = unique(mut.dat$Hugo_Symbol) #(839) all unique genes in the mut data
CNA.uni.genes = unique(CNA.dat$Hugo_Symbol) #(740) all unique genes in the CNA data (5 less than the total due to overlap)
fus.uni.genes = unique(fus.dat$Hugo_Symbol) #(1748) all unique genes in the fusion data
all.uni.genes = unique(c(CNA.uni.genes, mut.uni.genes, fus.uni.genes)) #(2153) all unique genes (mut and CNA)
pat.uni = unique(c(names(CNA.dat), mut.dat$Tumor_Sample_Barcode, fus.dat$Tumor_Sample_Barcode)) #(36492) all unique patient IDs (this includes "Hugo_symbol though!!!)
total = length(pat.uni) - 1 #(36491) total number of patients (mut, CNA, fus)

length(unique(mut.dat$Tumor_Sample_Barcode))/total * 100 #95.51122% of patients with mutation
length(unique(names(CNA.dat)))/total * 100 #77.82741% of patients with CNA
length(unique(fus.dat$Tumor_Sample_Barcode))/total * 100 #7.615576% of patients with fusion

#########################
### Import panel data ###
#########################
G_onc_act = as.character(read.csv("GENIE Onco_actionable_genes.csv")$Gene.Names)
G_MSK_act = as.character(unique(read.csv("GENIE MSK_actionable_genes.csv")$Gene))
G_ion_act = as.character(read.csv("GENIE Ion_actionable_genes.csv")$Gene.Names)
G_tru_act = as.character(read.csv("GENIE TruSeq_actionable_genes.csv")$Gene.Names)
G_FM_act = as.character(read.csv("GENIE Foundation Medicine_actionable_genes.csv")$Gene.Names)
ActSeq_act = as.character(read.csv("ActionSeq_actionable_genes.csv")$Gene.Names)
CKB_act = as.character(read.csv("CKB_actionable_genes.csv")$Gene.Names)
IPG_act = as.character(read.csv("IPG_actionable_genes.csv")$Gene.Names)

##########################################
### Correction for alternatiove naming ###
##########################################
"
1. Figure out which genes have alternate names
  > Account for alternative notations!!!!
    - See which pairs of alternate spellings are BOTH present in the dataset
    - Decide on a single spelling for each pair and change all occurrences in the dataframe  
"
# MLL (KMT2A)
# MLL2 (KMT2D)
# MYCL1 (MYCL)
# C11 or f30 (EMSY)
# AMER1 (FAM123B)
# KMT2C (MLL3)
# MRE11 (MRE11A)
"
2. Create a function that changes the genes names of a given list of genes
  > Make the function
  > Apply the function to all gene lists to normalize
"
change.name.func = function(names) { 
  #Set equal to original vector/column to replace all common alternative names 
  #with one name given a list of gene names. 
  for(i in 1:length(names)) { 
    if(names[i] == "KMT2A") {
      names[i] = "MLL"
    }
    if(names[i] == "KMT2D") {
      names[i] = "MLL2"
      print [i]
    }
    if(names[i] == "MYCL") {
      names[i] = "MYCL1"
    }
    if(names[i] == "EMSY") {
      names[i] = "C11"
    }
    if(names[i] == "FAM123B") {
      names[i] = "AMER1"
    }
    if(names[i] == "KMT2C") {
      names[i] = "MLL3"
    }
    if(names[i] == "MRE11A") {
      names[i] = "MRE11"
    }
  }
  names = return(names)
}

G_onc_act = change.name.func(G_onc_act)
G_MSK_act = change.name.func(G_MSK_act)
G_ion_act = change.name.func(G_ion_act)
G_tru_act = change.name.func(G_tru_act)
G_FM_act = change.name.func(G_FM_act)
ActSeq_act = change.name.func(ActSeq_act)
CKB_act = change.name.func(CKB_act)
IPG_act = change.name.func(IPG_act)

mut.dat$Hugo_Symbol = change.name.func(mut.dat$Hugo_Symbol)
CNA.dat$Hugo_Symbol = change.name.func(CNA.dat$Hugo_Symbol)
fus.dat$Hugo_Symbol = change.name.func(fus.dat$Hugo_Symbol)
CKB.dat$gene_symbol = change.name.func(CKB.dat$gene_symbol)

"
3. Make a dataframe with all of the gene names from each panel
  > Column 1 is the longest dataset (MSK)
  > Add each subsequent column one by one and fill all empty spaces with 'NA'
  > Can go back and check to make sure after that all of these names were changed
"
panels = data.frame(G_MSK_act)
new.col = G_onc_act
new.col2 = G_ion_act
new.col3 = G_tru_act
new.col4 = G_FM_act
new.col5 = ActSeq_act
new.col6 = CKB_act
new.col7 = IPG_act
panels$G_onc_act = c(new.col, rep(NA, nrow(panels)-length(new.col)))
panels$G_ion_act = c(new.col2, rep(NA, nrow(panels)-length(new.col2)))
panels$G_tru_act = c(new.col3, rep(NA, nrow(panels)-length(new.col3)))
panels$G_FM_act = c(new.col4, rep(NA, nrow(panels)-length(new.col4)))
panels$ActSeq_act = c(new.col5, rep(NA, nrow(panels)-length(new.col5)))
panels$CKB_act = c(new.col6, rep(NA, nrow(panels)-length(new.col6)))
panels$IPG_act = c(new.col7, rep(NA, nrow(panels)-length(new.col7)))
write.csv(panels, file = "ActionableGenePanels.csv")
#####################################
### Number of genes in each panel ###
#####################################
"
Summary Notes:
  > G_onc_act has 300 genes and is used by DFCI (12010) = 12010 patients
  > G_MSK_act has 344 genes and is used by MSK (15036) = 15036 patients 
  > G_ion_act has 50 genes and is used by GRCC (632), JHU (2647), and MDA (2605) = 5884 patients
  > G_tru_act has 48 genes and is used by NKI (223) and UHN (1406) = 1629 patients
  > G_FM_act has 322 genes and is used by VICC (1932) = 1932 patients
  > ActSeq_act has 212 genes
  > CKB_act has 86 genes
  > IPG_act has 162 genes
"
length(G_onc_act) #300-used by DFCI
length(G_MSK_act) #344-used by MSK
length(G_ion_act) #50-used by GRCC and JHU and MDA
length(G_tru_act) #48-used by NKI and UHN
length(G_FM_act) #322-used by VICC
length(ActSeq_act) #212
length(CKB_act) #86
length(IPG_act) #162

########################################
### Number of patients in each panel ###
########################################
"
1. Account for different panel sequencing different genes (changes the denominator when calculating frequency)
  > 8 total panels (NKI  MDA  UHN  VICC MSK  DFCI GRCC JHU)
  > Find the number of unique patients in each center (for mut.dat)
  > Find the number of unique patients in each center (for CNA.dat)
  > **Fusion data does not contribute any new patient data so it is not considered
  > Add both values when necessary and create variables for total number of patients in each institution
  > Check that totals add up.... THEY DO!
  > Make a function that returns the total patients sequenced for a given gene
"
length(unique(mut.dat$Center)) #8 different panel centers 
nrow(unique(mut.dat[,c("Center","Tumor_Sample_Barcode")])) #34853 unique patients/centers

# number of unique patients from each center with mutation data
p.m.NKI = length(which(unique(mut.dat[,c("Center","Tumor_Sample_Barcode")])$Center == "NKI")) #223
p.m.MDA = length(which(unique(mut.dat[,c("Center","Tumor_Sample_Barcode")])$Center == "MDA")) #2605
p.m.UHN = length(which(unique(mut.dat[,c("Center","Tumor_Sample_Barcode")])$Center == "UHN")) #1406
p.m.VICC = length(which(unique(mut.dat[,c("Center","Tumor_Sample_Barcode")])$Center == "VICC")) #1930
p.m.MSK = length(which(unique(mut.dat[,c("Center","Tumor_Sample_Barcode")])$Center == "MSK")) #13700
p.m.DFCI = length(which(unique(mut.dat[,c("Center","Tumor_Sample_Barcode")])$Center == "DFCI")) #11710
p.m.GRCC = length(which(unique(mut.dat[,c("Center","Tumor_Sample_Barcode")])$Center == "GRCC")) #632
p.m.JHU = length(which(unique(mut.dat[,c("Center","Tumor_Sample_Barcode")])$Center == "JHU")) #2647

total_mut_pt=p.m.NKI + p.m.MDA + p.m.UHN + p.m.VICC + p.m.MSK + p.m.DFCI + p.m.GRCC + p.m.JHU #34853 total (= unique(mut.dat$Tumor_Sample_Barcode))

library(data.table)
# number of unique patients from each center with CNV data
p.c.NKI = length(which(unique(names(CNA.dat)) %like% "GENIE.NKI." == TRUE)) #0
p.c.MDA = length(which(names(CNA.dat) %like% "GENIE.MDA." == TRUE)) #0
p.c.UHN = length(which(names(CNA.dat) %like% "GENIE.UHN." == TRUE)) #0
p.c.VICC = length(which(names(CNA.dat) %like% "GENIE.VICC." == TRUE)) #1353
p.c.MSK = length(which(names(CNA.dat) %like% "GENIE.MSK." == TRUE)) #15036
p.c.DFCI = length(which(names(CNA.dat) %like% "GENIE.DFCI." == TRUE)) #12010
p.c.GRCC = length(which(names(CNA.dat) %like% "GENIE.GRCC." == TRUE)) #0
p.c.JHU = length(which(names(CNA.dat) %like% "GENIE.JHU." == TRUE)) #0

total_cnv_pt=p.c.VICC + p.c.MSK + p.c.DFCI #28399 total (= length(names(CNA.dat)) - 1)

tot.NKI = p.m.NKI #223
tot.MDA = p.m.MDA #2605
tot.UHN = p.m.UHN #1406
tot.VICC = length(unique(c(unique(mut.dat$Tumor_Sample_Barcode[which(mut.dat$Tumor_Sample_Barcode %like% "GENIE.VICC." == TRUE)]), 
             names(CNA.dat)[which(unique(names(CNA.dat)) %like% "GENIE.VICC." == TRUE)]))) #1932
tot.MSK =  length(unique(c(unique(mut.dat$Tumor_Sample_Barcode[which(mut.dat$Tumor_Sample_Barcode %like% "GENIE.MSK." == TRUE)]), 
                           names(CNA.dat)[which(unique(names(CNA.dat)) %like% "GENIE.MSK." == TRUE)]))) #15036
tot.DFCI = length(unique(c(unique(mut.dat$Tumor_Sample_Barcode[which(mut.dat$Tumor_Sample_Barcode %like% "GENIE.DFCI." == TRUE)]), 
                           names(CNA.dat)[which(unique(names(CNA.dat)) %like% "GENIE.DFCI." == TRUE)]))) #12010
tot.GRCC = p.m.GRCC #632
tot.JHU = p.m.JHU #2647

length(unique(c(mut.dat$Tumor_Sample_Barcode, names(CNA.dat), fus.dat$Tumor_Sample_Barcode))) - 1 #36491
tot.NKI + tot.MDA + tot.UHN + tot.VICC + tot.MSK + tot.DFCI + tot.GRCC + tot.JHU #36491


total.func = function(gene.name) { 
  #Returns total number of patients sequenced for a particular gene
  tot = 0
  for(i in 1:length(names(panels))) {
    if(gene.name %in% panels[,i] & i == 1) {
      tot = tot + tot.MSK
    }
    if(gene.name %in% panels[,i] & i == 2) {
      tot = tot + tot.DFCI
    }
    if(gene.name %in% panels[,i] & i == 3) {
      tot = tot + tot.GRCC + tot.JHU + tot.MDA
    }
    if(gene.name %in% panels[,i] & i == 4) {
      tot = tot + tot.NKI + tot.UHN
    }
    if(gene.name %in% panels[,i] & i == 5) {
      tot = tot + tot.VICC
    }
  }
  return(tot)
}

######################################################
### Identify an optimized list of actionable genes ###
######################################################
"
## Compile all of the gene panels into a single list of all unique genes.
  > How many genes are actionable as declared by all 8 databases? 558
"
all.genes = c(G_onc_act,G_MSK_act,G_ion_act, G_tru_act, G_FM_act, ActSeq_act, CKB_act, IPG_act)
list.unique = unique(all.genes)
length(list.unique) #552 unique actionable genes from all the panels

"
## Create a dataframe with each gene and whether or not it appears in a specific database, 
then sum them up in the last column
  > Create an empty matrix with length(list.unique) rows and length(names(panels)) + 2 columns
  > Column names correspond to different gene panels
  > Put gene names in column 1
  > Create a for loop that inserts a 1 whenever a gene is present in that database
  > Sum the numbers in each row 
  > Order the dataframe based on the 'Score' column (decreasing = TRUE)
  > Determine the number of genes 8, >7, >6, >5 databases, etc...
  > Make a comprehensive list of the final actionable gene panel (act.final)
"
panel.dat = matrix(0, nrow = length(list.unique), ncol = length(names(panels)) + 2)
colnames(panel.dat) = c("Gene", "G_MSK_act", "G_onc_act", "G_ion_act", "G_tru_act", "G_FM_act",
                        "ActSeq_act", "CKB_act", "IPG_act", "Score")
panel.dat = as.data.frame(panel.dat)
panel.dat[,1] = as.character(list.unique)

for(j in 1:length(list.unique)) { #558
  for(i in 1:length(names(panels))) {
    if(panel.dat$Gene[j] %in% as.character(unlist(panels[i]))) {
      panel.dat[j,i+1] = 1
    }
  }
  panel.dat[j,10] = sum(panel.dat[j,2:9])
}

panel.dat = panel.dat[order(panel.dat$Score, decreasing = TRUE),]

length(which(panel.dat$Score == 8)) #41
length(which(panel.dat$Score >= 7)) #49
length(which(panel.dat$Score >= 6)) #58
length(which(panel.dat$Score >= 5)) #101
length(which(panel.dat$Score >= 4)) #155
length(which(panel.dat[,2] == 1 & panel.dat[,3] == 1 & panel.dat[,4] == 1 & panel.dat[,5] == 1)) #48 Final

act.final = panel.dat$Gene[which(panel.dat$G_MSK_act == 1 & panel.dat$G_onc_act == 1 & panel.dat$G_ion_act & 
                       panel.dat$G_tru_act == 1 & panel.dat$G_FM_act == 1)]
act.final #final 48 gene panel 

########################################################################################################################################################################


#################################
### GENIE Pan-cancer analysis: co-occurring actionable genes###
#################################
"
1. Subset the actionable genes (mutations data, CNA data, and fusion data) from the complete datasets
  > Create a which statement that will accumulate all the actionable genes in the GENIE datasets
  > Calculate the percentage of genes in the dataset that were actionable
    - Mutation data
    - CNA data
    - Fusion data
"
act.mut = mut.dat[which(mut.dat$Hugo_Symbol %in% act.final == TRUE), 1:19]
head(act.mut)

length(unique(mut.dat$Hugo_Symbol)) #839 different genes mutated
length(unique(act.mut$Hugo_Symbol)) #48 different actionable genes mutated
length(unique(act.mut$Hugo_Symbol))/length(unique(mut.dat$Hugo_Symbol)) * 100 #5.721097% of the mutated genes are actionable

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

act.CNA = CNA.dat[which(CNA.dat$Hugo_Symbol %in% act.final == TRUE),]

length(unique(CNA.dat$Hugo_Symbol)) #739
length(unique(act.CNA$Hugo_Symbol)) #48
length(unique(act.CNA$Hugo_Symbol))/length(unique(CNA.dat$Hugo_Symbol)) * 100 #6.495264% of the CNA genes are actionable

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

act.fus = fus.dat[which(fus.dat$Hugo_Symbol %in% act.final == TRUE),]

length(unique(fus.dat$Hugo_Symbol)) #1747
length(unique(act.fus$Hugo_Symbol)) #43
length(unique(act.fus$Hugo_Symbol))/length(unique(fus.dat$Hugo_Symbol)) * 100 #2.461362% of the fusion genes are actionable

########################################################################################################################################################################

"
2. Making variables... 
"
mut.uni.genes = unique(mut.dat$Hugo_Symbol) #(839) all unique genes in the mut data
CNA.uni.genes = unique(CNA.dat$Hugo_Symbol) #(740) all unique genes in the CNA data (5 less than the total due to overlap)
fus.uni.genes = unique(fus.dat$Hugo_Symbol) #(1748) all unique genes in the fusion data
all.uni.genes = unique(c(CNA.uni.genes, mut.uni.genes, fus.uni.genes)) #(2153) all unique genes (mut and CNA)
pat.uni = unique(c(names(CNA.dat), mut.dat$Tumor_Sample_Barcode, fus.dat$Tumor_Sample_Barcode)) #(36492) all unique patient IDs (this includes "Hugo_symbol though!!!)
total = length(pat.uni) - 1 #(36491) total number of patients (mut, CNA, fus)

all.act.uni = unique(c(act.mut$Hugo_Symbol, act.CNA$Hugo_Symbol, act.fus$Hugo_Symbol)) #(48) all unique actionable genes (same as act.final)
pts_mut.act.uni = unique(mut.dat$Tumor_Sample_Barcode[which(mut.dat$Hugo_Symbol %in% act.mut$Hugo_Symbol)]) #(30295) unique patient IDs for actionable mutations
pts_CNA.act.uni = unique(names(act.CNA)) #(28400) unique patient IDs for actionable CNA
pts_fus.act.uni = unique(fus.dat$Tumor_Sample_Barcode[which(fus.dat$Hugo_Symbol %in% act.fus$Hugo_Symbol)]) #(840) unique patient IDs for actionable fusions

pat.act.uni = unique(c(pts_CNA.act.uni,pts_mut.act.uni,pts_fus.act.uni)) #(36154) unique patient IDs for actionable mutations, CNA or fusion data (union of patients with one of the three datasets)

########################################################################################################################################################################

"
3. Finding which genes are co-occurring
  > Make an empty matrix/dataframe with actionable gene rows and patient columns
  > Make for loops that indicate which what mutation that gene has
    -KEY, recorded number for each of the pattern below, and the count of patients with the pattern
      -- Mutated = 1 (69155 patients with mutations in actionable genes)
      -- CNA = 2 (8477 patients with CNA in actionable genes)
      -- Fusion = 3 (812 patients with fusions in actionable genes)
      -- Mutated & CNA = 4 (737 patients with only mutations and CNA in the same actionable genes) 
      -- Mutated & Fusion = 5 (NONE. 0 patients with only mutations and fusion in the same actionable genes)
      -- CNA & Fusion = 6 (NONE. 0 patients with only CNA and fusion in the same actionable genes)
      -- Mutated & CNA & Fusion = 7 (48 patients with all three, mutations, CNA and fusion in the same actionable genes)
  > Make a for loop that identifies which columns have <= 1 mutation (patients do not have an co-occurring actionable gene mutation)
  > Remove all patient columns that don't have co-occurring mutations
  > KEEP all patients with co-occurring actionable genes in 'by.pat.mat' BEFORE REMOVAL STORED IN 'ALL.by.pat.mat', with both co-occurring and not-co-occurring actionable genes
"

num_actionableGenes=length(all.act.uni)
num_pts=length(pat.act.uni) # total number of unique patients with at lease one actionable gene out of the three, mutations, CNA or funsion data. 
by.pat.mat = matrix(0, nrow = num_actionableGenes, ncol = num_pts) #empty binomial matrix
colnames(by.pat.mat) = c(pat.act.uni)
by.pat.mat = as.data.frame(by.pat.mat) #dataframe containing actionable genes (rows) and patient IDs (cols)
by.pat.mat$Hugo_Symbol = all.act.uni

# For the actionable gene, patients' matrix, fill in mutation status for each gene and each patient.
# Used only mutations with at least one of actionable genes. For single gene ratio, it will be elevated due to selected patients.
# This is ok for ranking co-occurring actionable genes based on each actionable genes' statistics.
count_mut = 1 # (4 minutes)
total_mut_action=length(act.mut$Hugo_Symbol)
print ('Total number of patients with mutations covering actionable genes: ')
print (total_mut_action)
for(i in 1:total_mut_action) { #82168
  pat = act.mut$Tumor_Sample_Barcode[i]
  gene = act.mut$Hugo_Symbol[i]
  by.pat.mat[which(by.pat.mat$Hugo_Symbol == gene), which(names(by.pat.mat) == pat)] = 1
  if (count_mut==10000) {print(count_mut)}
  if (count_mut==50000) {print(count_mut)}
  count_mut = count_mut + 1
}

# For the actionable gene, patients' matrix, fill in CNA status for each gene and each patient.
# both mutation & CNA for the same gene and same patient, record as 4.
# CNA only for the patient, record as 2.
# CNA values: -2,-1,0,1,2; only consider -2 and 2 for this version. strict criteria for CNA.
# act.CNA[i,] == 1 | act.CNA[i,] == -1
count_CNA = 1
total_CNA_action=length(act.CNA$Hugo_Symbol)
print ('Total number of patients with CNA covering actionable genes: ')
print (total_CNA_action)
for(i in 1:total_CNA_action) { #48
  temp = which(act.CNA[i,] == 2 | act.CNA[i,] == -2)
  for(j in 1:length(temp)) {
    if(by.pat.mat[which(by.pat.mat$Hugo_Symbol == act.CNA$Hugo_Symbol[i]), temp[j]] == 1) {
      by.pat.mat[which(by.pat.mat$Hugo_Symbol == act.CNA$Hugo_Symbol[i]), temp[j]] = 4 # both mutation & CNA for the same gene and same patient, record as 4.
    } else{
      by.pat.mat[which(by.pat.mat$Hugo_Symbol == act.CNA$Hugo_Symbol[i]), temp[j]] = 2 # CNA only for the patient, record as 2.
    }
  }
  #print(count_CNA)
  count_CNA = count_CNA + 1
}

# For the actionable gene, patients' matrix, fill in fusion status for each gene and each patient.
# both mutation & fusion for the same gene and same patient, record as 5.
# both CNA & fusion for the same gene and same patient, record as 6.
# all three, mutation, CNA & fusion for the same gene and same patient, record as 7.
# fusion only for the patient, record as 3.
count_fusion = 1
total_fus_action=length(act.fus$Hugo_Symbol)
print ('Total number of patients with fusions covering actionable genes: ')
print (total_fus_action)
for(i in 1:total_fus_action) { #939
  pat = act.fus$Tumor_Sample_Barcode[i]
  gene = act.fus$Hugo_Symbol[i]
  if(by.pat.mat[which(by.pat.mat$Hugo_Symbol == gene), which(names(by.pat.mat) == pat)] == 1) {
    by.pat.mat[which(by.pat.mat$Hugo_Symbol == gene), which(names(by.pat.mat) == pat)] = 5
  }
  if(by.pat.mat[which(by.pat.mat$Hugo_Symbol == gene), which(names(by.pat.mat) == pat)] == 2) {
    by.pat.mat[which(by.pat.mat$Hugo_Symbol == gene), which(names(by.pat.mat) == pat)] = 6
  }
  if(by.pat.mat[which(by.pat.mat$Hugo_Symbol == gene), which(names(by.pat.mat) == pat)] == 4) {
    by.pat.mat[which(by.pat.mat$Hugo_Symbol == gene), which(names(by.pat.mat) == pat)] = 7
  } else {
    by.pat.mat[which(by.pat.mat$Hugo_Symbol == gene), which(names(by.pat.mat) == pat)] = 3
  }
  if (count_fusion==100) {print(count_fusion)}
  if (count_fusion==500) {print(count_fusion)}
  count_fusion = count_fusion + 1
}

# Summarize the number of patients with variant patterns:
'
-- Mutated = 1 (patients with mutations in actionable genes)
-- CNA = 2 (patients with CNA in actionable genes)
-- Fusion = 3 (patients with fusions in actionable genes)
-- Mutated & CNA = 4 (patients with only mutations and CNA in the same actionable genes) 
-- Mutated & Fusion = 5 (patients with only mutations and fusion in the same actionable genes)
-- CNA & Fusion = 6 (patients with only CNA and fusion in the same actionable genes)
-- Mutated & CNA & Fusion = 7 (patients with all three, mutations, CNA and fusion in the same actionable genes)
'
nrow(which(by.pat.mat == 1, arr.ind = TRUE)) #69155
nrow(which(by.pat.mat == 2, arr.ind = TRUE)) #8477
nrow(which(by.pat.mat == 3, arr.ind = TRUE)) #812
nrow(which(by.pat.mat == 4, arr.ind = TRUE)) #737
nrow(which(by.pat.mat == 5, arr.ind = TRUE)) #0
nrow(which(by.pat.mat == 6, arr.ind = TRUE)) #0
nrow(which(by.pat.mat == 7, arr.ind = TRUE)) #48

ALL.by.pat.mat = by.pat.mat

num_act_patients=length(names(by.pat.mat))
print ('Total number of patients with at least one actionable genes: ')
print (num_act_patients) ##36154
bad = c() # Not co-occurring patients
count_NotCoOccur = 1
for(i in 2:length(names(by.pat.mat))) { 
  if(length(which(by.pat.mat[,i] > 0)) < 2) {
    bad[i-1] = i
    count_NotCoOccur = count_NotCoOccur + 1
  }
}
print ('Number of patients without co-occurring actionable genes: ')
print(count_NotCoOccur) #15261

bad = bad[!is.na(bad)]
by.pat.mat[,bad] = NULL #now by.pat.mat has all the co-occurring actionable genes in it (20894)

by.pat.mat[,1:5] #has ONLY patients with co-occurring actionable genes
ALL.by.pat.mat[,1:5] #has patients with BOTH co-occurring and non-co-occurring actionable genes

########################################################################################################################################################################
"
## Create generic functions
  > Make a function that receives 2 gene names as input and returns the number of times they occur together
  > freq.gene.func find the frequency at which a given gene is mutated
  > exp.p is the expected probability of the two genes co-occurring
    - pA = instances of mutated A/total patients sequenced for that gene (using total.func())
    - pB = instances of mutated B/total patients sequenced for that gene (using total.func())
  > my.binom.func applies the binom.test of two genes
"
co.counts.func = function(gene1, gene2) { 
  #Returns the number of times 2 given genes occur together in the same patient
  co.num = 0
  co.num = length(which(which(by.pat.mat[which(by.pat.mat$Hugo_Symbol == gene1),] > 0) %in% which(
    by.pat.mat[which(by.pat.mat$Hugo_Symbol == gene2),] > 0) == TRUE))
  return(co.num)
}

freq.gene.func = function(g) { 
  #Returns the frequency (as a percent) that a given actionable gene is altered
  ((length(which(by.pat.mat[which(by.pat.mat$Hugo_Symbol == g),] == 1)) + 
      length(which(by.pat.mat[which(by.pat.mat$Hugo_Symbol == g),] == 2)) + 
      length(which(by.pat.mat[which(by.pat.mat$Hugo_Symbol == g),] == 3)) +
      2*length(which(by.pat.mat[which(by.pat.mat$Hugo_Symbol == g),] == 4)) +
      3*length(which(by.pat.mat[which(by.pat.mat$Hugo_Symbol == g),] == 7)))/
     total.func(g)) *100
}

exp.p.func = function(gene1, gene2) { 
  #Returns the expected probability of two genes to occur together randomly
  l.gene1 = length(which(by.pat.mat[which(ALL.by.pat.mat$Hugo_Symbol == gene1),] == 1)) + 
    length(which(ALL.by.pat.mat[which(ALL.by.pat.mat$Hugo_Symbol == gene1),] == 2)) + 
    length(which(ALL.by.pat.mat[which(ALL.by.pat.mat$Hugo_Symbol == gene1),] == 3)) +
    2*length(which(ALL.by.pat.mat[which(ALL.by.pat.mat$Hugo_Symbol == gene1),] == 4)) +
    2*length(which(ALL.by.pat.mat[which(ALL.by.pat.mat$Hugo_Symbol == gene1),] == 5)) +
    2*length(which(ALL.by.pat.mat[which(ALL.by.pat.mat$Hugo_Symbol == gene1),] == 6)) +
    3*length(which(ALL.by.pat.mat[which(ALL.by.pat.mat$Hugo_Symbol == gene1),] == 7))
  l.gene2 = length(which(ALL.by.pat.mat[which(ALL.by.pat.mat$Hugo_Symbol == gene2),] == 1)) + 
    length(which(ALL.by.pat.mat[which(ALL.by.pat.mat$Hugo_Symbol == gene2),] == 2)) + 
    length(which(ALL.by.pat.mat[which(ALL.by.pat.mat$Hugo_Symbol == gene2),] == 3)) +
    2*length(which(ALL.by.pat.mat[which(ALL.by.pat.mat$Hugo_Symbol == gene2),] == 4)) +
    2*length(which(ALL.by.pat.mat[which(ALL.by.pat.mat$Hugo_Symbol == gene1),] == 5)) +
    2*length(which(ALL.by.pat.mat[which(ALL.by.pat.mat$Hugo_Symbol == gene1),] == 6)) +
    3*length(which(ALL.by.pat.mat[which(ALL.by.pat.mat$Hugo_Symbol == gene2),] == 7))
  a1 = l.gene1/total.func(gene1) 
  b1 = l.gene2/total.func(gene2) 
  p1 = a1*b1 
  return(p1)
}

total.func2 = function(gene.name1, gene.name2) { 
  #Returns total number of patients sequenced for 2 genes
  tot = 0
  for(i in 1:length(names(panels))) {
    if(gene.name1 %in% panels[,i] & gene.name2 %in% panels[,i] & i == 1) {
      tot = tot + tot.MSK
    }
    if(gene.name1 %in% panels[,i] & gene.name2 %in% panels[,i] & i == 2) {
      tot = tot + tot.DFCI
    }
    if(gene.name1 %in% panels[,i] & gene.name2 %in% panels[,i] & i == 3) {
      tot = tot + tot.GRCC + tot.JHU + tot.MDA
    }
    if(gene.name1 %in% panels[,i] & gene.name2 %in% panels[,i] & i == 4) {
      tot = tot + tot.NKI + tot.UHN
    }
    if(gene.name1 %in% panels[,i] & gene.name2 %in% panels[,i] & i == 5) {
      tot = tot + tot.VICC
    }
  }
  return(tot)
}

my.binom.func = function(gene1, gene2) { #Returns the binomial test results for the two genes (whether or not they occur more)
  sums = length(which(which(by.pat.mat[which(by.pat.mat$Hugo_Symbol == gene1),] > 0) %in% 
                        which(by.pat.mat[which(by.pat.mat$Hugo_Symbol == gene2),] > 0) == TRUE))
  tots = total.func2(gene1, gene2) 
  return(binom.test(sums, n = tots, p = exp.p.func(gene1, gene2), alternative = "greater"))
}

########################################################################################################################################################################
# "
# 4. Make a binary matrix to input counts of gene co-occurrence
#   > Create an empty matrix/dataframe with rows and columns equaling actionable genes
#   > Make a for loop that considers every gene pair in the square matrix and puts its corresponding co-occurrence
#   > Create a table with all possible gene pairs
# "
my.gene.mat = matrix(0, ncol = (length(by.pat.mat$Hugo_Symbol) + 1), nrow = (length(by.pat.mat$Hugo_Symbol)))
my.gene.mat[,1] = c(by.pat.mat$Hugo_Symbol)
my.gene.mat = data.frame(my.gene.mat, stringsAsFactors = FALSE)
colnames(my.gene.mat) = c("Hugo_Symbol", by.pat.mat$Hugo_Symbol)

########################################################################################################################################################################
"
5. Make a table to hold all of the data
  > First 2 columns should be all possible pairs of genes
  > Remove repeats
"
rep1 = rep(my.gene.mat$Hugo_Symbol, times = length(my.gene.mat))
rep2 = rep(my.gene.mat$Hugo_Symbol, each = length(my.gene.mat))
wait = data.frame(rep1, rep2, stringsAsFactors = FALSE)
for (i in 1:nrow(wait)){
  wait[i,] = sort(wait[i,])
}
wait = wait[!duplicated(wait),]
rep1 = wait[,1]
rep2 = wait[,2]

table.dat = data.frame(rep1, rep2, mutrep1 = numeric(length(rep1)), CNArep1 = numeric(length(rep1)),
                       totrep1 = numeric(length(rep1)), mutrep2 = numeric(length(rep1)),
                       CNArep2 = numeric(length(rep1)), totrep2 = numeric(length(rep1)),
                       count.obs = numeric(length(rep1)), count.exp = numeric(length(rep1)),
                       pos = numeric(length(rep1)), p.val = numeric(length(rep1)), 
                       stringsAsFactors = FALSE)

########################################################################################################################################################################

"
6. Fill in the co-occurrent table with count of variants not-co-occurring and co-occurring, and p_values
"
count = 1 #VERY LONG RUN OVERNIGHT
length(rep1) #1176
for(i in 1:length(rep1)) { 
  gene1 = table.dat$rep1[i]
  gene2 = table.dat$rep2[i]
  table.dat$mutrep1[i] = length(which(ALL.by.pat.mat[which(ALL.by.pat.mat$Hugo_Symbol == gene1),] %in% c(1,4,5,7) == TRUE))
  table.dat$CNArep1[i] = length(which(ALL.by.pat.mat[which(ALL.by.pat.mat$Hugo_Symbol == gene1),] %in% c(2,4,6,7) == TRUE))
  table.dat$totrep1[i] = table.dat$mutrep1[i] + table.dat$CNArep1[i]
  table.dat$mutrep2[i] = length(which(ALL.by.pat.mat[which(ALL.by.pat.mat$Hugo_Symbol == gene2),] %in% c(1,4,5,7) == TRUE))
  table.dat$CNArep2[i] = length(which(ALL.by.pat.mat[which(ALL.by.pat.mat$Hugo_Symbol == gene2),] %in% c(2,4,6,7) == TRUE))
  table.dat$totrep2[i] = table.dat$mutrep2[i] + table.dat$CNArep2[i]
  table.dat$count.obs[i] = co.counts.func(gene1, gene2)
  table.dat$count.exp[i] = exp.p.func(gene1, gene2)*total
  table.dat$pos[i] = my.binom.func(gene1, gene2)$estimate
  table.dat$p.val[i] = my.binom.func(gene1, gene2)$p.value
  if (count==100) {print(count)}
  if (count==500) {print(count)}
  if (count==1000) {print(count)}
  count = count + 1
}

table.dat[order(table.dat$p.val),]

which(table.dat$rep1 == table.dat$rep2)
table.dat2 = table.dat[-which(table.dat$rep1 == table.dat$rep2),]
table.dat2[order(table.dat2$p.val),]
length(names(by.pat.mat))/total * 100 #57.25795% of patients had co-occurring alterations

########################################################################################################################################################################

"
7. Benjamini & Hochberg vs. Bonferroni correction of p-values
  > Helps avoid Type I errors (false positives)
  > Out of 100 tests, 5 are expected to randomly appear significant
  > Benjamini & Hochberg
    - p-values sorted then ranked
    - each p-value is multiplied by N (total number of comparisons) and divided by its rank to give the adjusted p-values
    - controls the false discovery rate (FDR)
  > Bonferroni
    - more conservative than FDR
    - controls experiment-wide false positive value
    - specifies a new alpha value for each individual test
"
table.dat2[,ncol(table.dat2) + 1] = numeric(nrow(table.dat2))
colnames(table.dat2)[ncol(table.dat2)] = "FDR p.val"
table.dat2$`FDR p.val` = p.adjust(table.dat2$p.val, method = "fdr")
table.dat2=table.dat2[order(table.dat2$`FDR p.val`),]

# Save the orderred co-occurring actionable genes
write.csv(table.dat2, file="Ordered_CoOccurringGenes.csv", sep = "\t")

# Top Five lowest p-values from table.dat
#         rep1    rep2 mutrep1 CNArep1 totrep1 mutrep2 CNArep2 totrep2 count.obs   count.exp      pos         p.val     FDR p.val
# 318      APC    KRAS    4033      48    4081    6187     430    6617      1547  704.689238 0.0423940150 2.839002e-169 3.202394e-166
# 324     KRAS   SMAD4    6187     430    6617    1432     196    1628       777  254.325286 0.0212929215 6.447123e-154 3.636177e-151
# 589      KIT  PDGFRA     966     191    1157     867     198    1065       245   27.269738 0.0067139843 1.350946e-141 5.079555e-139
# 282      KDR  PDGFRA    1015     151    1166     867     198    1065       253   32.547752 0.0069332164 3.884365e-132 1.095391e-129
# 289      KDR     KIT    1015     151    1166     966     191    1157       253   35.285413 0.0069332164 1.934869e-124 4.365065e-122

length(which(table.dat2$`FDR p.val` < 0.05)) #973
nrow(table.dat2) #out of 1128
length(which(table.dat2$`FDR p.val` < 0.05))/nrow(table.dat2) * 100 #86.26%

# Save the histogram of p values of co-occurring genes
png(file = "P_value_distribution_CooccurringActionableGenes.png")
histogram_p_dist=hist(table.dat2$`FDR p.val`, las = 1, xlab = "p-value", cex.lab = 1.5, main = "Distribution of p-values for co-occurring actionable genes")
dev.off()

########################################################################################################################################################################
########################################################################################################################################################################
########################################################################################################################################################################
########################################################################################################################################################################