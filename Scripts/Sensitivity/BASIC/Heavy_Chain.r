rm(list = ls())
library(dplyr)
library(tidyr)
library(stringr)

#Load Ground Truth data Heavy Chain
HC <- read.table("HC.txt",header=T)
dim(HC)
HC

#Exclude single end and low coverage samples which were not run
single_end <- c('PW1-A1_HC', 'PW1-A4_HC',"PW1-A6_HC","PW1-B1_HC","PW1-B3_HC","PW1-B4_HC","PW1-B6_HC","PW1-C1_HC","PW1-C4_HC","PW1-C5_HC","PW2_A5_HC","PW2_A6_HC","PW2_A9_HC","PW2_B12_HC","PW2_C5_HC","PW2_D11_HC","PW2_D12_HC","PW2_D7_HC","PW2_E12_HC","PW2_G10_HC","PW2_G12_HC","PW2_G4_HC","PW2_G5_HC","PW2_G6_HC","PW2_G7_HC","PW2_H2_HC","PW2_H3_HC","PW2_H4_HC","PW2_H5_HC","PW3_B11_HC","PW3_B9_HC","PW3_C1_HC","PW3_C10_HC","PW3_C11_HC","PW3_C12_HC","PW3_C2_HC","PW3_C3_HC","PW3_C4_HC","PW3_C5_HC","PW3_C7_HC","PW3_C9_HC","PW3_D1_HC","PW3_D2_HC","PW3_D4_HC","PW3_D5_HC","PW3_D6_HC","PW3_F1_HC","PW3_F10_HC","PW3_F12_HC","PW3_F8_HC","PW3_F9_HC","PW3_G12_HC","PW3_G6_HC","PW3_G8_HC","PW3_G9_HC","PW3_H1_HC","PW3_H10_HC","PW3_H2_HC","PW3_H4_HC","PW3_H5_HC","PW3_H6_HC","PW3_H7_HC","PW3_H8_HC","PW3_H9_HC")
HC_no_single_end <- HC[ !grepl(paste(single_end, collapse="|"), HC$sequence_id),]
dim(HC_no_single_end)

test_HC_no_single_end <- HC_no_single_end %>% separate(sequence_id, c("A","B","C"))
test_HC_no_single_end <- test_HC_no_single_end %>% unite("Id", A:B, remove = TRUE)
Ground_Truth <- subset(test_HC_no_single_end, select = c(1,3,4,5,6))

#count
tot_chains_productive <- Ground_Truth[ which( Ground_Truth$productive == "TRUE" ),]
length(tot_chains_productive$Id)
percentage_productive <- length(tot_chains_productive$Id)/length(Ground_Truth$Id)
print(percentage_productive)

#Load Basic Data output from a given tool
#coverage 50bp and 1mln.250k

setwd("/path/to/annoated_fasta_sequence_of_HC_obtained_with_immcantation_from_a_tool/")
dataset_LC_L002 <- read.table("Results_IGH_rank_all_no_empy_lines_L002.tsv",sep="\t",header=F)
dataset_LC_L002_subset <- subset(dataset_LC_L002, dataset_LC_L002$V2 == 1)

###count number of Light Chain L and how many are productive
dataset_LC_L002_subset_productive <- subset(dataset_LC_L002_subset, dataset_LC_L002_subset$V3 == "Yes")
print(paste0("Number Heavy Chain reconstructed   ",length(dataset_LC_L002_subset$V1)))
print(paste0("Number Heavy Chain reconstructed productive   ",length(dataset_LC_L002_subset_productive$V1)))


##Parse the name and subset the chains that are both in ground truth and assembled
test <- dataset_LC_L002_subset_productive %>% separate(V1, c("A","B","C","D","E"))
head(test)
test_LC_L002_subset <- test %>% unite("Id", A:B, remove = TRUE)
dim(test)
A <- intersect(as.character(Ground_Truth$Id),as.character(test_LC_L002_subset$Id))
length(A)

head(test_LC_L002_subset)
LC_L002_BASIC_1mln_250L <- subset(test_LC_L002_subset, select = c(1,6,7,8,9))
dim(LC_L002_BASIC_1mln_250L)
names(LC_L002_BASIC_1mln_250L)[names(LC_L002_BASIC_1mln_250L) == "V3"] <- "productive"
names(LC_L002_BASIC_1mln_250L)[names(LC_L002_BASIC_1mln_250L) == "V4"] <- "v_call"
names(LC_L002_BASIC_1mln_250L)[names(LC_L002_BASIC_1mln_250L) == "V5"] <- "d_call"
names(LC_L002_BASIC_1mln_250L)[names(LC_L002_BASIC_1mln_250L) == "V6"] <- "j_call"
dim(LC_L002_BASIC_1mln_250L)

#LC_L002_BASIC_1mln_250L$productive <- charr("Yes ", "TRUE", LC_L002_BASIC_1mln_250L$productive)
LC_L002_BASIC_1mln_250L$productive <- gsub("Yes","TRUE",LC_L002_BASIC_1mln_250L$productive)
LC_L002_BASIC_1mln_250L_productive_in_ground_truth <- LC_L002_BASIC_1mln_250L[ grepl(paste(tot_chains_productive$Id, collapse="|"), LC_L002_BASIC_1mln_250L$Id),]
head(LC_L002_BASIC_1mln_250L_productive_in_ground_truth)
dim(LC_L002_BASIC_1mln_250L_productive_in_ground_truth)

dim(Ground_Truth)
#Join Table basic 1mlnk 100bp
joined_tables_all_basic <- full_join(LC_L002_BASIC_1mln_250L_productive_in_ground_truth,Ground_Truth,by = c("Id"),copy = TRUE)
joined_tables_all_basic_in_ground_truth <- subset(joined_tables_all_basic, joined_tables_all_basic$productive.y == "TRUE" | joined_tables_all_basic$productive.y == "FALSE")
dim(joined_tables_all_basic_in_ground_truth)

#Perform the match
#Match V genes
v_genes <- joined_tables_all_basic_in_ground_truth %>% mutate(id = joined_tables_all_basic_in_ground_truth$Id) %>% 
  separate_rows(v_call.x, sep = ',') %>%
  separate_rows(v_call.y, sep = ',') %>%
  mutate(match = v_call.x == v_call.y) %>% 
  mutate(match = v_call.y == v_call.x) %>% 
  group_by(id) %>% mutate(v_call.x = toString(v_call.x)) %>% 
  group_by(id) %>% mutate(v_call.y = toString(v_call.y)) %>% 
  
  mutate(match = if_else(any(match == TRUE), TRUE, FALSE)) %>% 
  distinct() %>% ungroup() %>% select(-id)


v_genes$productive.x <- as.character(v_genes$productive.x)
v_genes$productive.y <- as.character(v_genes$productive.y)
v_genes$match <- as.character(v_genes$match)
v_genes$v_call.y <- as.character(v_genes$v_call.y)

v_genes$productive.y <- replace(v_genes$productive.y, is.na(v_genes$productive.y), FALSE)
v_genes$productive.x <- replace(v_genes$productive.x, is.na(v_genes$productive.x), FALSE)
v_genes$match <- replace(v_genes$match,is.na(v_genes$match), FALSE)
v_genes$productive.y <- replace(v_genes$productive.y, is.na(v_genes$productive.y), FALSE)
v_genes$j_call.y <- replace(v_genes$j_call.y,is.na(v_genes$j_call.y), FALSE)
v_genes$d_call.y <- replace(v_genes$d_call.y, is.na(v_genes$d_call.y), FALSE)
v_genes$v_call.y <- replace(v_genes$v_call.y, is.na(v_genes$v_call.y), FALSE)

names(v_genes)[10] <- "v.gene.match"
v_genes



####################################################################


#Perform the match considering multiple genes in the output of the tool or the sanger for every V-D-J
#Match D genes
d_genes <- joined_tables_all_basic_in_ground_truth %>% mutate(id = joined_tables_all_basic_in_ground_truth$Id) %>% 
  separate_rows(d_call.x, sep = ',') %>%
  separate_rows(d_call.y, sep = ',') %>%
  mutate(match = d_call.x == d_call.y) %>% 
  mutate(match = d_call.y == d_call.x) %>% 
  group_by(id) %>% mutate(d_call.x = toString(d_call.x)) %>% 
  group_by(id) %>% mutate(d_call.y = toString(d_call.y)) %>% 
  
  mutate(match = if_else(any(match == TRUE), TRUE, FALSE)) %>% 
  distinct() %>% ungroup() %>% select(-id)


d_genes$productive.x <- as.character(d_genes$productive.x)
d_genes$productive.y <- as.character(d_genes$productive.y)
d_genes$match <- as.character(d_genes$match)
d_genes$d_call.y <- as.character(d_genes$d_call.y)

d_genes$productive.y <- replace(d_genes$productive.y, is.na(d_genes$productive.y), FALSE)
d_genes$productive.x <- replace(d_genes$productive.x, is.na(d_genes$productive.x), FALSE)
d_genes$match <- replace(d_genes$match,is.na(d_genes$match), FALSE)
d_genes$productive.y <- replace(d_genes$productive.y, is.na(d_genes$productive.y), FALSE)
d_genes$j_call.y <- replace(d_genes$j_call.y,is.na(d_genes$j_call.y), FALSE)
d_genes$d_call.y <- replace(d_genes$d_call.y, is.na(d_genes$d_call.y), FALSE)
d_genes$v_call.y <- replace(d_genes$v_call.y, is.na(d_genes$v_call.y), FALSE)

names(d_genes)[10] <- "d.gene.match"
d_genes
###################################################################################
#Match J genes

j_genes <- joined_tables_all_basic_in_ground_truth %>% mutate(id = joined_tables_all_basic_in_ground_truth$Id) %>% 
  separate_rows(j_call.x, sep = ',') %>%
  separate_rows(j_call.y, sep = ',') %>%
  mutate(match = j_call.x == j_call.y) %>% 
  mutate(match = j_call.y == j_call.x) %>% 
  group_by(id) %>% mutate(j_call.x = toString(j_call.x)) %>% 
  group_by(id) %>% mutate(j_call.y = toString(j_call.y)) %>% 
  
  mutate(match = if_else(any(match == TRUE), TRUE, FALSE)) %>% 
  distinct() %>% ungroup() %>% select(-id)


j_genes$productive.x <- as.character(j_genes$productive.x)
j_genes$productive.y <- as.character(j_genes$productive.y)
j_genes$match <- as.character(j_genes$match)
j_genes$j_call.y <- as.character(j_genes$j_call.y)

j_genes$productive.y <- replace(j_genes$productive.y, is.na(j_genes$productive.y), FALSE)
j_genes$productive.x <- replace(j_genes$productive.x, is.na(j_genes$productive.x), FALSE)
j_genes$match <- replace(j_genes$match,is.na(j_genes$match), FALSE)
j_genes$productive.y <- replace(j_genes$productive.y, is.na(j_genes$productive.y), FALSE)
j_genes$j_call.y <- replace(j_genes$j_call.y,is.na(j_genes$j_call.y), FALSE)
j_genes$d_call.y <- replace(j_genes$d_call.y, is.na(j_genes$d_call.y), FALSE)
j_genes$v_call.y <- replace(j_genes$v_call.y, is.na(j_genes$v_call.y), FALSE)


names(j_genes)[10] <- "j.gene.match"
j_genes

###################################################################################
#number of heavy chain in Ground Truth and number of chain in ground truth productive
tot_chain <- 113
tot_chains_productive <- Ground_Truth[ which( Ground_Truth$productive == "TRUE" ),]
dim(tot_chains_productive)
head(tot_chains_productive)
length(tot_chains_productive$Id)
percentage_productive <- length(tot_chains_productive$Id)/tot_chain
print(percentage_productive)


#Match productivity
df <- as.data.frame(cbind(joined_tables_all_basic_in_ground_truth$Id,v_genes$v.gene.match,d_genes$d.gene.match,v_genes$v.gene.match))

#select only the chains that were productive in the ground truth
df_productive <- df[df$V1 %in% tot_chains_productive$Id,]

Id <- df_productive$V1
df_productive$V1 <- as.character(df_productive$V1) 
df_productive$V2 <- as.logical(df_productive$V2)
df_productive$V3 <- as.logical(df_productive$V3)
df_productive$V4 <- as.logical(df_productive$V4)
df_productive <- 1*df_productive[,2:4]
head(df_productive)
df_productive$sum <- rowSums(df_productive[,1:3])/3
Accuracy <- cbind(Id,df_productive)
Alll_Match <- subset(Accuracy, sum == 1)
value_all_match <- length(Alll_Match$Id)
Accuaracy_Value <- value_all_match/length(Accuracy$Id)
print(Accuaracy_Value)
print(value_all_match)
print(length(Accuracy$Id))


print(paste0("Percentage Productive in The Ground Truth ", length(tot_chains_productive$Id)/length(Ground_Truth$Id)))
print(paste0("Number Heavy Chain reconstructed   ",length(dataset_LC_L002_subset$V1)))
print(paste0("Number Heavy Chain reconstructed productive   ",length(dataset_LC_L002_subset_productive$V1)))
print(paste0("Accuracy  ", value_all_match/length(Accuracy$Id)))


