rm(list = ls())
library(dplyr)
library(tidyr)
library(stringr)
setwd("/root/cloud-data/snf-mgln-dds/AIDA/Bioinformatics/i0439277/RTcure/Benchmark/Baldr_data/ground_truth/fasta/")

#Load Ground Truth data Heavy Chain
HC <- read.table("heavy.Id.sra.sequence_igblast_db-pass.tsv",header=T,sep="\t")
dim(HC)
head(HC)

#Exclude single end and low coverage samples which were not run

Ground_Truth <- subset(HC, select = c(1,4,5,6,7))
head(Ground_Truth)
dim(Ground_Truth)
tot_chains_productive <- Ground_Truth[ which(Ground_Truth$productive == "TRUE" ),]
length(tot_chains_productive$sequence_id)
percentage_productive <- length(tot_chains_productive$sequence_id)/length(Ground_Truth$sequence_id)
print(percentage_productive)
Ground_Truth_split <- Ground_Truth %>% separate(sequence_id, c("A","B","C"))
Ground_Truth_split_subset <- subset(Ground_Truth_split, select = c(1,4,5,6,7))
head(Ground_Truth_split_subset)
names(Ground_Truth_split_subset)[names(Ground_Truth_split_subset) == "A"] <- "sequence_id"



#Load Baldr Data
#coverage 100 and 1mln.250k

setwd("/root/cloud-data/snf-mgln-dds/AIDA/Bioinformatics/i0439277/RTcure/Benchmark/Baldr_data/bracer/transfer/paired_end/75bp/run_50cov/test/tsv_files/IGH/")
dataset_LC_L002 <- read.table("IGH_header.tsv",sep="\t",header=T,skipNul = T)
head(dataset_LC_L002)
###count number of Light Chain L and how many are productive
dataset_LC_L002_productive <- subset(dataset_LC_L002, productive == "TRUE")
print(paste0("Number Heavy Chain reconstructed   ",length(dataset_LC_L002$sequence_id)))
print(paste0("Number Heavy Chain reconstructed productive   ",length(dataset_LC_L002_productive$sequence_id)))

##Parse the name and subset the chains that are both in ground truth and assembled
LC_L002_BASIC_1mln_250L_productive_in_ground_truth <- dataset_LC_L002[ grepl(paste(Ground_Truth_split_subset$sequence_id, collapse="|"), dataset_LC_L002$sequence_id),]
dim(dataset_LC_L002)
dim(Ground_Truth_split_subset)
dim(LC_L002_BASIC_1mln_250L_productive_in_ground_truth)
head(LC_L002_BASIC_1mln_250L_productive_in_ground_truth)
head(Ground_Truth_split_subset)


#Join Table basic 1mlnk 100bp
joined_tables_all_basic <- full_join(LC_L002_BASIC_1mln_250L_productive_in_ground_truth,Ground_Truth_split_subset,by = c("sequence_id"),copy = TRUE)

joined_tables_all_basic_in_ground_truth <- subset(joined_tables_all_basic, joined_tables_all_basic$productive.y == "TRUE" | joined_tables_all_basic$productive.y == "FALSE")
dim(joined_tables_all_basic_in_ground_truth)
head(joined_tables_all_basic_in_ground_truth)
#Perform the match
#Match V genes
v_genes <- joined_tables_all_basic_in_ground_truth %>% mutate(id = joined_tables_all_basic_in_ground_truth$sequence_id) %>% 
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
d_genes <- joined_tables_all_basic_in_ground_truth %>% mutate(id = joined_tables_all_basic_in_ground_truth$sequence_id) %>% 
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

j_genes <- joined_tables_all_basic_in_ground_truth %>% mutate(id = joined_tables_all_basic_in_ground_truth$sequence_id) %>% 
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


#Match productivity
df <- as.data.frame(cbind(joined_tables_all_basic_in_ground_truth$sequence_id,v_genes$v.gene.match,d_genes$d.gene.match,v_genes$v.gene.match))
head(df)
#select only the chains that were productive in the ground truth
df_productive <- df[df$V1 %in% Ground_Truth_split_subset$sequence_id,]
dim(df_productive)
Id <- df_productive$V1
df_productive$V1 <- as.character(df_productive$V1) 
df_productive$V2 <- as.logical(df_productive$V2)
df_productive$V3 <- as.logical(df_productive$V3)
df_productive$V4 <- as.logical(df_productive$V4)
df_productive <- 1*df_productive[,2:4]
dim(df_productive)
df_productive$sum <- rowSums(df_productive[,1:3])/3
Accuracy <- cbind(Id,df_productive)
dim(Accuracy)
Alll_Match <- subset(Accuracy, sum == 1)
value_all_match <- length(Alll_Match$Id)
Accuaracy_Value <- value_all_match/length(Accuracy$Id)
print(Accuaracy_Value)
print(value_all_match)
print(length(Accuracy$Id))


print(paste0("Number Heavy Chain reconstructed   ",length(dataset_LC_L002$sequence_id)))
print(paste0("Number Heavy Chain reconstructed productive   ",length(dataset_LC_L002_productive$sequence_id)))
print(paste0("Accuracy  ", value_all_match/length(Accuracy$Id)))


