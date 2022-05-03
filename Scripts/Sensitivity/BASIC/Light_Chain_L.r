#Load Ground Truth data Light Chain data
rm(list = ls())
LC <- read.table("LC.txt",header=T)
head(LC)
dim(LC)

#Exclude single end and low coverage samples which were not run
single_end_LC <- c('PW1-A1_LC', 'PW1-A4_LC',"PW1-A6_LC","PW1-B1_LC","PW1-B3_LC","PW1-B4_LC","PW1-B6_LC","PW1-C1_LC","PW1-C4_LC","PW1-C5_LC","PW2_A5_LC","PW2_A6_LC","PW2_A9_LC","PW2_B12_LC","PW2_C5_LC","PW2_D11_LC","PW2_D12_LC","PW2_D7_LC","PW2_E12_LC","PW2_G10_LC","PW2_G12_LC","PW2_G4_LC","PW2_G5_LC","PW2_G6_LC","PW2_G7_LC","PW2_H2_LC","PW2_H3_LC","PW2_H4_LC","PW2_H5_LC","PW3_B11_LC","PW3_B9_LC","PW3_C1_LC","PW3_C10_LC","PW3_C11_LC","PW3_C12_LC","PW3_C2_LC","PW3_C3_LC","PW3_C4_LC","PW3_C5_LC","PW3_C7_LC","PW3_C9_LC","PW3_D1_LC","PW3_D2_LC","PW3_D4_LC","PW3_D5_LC","PW3_D6_LC","PW3_F1_LC","PW3_F10_LC","PW3_F12_LC","PW3_F8_LC","PW3_F9_LC","PW3_G12_LC","PW3_G6_LC","PW3_G8_LC","PW3_G9_LC","PW3_H1_LC","PW3_H10_LC","PW3_H2_LC","PW3_H4_LC","PW3_H5_LC","PW3_H6_LC","PW3_H7_LC","PW3_H8_LC","PW3_H9_LC")
LC_no_single_end <- LC[ !grepl(paste(single_end_LC, collapse="|"), LC$sequence_id),]
test_LC_no_single_end <- LC_no_single_end %>% separate(sequence_id, c("A","B","C"))
head(test_LC_no_single_end)
test_LC_no_single_end <- test_LC_no_single_end %>% unite("Id", A:B, remove = TRUE)
Ground_Truth_light_chain_L <- subset(test_LC_no_single_end, select = c(1,3,4,5,6))
print(length(Ground_Truth_light_chain_L$Id))
Ground_Truth_Productive <- subset(Ground_Truth_light_chain_L,productive=="TRUE")
print(length(Ground_Truth_light_chain_L$Id))
print(length(Ground_Truth_Productive$Id)/length(Ground_Truth_light_chain_L$Id))

#Load Basic Data
#coverage 1mln & 250L
#Load Basic Data
#coverage 1mln & 250k
dataset_LC_L002 <- read.table("Results_IGKL_from_BALDR.tsv",sep="\t",header=F)
dataset_LC_L002_IGL_Productive <- subset(dataset_LC_L002,V50=="VL" & V4==1)
dataset_LC_L002_subset <- subset(dataset_LC_L002_IGL_Productive, select = c(1,53,48,49))
head(dataset_LC_L002_subset)
dim(dataset_LC_L002_subset)
###count number of Light Chain Kappa and how many are productive
dataset_LC_L002_subset_productive <- subset(dataset_LC_L002_subset, V53 == "Yes")
length(dataset_LC_L002_subset_productive$V1)/length(dataset_LC_L002_subset$V1)
dataset_LC_L002_subset_productive_L <- subset(dataset_LC_L002_subset_productive, select = c(1,2,3,4))
head(dataset_LC_L002_subset_productive_L)

##Parse the name and subset the chains that are both in ground truth and assembled
test <- dataset_LC_L002_subset_productive_L %>% separate(sequence_id, c("A","B","C","D","E"))
test_LC_L002_subset <- test %>% unite("Id", C:D, remove = TRUE)
head(test_LC_L002_subset)
LC_L002_BASIC_1mln_250L <- subset(test_LC_L002_subset, select = c(3,5,6,7,8))
head(LC_L002_BASIC_1mln_250L)
dim(LC_L002_BASIC_1mln_250L)
LC_L002_BASIC_1mln_250L_productive_in_ground_truth <- LC_L002_BASIC_1mln_250L[ grepl(paste(Ground_Truth_light_chain_L$Id, collapse="|"), LC_L002_BASIC_1mln_250L$Id),]
head(LC_L002_BASIC_1mln_250L_productive_in_ground_truth)
dim(LC_L002_BASIC_1mln_250L_productive_in_ground_truth)

#Join Table basic 1mlnL 100bp
joined_tables_all_basic <- full_join(LC_L002_BASIC_1mln_250L_productive_in_ground_truth,Ground_Truth_light_chain_L,by = c("Id"),copy = TRUE)
joined_tables_all_basic_in_ground_truth <- subset(joined_tables_all_basic, joined_tables_all_basic$productive.y == "TRUE" | joined_tables_all_basic$productive.y == "FALSE")
dim(joined_tables_all_basic_in_ground_truth)

#Perform the match
#Match V genes
v_genes <- joined_tables_all_basic %>% mutate(id = joined_tables_all_basic$Id) %>% 
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
v_genes$v_call.y <- replace(v_genes$v_call.y, is.na(v_genes$v_call.y), FALSE)

names(v_genes)[10] <- "v.gene.match"
head(v_genes)



####################################################################
#Perform the match considering multiple genes in the output of the tool or the sanger for every V-D-J
###################################################################################

#Match J genes

j_genes <- joined_tables_all_basic %>% mutate(id = joined_tables_all_basic$Id) %>% 
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
j_genes$v_call.y <- replace(j_genes$v_call.y, is.na(j_genes$v_call.y), FALSE)


names(j_genes)[10] <- "j.gene.match"
j_genes

###################################################################################
#number of heavy chain in Ground Truth and number of chain in ground truth productive
dim(Ground_Truth_light_chain_L)
head(Ground_Truth_light_chain_L)
dim(Ground_Truth_light_chain_L)
tot_chain_L <- length(Ground_Truth_light_chain_L$Id)
tot_chains_productive <- Ground_Truth_light_chain_L[which(Ground_Truth_light_chain_L$productive == "TRUE" ),]
length(tot_chains_productive$Id)
percentage_productive <- length(tot_chains_productive$Id)/tot_chain_L
print(percentage_productive)
head(Ground_Truth_light_chain_L)

#Match productivity
j_genes_L <- subset(j_genes,Chain == "LC")
head(j_genes_L)
match_productivity <- j_genes_L$productive.x == j_genes_L$productive.y
match_productivity$productivity.match <- match_productivity
length(match_productivity)
dim(j_genes_L)

v_genes_L <- subset(v_genes,Chain == "LC")
head(v_genes_L)
match_productivity <- v_genes_L$productive.x == v_genes_L$productive.y
match_productivity$productivity.match <- match_productivity
length(match_productivity)
v_genes_L$productive.y

joined_tables_all_basic_L <- subset(joined_tables_all_basic,Chain == "LC")
head(joined_tables_all_basic_L)
dim(joined_tables_all_basic_L)

df <- as.data.frame(cbind(joined_tables_all_basic_L$Id,v_genes_L$v.gene.match,j_genes_L$j.gene.match))
head(df)
productive_L <- length(df$V1)
productive_L
list_lightL_id <- as.data.frame(Ground_Truth_light_chain_L$Id)
head(Ground_Truth_light_chain_L)
dim(Ground_Truth_light_chain_L)


#select only the chains that were productive in the ground truth
df_productive <- df[ grepl(paste(tot_chains_productive$Id, collapse="|"), df$V1),]
head(df_productive)
Id <- df_productive$V1
df_productive$V1 <- as.character(df_productive$V1) 
df_productive$V2 <- as.logical(df_productive$V2)
df_productive$V3 <- as.logical(df_productive$V3)
df_productive <- 1*df_productive[,2:3]
head(df_productive)
df_productive$sum <- rowSums(df_productive[,1:2])/2
Accuracy <- cbind(Id,df_productive)
Alll_Match <- subset(Accuracy, sum == 1)
value_all_match <- length(Alll_Match$Id)
Accuaracy_Value <- value_all_match/length(Accuracy$Id)
print(Accuaracy_Value)
print(value_all_match)



print(paste0("Percentage Productive in The Ground Truth  ", length(Ground_Truth_Productive$Id)/length(Ground_Truth_light_chain_L$Id)))
print(paste0("Percentage Reconstructed Lambda  ",length(dataset_LC_L002_subset$sequence_id)))
print(paste0("Percentage Light Chain Lambda Productive ",length(dataset_LC_L002_subset_productive$sequence_id)))      
print(paste0("Percentage Productive in The Total Assembled  ",length(dataset_LC_L002_subset_productive$sequence_id)/length(dataset_LC_L002_subset$sequence_id)))
print(paste0("Accuracy  " ,value_all_match/length(Accuracy$Id)))



