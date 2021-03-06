#Load Ground Truth data Light Chain data
rm(list = ls())
KC <- read.table("KC.txt",header=T)
head(KC)
dim(KC)
#Exclude single end and low coverage samples which were not run
single_end_KC <- c('PW1-A1_KC', 'PW1-A4_KC',"PW1-A6_KC","PW1-B1_KC","PW1-B3_KC","PW1-B4_KC","PW1-B6_KC","PW1-C1_KC","PW1-C4_KC","PW1-C5_KC","PW2_A5_KC","PW2_A6_KC","PW2_A9_KC","PW2_B12_KC","PW2_C5_KC","PW2_D11_KC","PW2_D12_KC","PW2_D7_KC","PW2_E12_KC","PW2_G10_KC","PW2_G12_KC","PW2_G4_KC","PW2_G5_KC","PW2_G6_KC","PW2_G7_KC","PW2_H2_KC","PW2_H3_KC","PW2_H4_KC","PW2_H5_KC","PW3_B11_KC","PW3_B9_KC","PW3_C1_KC","PW3_C10_KC","PW3_C11_KC","PW3_C12_KC","PW3_C2_KC","PW3_C3_KC","PW3_C4_KC","PW3_C5_KC","PW3_C7_KC","PW3_C9_KC","PW3_D1_KC","PW3_D2_KC","PW3_D4_KC","PW3_D5_KC","PW3_D6_KC","PW3_F1_KC","PW3_F10_KC","PW3_F12_KC","PW3_F8_KC","PW3_F9_KC","PW3_G12_KC","PW3_G6_KC","PW3_G8_KC","PW3_G9_KC","PW3_H1_KC","PW3_H10_KC","PW3_H2_KC","PW3_H4_KC","PW3_H5_KC","PW3_H6_KC","PW3_H7_KC","PW3_H8_KC","PW3_H9_KC")
KC_no_single_end <- KC[ !grepl(paste(single_end_KC, collapse="|"), KC$sequence_id),]
test_KC_no_single_end <- KC_no_single_end %>% separate(sequence_id, c("A","B","C"))
test_KC_no_single_end <- test_KC_no_single_end %>% unite("Id", A:B, remove = TRUE)
Ground_Truth_light_chain_K <- subset(test_KC_no_single_end, select = c(1,3,4,5,6))
dim(Ground_Truth_light_chain_K)
Ground_Truth_Productive <- subset(Ground_Truth_light_chain_K,productive=="TRUE")
dim(Ground_Truth_Productive)
print(length(Ground_Truth_light_chain_K$Id))
print(length(Ground_Truth_Productive$Id)/length(Ground_Truth_light_chain_K$Id))
dim(Ground_Truth_light_chain_K)
head(Ground_Truth_light_chain_K)

#Load Basic Data
#coverage 1mln & 250k
dataset_LC_L002 <- read.table("Results_IGKL_from_BALDR.tsv",sep="\t",header=F)
dataset_LC_L002_IGL_Productive <- subset(dataset_LC_L002,V50=="VK" & V4==1)
dataset_LC_L002_subset <- subset(dataset_LC_L002_IGL_Productive, select = c(1,53,48,49))
head(dataset_LC_L002_subset)
dim(dataset_LC_L002_subset)
###count number of Light Chain Kappa and how many are productive
dataset_LC_L002_subset_productive <- subset(dataset_LC_L002_subset, V53 == "Yes")
length(dataset_LC_L002_subset_productive$V1)/length(dataset_LC_L002_subset$V1)
dataset_LC_L002_subset_productive_L <- subset(dataset_LC_L002_subset_productive, select = c(1,2,3,4))
head(dataset_LC_L002_subset_productive_L)

##Parse the name and subset the chains that are both in ground truth and assembled
test <- dataset_LC_L002_subset_productive_L %>% separate(V1, c("A","B","C","D","E"))
head(test)
test_LC_L002_subset <- test %>% unite("Id", A:B, remove = TRUE)
dim(test_LC_L002_subset)
A <- intersect(as.character(Ground_Truth_light_chain_K$Id),as.character(test_LC_L002_subset$Id))
length(A)
print(paste0("Number Light Chain Kappa Reconstructed ",length(A)))
LC_L002_BASIC_1mln_250L <- subset(test_LC_L002_subset, select = c(1,5,6,7))

LC_L002_BASIC_1mln_250L_productive_in_ground_truth <- LC_L002_BASIC_1mln_250L[ grepl(paste(Ground_Truth_light_chain_K$Id, collapse="|"), LC_L002_BASIC_1mln_250L$Id),]
head(LC_L002_BASIC_1mln_250L_productive_in_ground_truth)

#Rename the columns of the files
names(LC_L002_BASIC_1mln_250L_productive_in_ground_truth)[names(LC_L002_BASIC_1mln_250L_productive_in_ground_truth) == "V53"] <- "productive"
names(LC_L002_BASIC_1mln_250L_productive_in_ground_truth)[names(LC_L002_BASIC_1mln_250L_productive_in_ground_truth) == "V48"] <- "v_call"
names(LC_L002_BASIC_1mln_250L_productive_in_ground_truth)[names(LC_L002_BASIC_1mln_250L_productive_in_ground_truth) == "V49"] <- "j_call"
LC_L002_BASIC_1mln_250L_productive_in_ground_truth$productive <- gsub("Yes","TRUE",LC_L002_BASIC_1mln_250L_productive_in_ground_truth$productive)
head(LC_L002_BASIC_1mln_250L_productive_in_ground_truth)

#Join Table basic 1mlnk 100bp
joined_tables_all_basic <- full_join(LC_L002_BASIC_1mln_250L_productive_in_ground_truth,Ground_Truth_light_chain_K,by = c("Id"),copy = TRUE)
joined_tables_all_basic_in_ground_truth <- subset(joined_tables_all_basic, joined_tables_all_basic$productive.y == "TRUE" | joined_tables_all_basic$productive.y == "FALSE")
dim(joined_tables_all_basic_in_ground_truth)
joined_tables_all_basic_in_ground_truth
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
#Perform the match considering multiple genes in the output of the tool or the sanger for every V-J
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
dim(Ground_Truth_light_chain_K)
head(Ground_Truth_light_chain_K)
tot_chain_k <- 62
tot_chains_productive <- Ground_Truth_light_chain_K[which(Ground_Truth_light_chain_K$productive == "TRUE" ),]
length(tot_chains_productive$Id)
percentage_productive <- length(tot_chains_productive$Id)/tot_chain_k
print(percentage_productive)
head(Ground_Truth_light_chain_K)

#Match productivity
j_genes_k <- subset(j_genes,Chain == "KC")
head(j_genes_k)
match_productivity <- j_genes_k$productive.x == j_genes_k$productive.y
match_productivity$productivity.match <- match_productivity
head(match_productivity)
dim(j_genes_k)

v_genes_k <- subset(v_genes,Chain == "KC")
head(v_genes_k)
match_productivity <- v_genes_k$productive.x == v_genes_k$productive.y
match_productivity$productivity.match <- match_productivity
length(match_productivity)
dim(v_genes_k)

joined_tables_all_basic_k <- subset(joined_tables_all_basic,Chain == "KC")
head(joined_tables_all_basic_k)
dim(joined_tables_all_basic_k)

df <- as.data.frame(cbind(joined_tables_all_basic_k$Id,v_genes_k$match,j_genes_k$match))
head(df)
productive_k <- length(df$V1)
productive_k
head(productive_k)
list_lightk_id <- as.data.frame(Ground_Truth_light_chain_K$Id)
head(Ground_Truth_light_chain_K)
dim(Ground_Truth_light_chain_K)


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
Sensitivity <- cbind(Id,df_productive)
Alll_Match <- subset(Sensitivity, sum == 1)
value_all_match <- length(Alll_Match$Id)
Sensitivity_Value <- value_all_match/length(Sensitivity$Id)
print(Sensitivity_Value)
print(value_all_match)
print(length(Sensitivity$Id))



print(paste0("Percentage Reconstructed Kappa  ",length(dataset_LC_L002_subset$V1)))
print(paste0("Percentage Light Chain Kappa Productive ",length(dataset_LC_L002_subset_productive$V1)))      
print(paste0("Percentage Productive in The Total Assembled  ",length(dataset_LC_L002_subset_productive$V1)/length(dataset_LC_L002_subset$V1)))
print(paste0("Sensitivity  " ,value_all_match/length(Accuracy$Id)))

