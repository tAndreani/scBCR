R
library(immuneSIM)

# first obtain BCR fasta sequences with different amount of mutations in the variable domains (15, 30, 45, 60).
# here an example for IGK BCR with 15 somatic hypermutations (SHM):

sim_repertoire_IGK_shm15 <- immuneSIM(number_of_seqs = 100,
                                        species = "hs",
                                        receptor = "ig",
                                        chain = "k",
                                        shm.mode="data",
                                        shm.prob=15/350)
                     
write.table(sim_repertoire_IGK_shm15,"igk_shm15.tsv",quote=F,col.names=T,row.names=F,sep="\t")


sim_repertoire_IGK_shm30 <- immuneSIM(number_of_seqs = 100,
                                        species = "hs",
                                        receptor = "ig",
                                        chain = "k",
                                        shm.mode="data",
                                        shm.prob=30/350)
                     
write.table(sim_repertoire_IGK_shm30,"igk_shm30.tsv",quote=F,col.names=T,row.names=F,sep="\t")

sim_repertoire_IGK_shm45 <- immuneSIM(number_of_seqs = 100,
                                        species = "hs",
                                        receptor = "ig",
                                        chain = "k",
                                        shm.mode="data",
                                        shm.prob=45/350)
                     
write.table(sim_repertoire_IGK_shm45,"igk_shm45.tsv",quote=F,col.names=T,row.names=F,sep="\t")

sim_repertoire_IGK_shm60 <- immuneSIM(number_of_seqs = 100,
                                        species = "hs",
                                        receptor = "ig",
                                        chain = "k",
                                        shm.mode="data",
                                        shm.prob=60/350)
                     
write.table(sim_repertoire_IGK_shm60,"igk_shm60.tsv",quote=F,col.names=T,row.names=F,sep="\t")
q()
