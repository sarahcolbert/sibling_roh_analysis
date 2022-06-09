library(data.table)
library(igraph)
library(plyr)

# Read the kinship file produced by KING
df_kinship <- fread("/data/king.kin0")

# Ensure IDs are characters not integers
df_kinship$ID1 <- as.character(df_kinship$ID1)
df_kinship$ID2 <- as.character(df_kinship$ID2)

# Find and select pairs with kinship consistent with full siblings
iSibs <- which(((df_kinship$Kinship>0.17)&(df_kinship$IBS0>0.001))|(df_kinship$Kinship>0.4))  # Note this includes identical twins
df_sibs <- df_kinship[iSibs,]

# Use igraph to find groups where all members appear to be siblings
g <- graph.data.frame(df_sibs[,c("ID1","ID2")],directed = F)
cliques <- max_cliques(g, min=2)

# make a data table where family is recorded for each ID
df_family <- data.table(id=names(unlist(cliques)),sib_family=rep(paste0("f", seq_along(cliques)), sapply(cliques, length)))
# and use this data table to add sib_family to each sibling pair
setkey(df_family,id)
df_sibs$sib_family <- unlist(df_family[df_sibs$ID1,"sib_family"])

# to have greater confidence that groups are truely all siblings I only keep families in which
# the mean kinship between all members is greater than 0.19 - in UKB this only removes ~30 families of size 2.
# Obviously, a similar result could be achieved by changing the threshold in line 13 above, but this removes a few individuals from larger
# families who happen to have lower kinship (between 0.17 and 0.19) with just one of their siblings.

# calculate the mean kinships of all members for each family
df_lookup <- as.data.table(ddply(df_sibs, .(sib_family), summarise,  mean_kinship=mean(Kinship) ))
# and only keep individuals in families with mean kinship > 0.19
i <- which(df_lookup$mean_kinship<0.19)
j <- which(df_family$sib_family %in% df_lookup$sib_family[i])
if(length(j)>0){
	df_family <- df_family[-j,]
}

# order by sib_family
df_family <- df_family[order(df_family$sib_family),]

# finally write out the family of each sibling
write.table(df_family,"/data/cliques_17_19.tsv",col.names = T, row.names = F, sep = "\t",quote = F)
