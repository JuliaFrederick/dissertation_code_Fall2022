## Tree manipulation for JCF Ixodes ##

## Install all necessary packages ##

# install.packages('ggplot2')
# install.packages('ggrepel')
# install.packages('phangorn')
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("ggtree")
# BiocManager::install("Biostrings")
## BiocManager provides tools for the analysis and comprehension of high-throughput genomic data. 
## Bioconductor uses the R statistical programming language, and is open source and open 
## development. Use the 'biocLite()' script to install Bioconductor packages.

## Load all necessary packages ##
library('ape')
library('ggplot2')
library('ggrepel')
library('Biostrings')
library('ggtree')
library('phangorn')

## GETTING STARTED ##
Iscap <- read.tree("./data/Tree_alignment_v5_realigned_modified_consensus_tree.newick")
#Iscap <- read.tree("/Users/alecthompson/Desktop/Tree\ alignment\ v5\ -\ realigned\ modified\ consensus\ tree.newick") 
  ## load tree into R, tree format has to be either Newick or Nexus.

Iscap[["tip.label"]]
str(Iscap)
  # tree$edge provides a matrix of node and tip numbers. It is numbered as follows:
    # tips: 1 to the total number of ‘n’ tips (extant species)
    # nodes: ‘n+1’ to the total number of nodes  
  # tree$Nnode indicates that the number of internal nodes in this case NULL (tree is unrooted)\
  # tree$tip.label names of taxa at the tips
  # tree$edge.length branch length values

is.rooted(Iscap)
outgroup <- c("L34302_-_Rhipicephalus_sanguineus") 
  # depending on what is listed an entire clade can be an outgroup (use tip labels)
rooted_tree <- root(Iscap, outgroup, resolve.root = TRUE)
is.rooted(rooted_tree)

plot(rooted_tree, type = 'phylogram', label.offset = 0.005, edge.width = 2) # change "type" to alter format of tree
nodelabels() # shows tip labels in blue
tiplabels() # shows node labels in yellow

tree.x <- ggtree(rooted_tree)

tree.x + geom_text(aes(label=node), hjust=-.3)
## Label the nodes so we can play around with the format of the tree.

tree.x + geom_tiplab(size=5, label = gsub("_", " ", tree.x$data$label
                                        [which(tree.x$data$isTip == TRUE)]), fontface="italic")
## Add the labels to tips of the tree and format them to remove underscores

bootstrap <- tree.x$data
bootstrap <- bootstrap[!bootstrap$isTip,]
bootstrap$label <- as.numeric(bootstrap$label)
bootstrap$label <- round(bootstrap$label)
bootstrap <- bootstrap[bootstrap$label > 65,]

tree.x + geom_label_repel(data=bootstrap, aes(label=label), size = 3, check_overlap = TRUE)
## Get the bootstrap data separated from the node data

## Put it all together
tree.y <- tree.x + 
  geom_tiplab(size=5, aes(subset = (node %in% c(18,20,22:24,27)),
                          label = gsub("_", " ", tree.x$data$label), size = 3), fontface = 4) +
  geom_tiplab(size=5, aes(subset = (node %in% c(1:17,19,21,25:26,28:30)),
                          label = gsub("_", " ", tree.x$data$label), size = 3), fontface = 3) +
  geom_cladelabel(node = 45, label = "Southern Clade", align = TRUE, color = 'red') +
  geom_cladelabel(node = 41, label = "American Clade", align = TRUE, color= 'blue') +
  # geom_treescale(x = 0, y = 0, width = 0.05, fontsize = 3) + 
  geom_label_repel(data = bootstrap, aes(label = label), size = 4, nudge_x = -0.0016) +
  coord_cartesian(clip = "off") +
  theme_tree(plot.margin = margin(0, 200, 0, 0))
tree.y

ggsave(file = "assets/test.pdf", tree.y, width = 50, height = 50, units = "cm", limitsize = FALSE)


