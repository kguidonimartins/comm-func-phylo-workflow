###############################################################################
################### PHYLOCOM WORKFLOW #########################################
###############################################################################

# Uma pequena adaptação baseada em: 
# https://greggilbertlab.sites.ucsc.edu/wp-content/uploads/sites/276/2015/10/R_Class10b_PicatePhylomatic.pdf

temp <- tempfile()
download.file("http://phylodiversity.net/phylocom/phylocom-4.2.zip", temp, mode = "wb")
unzip(temp)

setwd("phylocom-4.2/w32")
getwd()

write.tree(get.tree.phylomatic, "tree-phylomatic")

file.rename(from = "../example_data/bladj_example/wikstrom.ages",
           to = "ages")

dated.tree <- system("phylocom bladj -f tree-phylomatic > dated-tree",
                     intern = TRUE)

write(dated.tree, "dated-tree")

dated.clean <- scan("dated-tree", what = "character", nmax = -1, sep = "]")

dated.clean <-
  sub(pattern = ")euphyllophyte:1.000000",
      replacement = "",
      x = dated.clean)

dated.clean <- sub(pattern = "\\(",
                   replacement = "",
                   x = dated.clean)

dated.clean <- gsub(pattern = " ",
                    replacement = "",
                    x = dated.clean)

dated.clean <- read.newick(text = dated.clean)

setwd("../../")
getwd()

write.tree(dated.clean, "data/dated.clean.txt")
