###############################################################################
################### PHYLOCOM WORKFLOW #########################################
###############################################################################

# Uma pequena adaptação baseada em: 
# https://greggilbertlab.sites.ucsc.edu/wp-content/uploads/sites/276/2015/10/R_Class10b_PicatePhylomatic.pdf

temp <- tempfile()
download.file("https://github.com/phylocom/phylocom/releases/download/4.2/phylocom-4.2.zip", temp, method = "wget")
unzip(temp)

setwd("phylocom-4.2/src")
system("make", invisible = FALSE)

write.tree(get.tree.phylomatic, "tree-phylomatic")

file.copy(from = "../example_data/bladj_example/wikstrom.ages", to = ".")

file.rename(from = "wikstrom.ages", to = "ages")

dated.tree <- system("./phylocom bladj -f tree-phylomatic > dated-tree")

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

write.tree(dated.clean, "data/dated.clean.txt")
