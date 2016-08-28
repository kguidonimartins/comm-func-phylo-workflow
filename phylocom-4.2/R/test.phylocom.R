# Test suite for phylocom.R
# Cam Webb <cwebb@oeb.harvard.edu>


# read the utility file:
source("phylocom.R")
par(ask=TRUE)

# read sample file, default is `sample' (in the active directory)
readsample()
a <- readsample()

# read the traits file, default is `traits' (in the active directory)
p <- readphylo()
t <- readtraits()
writetraits(t, "traits.out", phylo=p)

# draw the phylo file, default is `phylo' (in the active directory)
drawphylo()

# do a simple mds on the community data
drawmds(a)

# vegan table on the community data
vegetab(a)

# phylocom.  Phylocom must be either in the active directory (and the path includes the . directory) or within the path.  I.e. will not run if you still need to type ./phylocom.  If you still are not going to fix your paths (basic UNIX), you'll have to edit the phylocom.R file so that every occurence of phylocom is changed to ./phylocom

phylocom.comstruct(a)

#more complex - a phylo file with branch lengths

system("phylocom bladj > phylo.bl")
p <- readphylo("phylo.bl")

phylocom.comstruct(a, phylo=p, runs=50)$MPD

phylocom.comdist(a, phylo=p, method="nearest")

phylocom.pd(a)

# phylordination

phylord(a, phylo=p)

phylord.naf(a, phylo = p)

phyclust(a, phylo=p)

# using the nearest taxon phy dist
phyclust(a, phylo=p, phydist="nearest")



