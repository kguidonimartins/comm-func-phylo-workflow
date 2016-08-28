# R-functions to interface with phylocom
# Cam Webb <cwebb@oeb.harvard.edu>
# and later Steve Kembel and David Ackerly


#### 1. Internal utility functions ####

# write out the sample file
writesample <- function(community, file="sample.tmp")
  {
    # turn a vegan-style community matrix into a condensed table
    stack.phylocom <- function(z)
      {
        temp <- data.frame(expand.grid(dimnames(z))[1:2], as.vector(as.matrix(z)))
        temp <- temp[(temp[,3]>0) & !is.na(temp[,3]), ]
        temp <- temp[sort.list(temp[,1]),]
        data.frame(temp[,1], temp[,3], temp[,2])
      }

    write.table(stack.phylocom(community), file = file,
                append = FALSE, sep = "\t",
                eol = "\n", quote=F, row.names = F, col.names = F)
  }

readtraits <- function(file = "traits")
  {
    x <- read.table(file=file, header=F, sep="\t")
    type <- as.numeric(as.matrix(x[1,-1]))
    cnames <- as.character(as.matrix(x[2,]))
    if (identical("name", cnames[1]))
      {
        traits <- x[c(-1,-2),-1]
        rownames(traits) <- x[c(-1,-2),1]
        colnames(traits) <- cnames[-1]
      }
    else
      {
        traits <- x[-1,-1]
        rownames(traits) <- x[-1,1]
        colnames(traits) <- paste("trait", 1:(ncol(x)-1), sep="")
      }
    for (i in 1:dim(traits)[2])
      {
        if (type[i] == 0) traits[,i] <- as.factor(traits[,i])
        else if (type[i] == 3) traits[,i] <- as.numeric(traits[,i])
        else stop("cannot handle trait types 1 or 2")
      }
    drop.levels <- function(dat){
    # Drop unused factor levels from all factors in a data.frame
    # Author: Kevin Wright.  Idea by Brian Ripley.
      dat[] <- lapply(dat, function(x) x[,drop=TRUE])
      return(dat)
    }
    return(drop.levels(traits))
  }

writetraits <- function(traits, file="traits.tmp", phylo=NULL)
  {
    # Reads a dataframe, where taxa names are row names, and columns are
    # either two-level factors, or continuous variables.  A phylo file can be
    # provided to check that the rows of the traits object are the same as the 
    # terminals on the phylo

    # Some tests
    if (!is.data.frame(traits)) stop("traits object is not a data frame")
    type <- rep(NA, ncol(traits))
    for (i in 1:ncol(traits))
      {
        if (is.factor(traits[,i]))
          {
            type[i] <- 0
            if (nlevels(traits[,i]) > 2) stop(paste("traits col", i,
                                                    "has too many levels"))
          }
        if (is.numeric(traits[,i])) type[i] <- 3
        if (is.na(type[i])) stop(paste("traits col", i,
                                      "is neither a factor nor numeric"))
      }

    if (!is.null(phylo))
      {
        library(ape)
        taxa <- sort(phylo$tip.label)
        if (!identical(taxa, sort(rownames(traits)))) stop(
              "traits taxon names are not the same as phylo terminals")
      }

#    system("rm -f traits.tmp >& /dev/null; echo \"type\t0\nname\ttrait0\" > traits.tmp")
    type <- paste(c("type", type), collapse="\t")
    write(type, file = file, append = FALSE)
    name <- paste(c("name", colnames(traits)), collapse="\t")
    write(name, file = file, append = TRUE)
    write.table(traits, file = file,
                sep="\t", append = TRUE, eol = "\n", quote=F,
                row.names = T, col.names = F)
  }


#### 2. I/O functions ####

# read in a sample file (default name = `sample') 
readsample <- function(file = "sample")
  {
    x <- read.table(file=file, header=F, sep="\t",
                    col.names=c("plot", "abund", "id"))
    y <- tapply(x$abund, list(x$plot, x$id), sum)
    y[is.na(y)] <- 0
    y
  }

# read in a phylo file (default name = `phylo') 
readphylo <- function(file = "phylo")
  {
    library(ape)
    x <- read.tree(file=file)
    x
  }

# draw a phylogeny from a newick file using the `ape' package.  Default file is `phylo' 
drawphylo <- function(file = "phylo")
  {
    library(ape)
    t <- read.tree(file = file)
    plot(t)
  }


#### 3. Community functions using vegan ####

# simple MDS, with my favourite parameters
drawmds <- function(community)
  {
    library(vegan)
    community[community > 0] <- 1 # turn into p/a
    com.mds <- metaMDS(community, distance="mountford")
    plot(com.mds, display="species", type="n")
    points(com.mds, display="species", col="red", cex=0.7)
    text(com.mds, display="sites", col="blue")
  }

# the vegemite function, gotta love it!
vegetab <- function(community)
  {
    library(vegan)
    dca <- decorana(community)
    vegemite(community, dca, scale="log")
  }


#### 4. Phylocom functions #### 

# run comdist on a community matrix, and optional phylogeny (file `phylo' is default)
phylocom.comdist <- function(community, phylo=NULL, method="mean")
  {
    if (!is.null(phylo))
      {
        library(ape)
        write.tree(phylo, file = "phylo.tmp", append = FALSE)
      }
    else
      {
        system("cp -f phylo phylo.tmp")
      }
    writesample(community)
    
    if (method == "nearest") {pc <- "comdistnn"}
    else {pc <- "comdist"}

    system(paste("phylocom", pc, "-s sample.tmp -f phylo.tmp > comdist.tmp"))
    cd <- read.table("comdist.tmp", header=T, row.names=1, sep="\t")
    system("rm -f comdist.tmp sample.tmp phylo.tmp")
    as.dist(cd)
  }

# run comstruct on a community matrix, and optional phylogeny (file `phylo' is default)
phylocom.comstruct <- function(community, swapmethod = 2,
                               runs = 999, swaps = 1000, phylo=NULL)
  {
    if (!is.null(phylo))
      {
        library(ape)
        write.tree(phylo, file = "phylo.tmp", append = FALSE)
      }
    else
      {
        system("cp -f phylo phylo.tmp")
      }
    writesample(community)
    command <- paste("phylocom comstruct -r", runs, "-m", swapmethod,
                     "-w", swaps, "-s sample.tmp -f phylo.tmp > comstruct.tmp")
    system(command)
    cs <- read.table("comstruct.tmp", header=T, row.names=1, skip=1, sep="\t")
    system("rm -f comstruct.tmp sample.tmp phylo.tmp")
    #structure(list(nri = cs$NRI, nti = cs$NTI, mpd = cs$MPD, mntd = cs$MNND,
    #               all = cs), class="phylocom.construct")
    cs
  }

# run pd on a community matrix, and optional phylogeny (file `phylo' is default)
phylocom.pd <- function(community, phylo=NULL)
  {
    if (!is.null(phylo))
      {
        library(ape)
        write.tree(phylo, file = "phylo.tmp", append = FALSE)
      }
    else
      {
        system("cp -f phylo phylo.tmp")
      }
    writesample(community)
    command <- "phylocom pd -s sample.tmp -f phylo.tmp > pd.tmp"
    system(command)
    pd <- read.table("pd.tmp", header=T, row.names=1)
    system("rm -f pd.tmp sample.tmp phylo.tmp")
    pd$PD
  }

# Run phylocom node-as-factor (needs phylocom v4.0)
phylocom.naf <- function(community, phylo=NULL)
 {
    if (!is.null(phylo))
      {
        library(ape)
        write.tree(phylo, file = "phylo.tmp", append = FALSE)
        taxa <- phylo$tip.label
      }
    else
      {
        library(ape)
        system("cp -f phylo phylo.tmp")
        taxa <- read.tree("phylo")$tip.label
      }
#    # make a dummy traits file (this should be fixed in phylocom so that a traits file is not needed
    system("rm -f traits.tmp >& /dev/null; echo \"type\t0\nname\ttrait0\" > traits.tmp")
    write.table(data.frame(taxa, rep(0, length(taxa))),
                file = "traits.tmp", append = TRUE, sep = "\t",
                eol = "\n", quote=F, row.names = F, col.names = F)
    writesample(community)
    system("phylocom naf -f phylo.tmp -t traits.tmp -s sample.tmp > naf.tmp")

    naf <- read.table("naf.tmp", header=T, row.names=1, skip=0,
                      sep="\t", na.strings = ".")
    system("rm -f naf.tmp phylo.tmp sample.tmp traits.tmp")
    #structure(list(nri = cs$NRI, nti = cs$NTI, mpd = cs$MPD, mntd = cs$MNND,
    #               all = cs), class="phylocom.construct")
    naf <- naf[, (substr(colnames(naf),1,1)=="N")]
    colnames(naf) <- substr(colnames(naf),3,1000)
    as.matrix(naf)
  }


#### 5. Phylordination functions ####

# clustering based on phylogenetic distance
phyclust <- function(community, clustmethod="ward", phydist = "mean",
                     phylo=NULL)
  {
    x <- phylocom.comdist(community, phylo=phylo, method=phydist)
    plot(hclust(x, method=clustmethod),
         sub=paste("Cluster method is", clustmethod),
         main="Clustering based on inter-sample phylogenetic distance",
         xlab="")
  }

# ordination based on phylogenetic distance (cf. drawmds)
phylord <- function(community, phydist = "mean", phylo=NULL)
  {
    library(MASS)
    x <- phylocom.comdist(community, phylo=phylo, method=phydist)
    x.mds <- isoMDS(x)
    plot(x.mds$points, type = "n", xlab="MDS 1", ylab="MDS 2", main="NMDS based on inter-sample phylogenetic distance")
    text(x.mds$points, labels = rownames(community))
  }

# ordination based on internal node structure (phylocom naf)
phylord.naf <- function(community, phylo=NULL)
  {
    naf <- phylocom.naf(community, phylo)
    naf2 <- matrix(data = NA, nrow = nrow(naf), ncol=ncol(naf))
    rownames(naf2) <- rownames(naf)
    colnames(naf2) <- colnames(naf)
    naf2[!is.na(naf)] <- 1
    naf2[is.na(naf)] <- 0
    # naf2

    # this is an ugly loop method - please find elegant R-ish method if you can - but note that the species lists need not match up between the naf and community matrices
    comnaf <- matrix(data = 0, nrow = nrow(community), ncol = ncol(naf2))
    rownames(comnaf) <- rownames(community)
    colnames(comnaf) <- colnames(naf2)
    for (i in 1:nrow(community))
      {
        for (j in 1:ncol(community))
          {
            if (community[i,j] > 0)
              {
                for (k in 1:nrow(naf2))
                  {
                    if( colnames(community)[j] == rownames(naf2)[k])
                      {
                        comnaf[i,] = comnaf[i,] + naf2[k,]
                      }
                  }
              }
          }
      }
    # drop any 0 nodes
    comnaf <- comnaf[,colSums(comnaf) > 0]

    library(vegan)
    comnaf.mds <- metaMDS(comnaf)
    plot(comnaf.mds, display="species", type="n")
    points(comnaf.mds, display="species", col="red", cex=0.7)
    text(comnaf.mds, display="sites", col="blue")

    structure(list(comnaf = comnaf))

  }


