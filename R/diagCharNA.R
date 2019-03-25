#' Identification of diagnostic molecular characters in nucleotide alignments
#' for the delineation of taxa
#'
#' This function is a tool for an automated identification of diagnostic
#' molecular characters that allow to distinguish taxa within a nucleotide
#' alignment. For each taxon given in \code{taxOfInt}, it identifies the
#' diagnostic characters and returns their alignment positions, their types,
#' the states that are characteristic for the taxon of interest and
#' (in case of type 4 characters) the taxon that it was compared with.
#'
#' @param DNAbin An object (the nucleotide alignment) of class 'DNAbin'.
#' @param taxVector The taxon vector. Default assumes that each row in the
#' alignment belongs to a different taxon.
#' @param taxOfInt A vector containing the taxa for which diagnostic
#' molecular characters shall be extracted. Default is "all".
#' @param types A vector containing the types of diagnostic molecular
#' characters that shall be extracted. The types can be "all" or any
#' combination of "type1", "type2", "type3" and "type4". Default is "type1",
#' "type 2" and "type3".
#' @param gapValid Boolean variable denoting if a gap can be a characteristic
#' state for taxon i (and taxon l in case of type 4). Default is \code{TRUE}.
#'
#' @return \code{diagCharNA} returns a list, where each entry belongs to one
#' taxon of interest. Each taxon of interest has a set of diagnostic
#' molecular characters (position, type, characteristic states for taxon of
#' interest, compared taxa) assigned to it.
#' \code{type1} means that the character is suitable to distinguish each
#' individual of the taxon of interest from all individuals of the remaining
#' taxa, and that it is fixed for one state in the taxon of interest.
#' \code{type2} means that the character is suitable to distinguish each
#' individual of the taxon of interest from all individuals of the remaining
#' taxa, and that it is not fixed for one state in the taxon of interest.
#' \code{type3} means that the character is suitable to distinguish some (but
#' not all) individuals of the taxon of interest from all individuals of the
#' remaining taxa.
#' \code{type4} means that the character is suitable to distinguish each
#' individual of the taxon of interest from all individuals of at least one
#' (but not all) other taxon while being fixed in both the taxon of interest
#' and the compared taxa.
#'
#' \code{diagCharNA} returns for each taxon of interest the following elements:
#' \item{position}{The positions of its diagnostic molecular characters.}
#' \item{type}{The types of the diagnostic molecular characters.}
#' \item{states}{The states that are characteristic for the taxon i of interest,
#' i.e. states that are distinct from "n" and unique to the taxon of interest
#' (in case of type 1, 2 or 3), or fixed in the taxon of interest (type 4).}
#' \item{compared taxa}{Only relevant for type 4 characters. It
#' contains the name x if the character is found to be a type 4 character of
#' the taxon of interest when being compared to taxon x.}
#'
#' @author A. Luise Kuehn <luise.kuehn@@uni-greifswald.de>
#'
#' @references Kuehn, A.L., Haase, M. 2019. QUIDDICH: QUick IDentification of
#' DIagnostic CHaracters.

#' @examples
#' #using a dataset from spider
#' #install.packages("spider")
#' library(spider)
#' data("anoteropsis")
#' anoTax <- sapply(strsplit(dimnames(anoteropsis)[[1]], split="_"),
#' 	function(x) paste(x[1], x[2], sep="_"))
#' diagCharNA(anoteropsis, anoTax, taxOfInt="all")
#' diagCharNA(anoteropsis, anoTax, taxOfInt="all", types=c("type1","type2"))
#' #
#' #with loading of a fasta file
#' #install.packages("adegenet")
#' library(adegenet)
#' alignment <- fasta2DNAbin(paste0(find.package("quiddich"), "/extData/example.fasta"))
#' taxonVector <- as.vector(sapply(dimnames(alignment)[[1]], function(x) substr(x,1,4)))
#' diagCharNA(alignment, taxonVector)

#' @importFrom ape seg.sites

#' @export
diagCharNA <- function(DNAbin, taxVector=dimnames(DNAbin)[[1]], taxOfInt="all", types=c("type1","type2","type3"), gapValid=TRUE)
{
  # transform dataset (alignment) into matrix
  DNAbin <- as.matrix(DNAbin)
  # extract polymorphic sites in the alignment (positions of interest)
  inform <- ape::seg.sites(DNAbin)

  # specify the taxa for which diagnostic characters shall be extracted (taxa of interest)
  if(taxOfInt[1]=="all")
  {
    taxOfInt <- taxVector
  }
  # for each taxon: store the assignment (taxon number -> sequence numbers)
  taxSeqs <- lapply(unique(taxVector), function(x) which(taxVector==x))
  # for each taxon of interest: store at which entry in taxSeq its sequences can be found
  taxSeqsPos <- sapply(unique(taxOfInt), function(x) which(unique(taxVector)==x))

  # specify the types of diagnostic molecular characters that shall be extracted
  if(types[1]=="all")
  {
    types <- c("type1","type2","type3","type4")
  }
  # give a warning if type 4 characters shall be extracted and the dataset
  # contains singletons
  if("type4" %in% types && 1 %in% table(taxVector))
  {
    writeLines("WARNING: \n   The dataset contains singletons! \n   An extraction of type 4 characters might not be sensible! \n...continuing calculations...")
  }

  # result stores the diagnostic characters for each taxon of interest
  result <- list()

  # consider all taxa of interest successively
  for(i in taxSeqsPos)
  {
    # newtax stores the diagnostic characters for taxon i
    newtax <- matrix(c("position","type","states","compared taxa"), ncol=4)
    # go through all positions of interest
    for(j in 1:length(inform))
    {
      # extract states(i,j) and states(rest,j) from the alignment
      states_i <- unique(as.character(DNAbin[taxSeqs[[i]],inform[j]]))
      states_rest <- unique(as.character(DNAbin[-taxSeqs[[i]],inform[j]]))

      # check if j is diagnostic for taxon i
      # newpos = (position, type, states that make it diagnostic, compared taxa)
      newpos <- vector(mode="character",length=4)
      newpos[1] <- inform[j]
      newpos[2] <- "not type1/2/3"
      # for j to be diagnostic (type1/2/3), states(rest,j) cannot contain masked entries
      if(!("n" %in% states_rest))
      {
        for(k in 1:length(states_i))
        {
          a <- states_i[k]
          # for j to be of type1/2/3, state a cannot be masked or in states(rest,j)
          # (or a gap if gaps are not valid states)
          if(a != "n" && !(a %in% states_rest) && isTRUE(gapValid | a != "-"))
          {
            newpos[2] <- "type3"
            newpos[3] <- paste(newpos[3], a, sep="")
            # for j to be of type1/2, states(i,j) and states(rest,j) cannot share any states
            # and states(i,j) cannot contain N
            if(!(any(states_i %in% states_rest)) && !("n" %in% states_i))
            {
              newpos[2] <- "type2"
              # for j to be of type1, a must be the only state in states(i,j)
              if(length(states_i)==1)
              {
                newpos[2] <- "type1"
              }
            }
          }
        }
      }
      # if j is identified as type1/2/3 and the type is wanted -> store position j
      if(newpos[2] %in% types)
      {
        newtax <- rbind(newtax, newpos)
      }
      # if j is not identified as type1/2/3, it might be type 4
      else if("type4" %in% types && newpos[2] == "not type1/2/3")
      {
        # for j to be of type4, it must be fixed in taxon i for any state except "n"
        # (and "-" if gaps are not valid states)
        if(length(states_i)==1 && states_i[1]!="n" && isTRUE(gapValid | states_i[1]!="-"))
        {
          # go through all comparable taxa l != i
          for(l in 1:length(taxSeqs))
          {
            if(l != i)
            {
              # extract states(l,j) from the alignment
              states_l <- unique(as.character(DNAbin[taxSeqs[[l]],inform[j]]))
              # for j to be of type4, taxon l must be fixed for one
              # state that is different from "n" and different from the state of taxon i
              # (and different from "-" if gaps are not valid states)
              if(length(states_l)==1 && states_l[1]!="n" && states_i[1]!=states_l[1] && isTRUE(gapValid | states_l[1]!="-"))
              {
                newpos[1] <- inform[j]
                newpos[2] <- "type4"
                newpos[3] <- states_i[1]
                newpos[4] <- paste(newpos[4], unique(taxVector)[l], sep=" ")
              }
            }
          }
          # if j is identified as type4 and the type is wanted -> store position j
          if(newpos[2] == "type4")
          {
            newpos[4] <- substring(newpos[4], 2)
            newtax <- rbind(newtax, newpos)
          }
        }
      }
    }

    # make the output readable
    colnames(newtax) <- newtax[1,]
    newtax <- newtax[-1,,drop=FALSE]
    rownames(newtax) <- NULL
    taxname <- unique(taxVector)[i]
    result[[taxname]] <- newtax
  }

  # return a list of the taxa of interest and their diagnostic characters
  result
}
