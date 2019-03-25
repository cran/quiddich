#' Identification of diagnostic molecular characters in a nucleotide alignment
#' that also cause diagnostic characters in the corresponding amino acid
#' alignment
#'
#' This function is a tool for an automated identification of such diagnostic
#' molecular characters in a nucleotide alignment that also cause diagnostic
#' characters in the corresponding amino acid alignment. For each taxon given
#' in \code{taxOfInt}, it identifies the mentioned characters and returns their
#' positions and types in the nucleotide alignment and their positions and
#' types in the corresponding amino acid alignment.
#'
#' @param DNAbin An object (the nucleotide alignment) of class 'DNAbin'.
#' @param taxVector The taxon vector. Default assumes that each row in the
#' alignment belongs to a different taxon.
#' @param codonstart An integer defining where to start the translation of the
#' nucleotide alignment into the corresp. amino acid alignment. Default is 1.
#' @param taxOfInt A vector containing the taxa for which the above mentioned
#' characters shall be extracted. Default is "all".
#' @param types A vector containing the types of diagnostic characters
#' (in the nucleotide alignment) that shall be considered. The types
#' can be "all" or any combination of "type1", "type2", "type3" and "type4".
#' Default is "type1", "type 2" and "type3".
#'
#' @return \code{changesAA} returns a list, where each entry belongs to one
#' taxon of interest. Each taxon of interest has a matrix assigned to it, in
#' which the first and second column contain the positions and types of the
#' identified diagnostic characters in the nucleotide alignment and the third
#' and fourth column contain the positions and types of the diagnostic
#' characters in the corresponding amino acid alignment.
#'
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
#' \code{changesAA} returns for each taxon of interest the following elements:
#' \item{posDiNu}{The positions of the identified diagnostic characters in the
#' nucleotide alignment.}
#' \item{typeDiNu}{The types of the identified diagnostic characters in the
#' nucleotide alignment.}
#' \item{posDiAA}{The positions of the corresponding diagnostic characters in
#' the amino acid alignment.}
#' \item{typeDiAA}{The types of the corresponding diagnostic characters in the
#' amino acid alignment.}
#'
#' @author A. Luise Kuehn <luise.kuehn@@uni-greifswald.de>
#'
#' @references Kuehn, A.L., Haase, M. (2019) QUIDDICH: QUick IDentification of
#' DIagnostic CHaracters.
#'
#' @examples
#' #using a dataset from spider
#' #install.packages("spider")
#' library(spider)
#' data("anoteropsis")
#' anoTax <- sapply(strsplit(dimnames(anoteropsis)[[1]], split="_"),
#'  function(x) paste(x[1], x[2], sep="_"))
#' changesAA(anoteropsis, anoTax, codonstart=2, taxOfInt="all")
#' changesAA(anoteropsis, anoTax, codonstart=2, taxOfInt="Artoria_flavimanus",
#'  types="type4")
#'
#' @importFrom ape trans
#'
#' @export
changesAA <- function(DNAbin, taxVector=dimnames(DNAbin)[[1]], codonstart=1, taxOfInt="all", types=c("type1","type2","type3"))
{
  # specify taxa of interest
  if(taxOfInt[1]=="all")
  {
    taxOfInt <- taxVector
  }
  # specify types of interest (in nucleotide alignment)
  if(types[1]=="all")
  {
    types <- c("type1","type2","type3","type4")
  }

  # translate nucleotide alignment into amino acid alignment using standard genetic code
  # (start of the translation is codonstart)
  AAbin <- ape::trans(DNAbin,code=1,codonstart=codonstart)
  # extract diag. char. from nucleotide and amino acid alignment
  diNu <- diagCharNA(DNAbin,taxVector,taxOfInt=taxOfInt,types=types)
  diAA <- diagCharAA(AAbin,taxVector,taxOfInt=taxOfInt,types="all")

  # result stores for each taxon of interest the diag. char. of the nucleotide alignment that
  # also cause diag. char. in the amino acid alignment
  result <- list()

  # consider all taxa of interest successively
  for(i in 1:length(unique(taxOfInt)))
  {
    # get the name of the current taxon of interest
    curTax <- unique(taxOfInt)[i]
    # create new list element for curTax
    newtax <- matrix(c("posDiNu","typeDiNu","posDiAA","typeDiAA"), ncol=4)
    # consider all identified diagn. char. from nucleotide alignment successively
    if(length(diNu[[curTax]])!=0)
    {
      for(j in 1:length(diNu[[curTax]][,1]))
      {
        curNuc <- as.numeric(diNu[[curTax]][j,1])
        # check if current diag. char. of nucleotide alignment is also diag. char. in amino acid alignment
        curAA <- which(diAA[[curTax]][,1]==as.character(ceiling((curNuc-codonstart+1)/3)))
        if(length(curAA)!=0)
        {
          # if current diag. char. of nucleotide alignment causes diag. char. in amino acid alignment, store:
          # position in nucleotide alignment (posDiNu)
          # type of char. in nucleotide alignment (typeDiNu)
          # position in amino acid alignment (posDiAS)
          # type of char. in amino acid alignment (typeDiAS)
          newtax <- rbind(newtax,c(diNu[[curTax]][j,1], diNu[[curTax]][j,2], diAA[[curTax]][curAA,1], diAA[[curTax]][curAA,2]))
        }
      }
    }
    # make the output readable
    colnames(newtax) <- newtax[1,]
    newtax <- newtax[-1,,drop=FALSE]
    rownames(newtax) <- NULL
    result[[curTax]] <- newtax
  }
  # return the result
  result
}
