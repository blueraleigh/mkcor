#' Read a tree from a Newick file or string
#'
#' @param file The path to the Newick file.
#' @param text The Newick string. Ignored if \code{file} is specified.
#' @param num_tips The number of tips in the phylogeny (optional). 
#' Specifying the number of tips in advance will speed the parsing
#' of large trees by eliminating the need for memory reallocations.
#' @return An S3 object with class attribute \code{c("ephylo", "phylo")}.
#' @seealso \code{\link{ephylo}}
read_newick = function(file, text, num_tips)
{
    bufsize = if (missing(num_tips)) 8192L else as.integer(2L*num_tips)
    if (!missing(file)) {
        filename = normalizePath(file, mustWork=TRUE)
    } else {
        stopifnot(!missing(text))
        filename = tempfile()
        cat(text, file=filename)
        on.exit(unlink(filename))
    }
    tables = .Call(C_ephylo_read_newick, filename, bufsize)
    ephy = structure(list(
          edge=tables[[8]]
        , edge.length=tables[[6]][tables[[8]][,2L]]
        , Nnode=tables[[10]]
        , node.label=tables[[7]][-(1:tables[[9]])]
        , tip.label=tables[[7]][1:tables[[9]]]
        , num.nodes=tables[[9]] + tables[[10]]
        , num.tips=tables[[9]]
        , root=tables[[9]]+1L
        , parent=tables[[1]]
        , left.child=tables[[2]]
        , right.child=tables[[3]]
        , left.sib=tables[[4]]
        , right.sib=tables[[5]]
        , brlen=tables[[6]]
        , time=tables[[11]]
        , label=tables[[7]])
    , class=c("ephylo", "phylo"), order="cladewise")
    ephy$preorder.nodes = ephylo_preorder(ephy)
    ephy$postorder.nodes = ephylo_postorder(ephy)
    if (all(is.na(ephy$node.label))) ephy$node.label = NULL
    return (ephy)
}
