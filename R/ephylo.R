#' Preorder tree traversal
#'
#' Visit the nodes of a phylogeny in preorder traversal sequence.
#'
#' @param phy An object inheriting from class \code{ephylo}.
#' @param root Index of the node that will be the root of the traversal.
#' @return An integer vector of node indices in preorder traversal
#' sequence from the given root.
ephylo_preorder = function(phy, root)
{
    stopifnot(inherits(phy, "ephylo"))
    if (missing(root))
    {
        root = phy$root
    }
    else
    {
        storage.mode(root) = "integer"
        stopifnot(root > phy$num.tips && root <= phy$num.nodes)
    }
    .Call(C_ephylo_preorder, phy$num.nodes, root, phy$right.child, 
        phy$left.sib)
}

#' Return the tip indices descended from a given node
#'
#' @param phy An object inheriting from class \code{ephylo}.
#' @param root Index of the node that will be the root of the traversal.
#' @return An integer vector of tip indices descended from the given root.
ephylo_tips = function(phy, root)
{
    stopifnot(inherits(phy, "ephylo"))
    if (missing(root))
    {
        root = phy$root
    }
    else
    {
        storage.mode(root) = "integer"
        stopifnot(root > phy$num.tips && root <= phy$num.nodes)
    }
    .Call(C_ephylo_preorder_tips, phy$num.tips, phy$num.nodes, root, 
        phy$right.child, phy$left.sib)
}

#' Postorder tree traversal
#'
#' Visit the nodes of a phylogeny in postorder traversal sequence.
#'
#' @param phy An object inheriting from class \code{ephylo}.
#' @param root Index of the node that will be the root of the traversal.
#' @return An integer vector of node indices in postorder traversal
#' sequence from the given root.
ephylo_postorder = function(phy, root)
{
    stopifnot(inherits(phy, "ephylo"))
    if (missing(root))
    {
        root = phy$root
    }
    else
    {
        storage.mode(root) = "integer"
        stopifnot(root > phy$num.tips && root <= phy$num.nodes)
    }
    .Call(C_ephylo_postorder, phy$num.nodes, root, phy$parent, 
        phy$right.child, phy$left.sib)
}

#' Sample random clades from a phylogeny
#'
#' @param phy An object inheriting from class \code{ephylo}.
#' @param min.size The minimum clade size to sample.
#' @param max.size The maximum clade size to sample.
#' @param replace Should sampling be conducted with replacement?
#' @param include.root Should the root clade be a candidate for sampling?
#' @return A function to draw repeated random samples. It takes a single
#' parameter \code{n} specifying the number of random samples to take.
#' @details Sampling is conducted so that every size of clade (as measured by
#' the number of tips it contains) has an equal probability of being drawn.
#' This does not mean that you will see a uniform distribution of clade sizes
#' when drawing repeated random samples because there are still many more
#' small clades than large clades.
ephylo_rclade = function(phy, min.size, max.size, replace=FALSE,
    include.root=FALSE)
{
    stopifnot(inherits(phy, "ephylo"))
    num_tips = phy$num.tips
    num_nodes = phy$num.nodes
    ancestors = phy$root:num_nodes
    if (!include.root)
        ancestors = ancestors[-1L]
    ndesc = sapply(ancestors, function(n) length(ephylo_tips(phy, n)))
    keep = which(ndesc >= min.size & ndesc <= max.size)
    ndesc = ndesc[keep]
    ancestors = ancestors[keep]
    w = 1 / table(ndesc)
    sample_weights = unname(w[match(ndesc, as.integer(names(w)))])
    function(n) {
        stopifnot(n >= 1)
        subtrees = sample(
            ancestors, size=n, replace=replace, prob=sample_weights)
        return (subtrees)
    }
}

#' Extended \code{phylo} object
#'
#' Turn a \code{phylo} object into a quintuply linked tree.
#'
#' @param phy An object inheriting from class \code{phylo}.
#' @return An S3 object with class attribute \code{c("ephylo","phylo")}. This
#' has the usual components and node indexing schema of a \code{phylo} object
#'  and adds the following additional components:
#' \describe{
#'   \item{parent}{An integer vector with the index of every node's parent.
#'     \code{parent[u]} returns the index of node \code{u}'s immediate ancestor
#'     or \code{0} if \code{u} is the root node.}
#'   \item{left.child}{An integer vector with the index of every node's left
#'     child. \code{left.child[u]} returns the index of node \code{u}'s left
#'     child or \code{0} if \code{u} has no left child.}
#'   \item{right.child}{An integer vector with the index of every node's right
#'     child. \code{right.child[u]} returns the index of node \code{u}'s right
#'     child or \code{0} if \code{u} has no right child.}
#'   \item{left.sib}{An integer vector with the index of every node's left
#'     sibling. \code{left.sib[u]} returns the index of node \code{u}'s left
#'     sibling or \code{0} if \code{u} has no left sibling.}
#'   \item{right.sib}{An integer vector with the index of every node's right
#'     sibling. \code{right.sib[u]} returns the index of node \code{u}'s right
#'     sibling or \code{0} if \code{u} has no right sibling.}
#' }
ephylo = function(phy)
{
    stopifnot(inherits(phy, "phylo"))
    if (inherits(phy, "ephylo")) return (phy)
    edge_parent = phy$edge[,1] 
    edge_child = phy$edge[,2]
    edge_length = phy$edge.length
    num_edges = nrow(phy$edge)
    num_nodes = ape::Ntip(phy) + ape::Nnode(phy)
    tables = .Call(C_ephylo, num_nodes, num_edges, edge_parent, edge_child,
        edge_length)
    ephy = structure(
        c(unclass(phy),
            list(
              num.nodes=num_nodes
            , num.tips=ape::Ntip(phy)
            , root=ape::Ntip(phy)+1L
            , parent=tables[[1]]
            , left.child=tables[[2]]
            , right.child=tables[[3]]
            , left.sib=tables[[4]]
            , right.sib=tables[[5]]
            , brlen=tables[[6]]
            , time=ape::node.depth.edgelength(phy)
            , label=c(phy$tip.label, phy$node.label))
        )
    , class=c("ephylo", "phylo"), order=attr(phy, "order"))
    ephy$preorder.nodes = ephylo_preorder(ephy)
    ephy$postorder.nodes = ephylo_postorder(ephy)
    return (ephy)
}

#' Conversion among tree objects
#'
#' \code{as.ephylo} is a generic function that converts \code{phylo} and
#' \code{multiPhylo} objects to \code{ephylo} objects.
#'
#' @param x A \code{phylo} or \code{multiPhylo} object.
#' @param ... Additional arguments.
#' @rdname as.ephylo
as.ephylo = function(x, ...)
{
    if (inherits(x, "ephylo")) return (x)
    UseMethod("as.ephylo")
}

#' @rdname as.ephylo
as.ephylo.phylo = function(x, ...)
{
    return (ephylo(x))
}

#' @rdname as.ephylo
as.ephylo.multiPhylo = function(x, ...)
{
    for (i in 1:length(x))
        x[[i]] = as.ephylo(x[[i]])
    return (x)
}

#' S3 method for generic ape::as.phylo
#' @param x An \code{ephylo} object.
#' @param ... Additional arguments.
as.phylo.ephylo = function(x, ...)
{
    phylo_names = c("edge", "edge.length", "tip.label", "Nnode", "node.label", 
        "root.edge")
    structure(
        unclass(x)[match(phylo_names, names(x), nomatch=0L)],
        class = "phylo",
        order = attr(x, "order")
    )
}

#' S3 method for generic ape::keep.tip
#' @param phy An \code{ephylo} object.
#' @param tip Integer vector with indices of tips to keep.
#' @param ... Additional arguments.
keep.tip.ephylo = function(phy, tip, ...)
{
    as.ephylo(ape::keep.tip(ape::as.phylo(phy), tip))
}

#' S3 method for generic ape::drop.tip
#' @param phy An \code{ephylo} object.
#' @param tip Integer vector with indices of tips to drop.
#' @param ... Additional arguments.
drop.tip.ephylo = function(phy, tip, ...)
{
    as.ephylo(ape::drop.tip(ape::as.phylo(phy), tip, ...))
}

#' S3 method for generic stats::reorder that hooks into ape::reorder.phylo
#' @param x An \code{ephylo} object.
#' @param ... Additional arguments.
reorder.ephylo = function(x, ...)
{
    as.ephylo(stats::reorder(ape::as.phylo(x), ...))
}
