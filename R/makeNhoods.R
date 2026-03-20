#' Define neighbourhoods on a graph (fast)
#'
#' This function randomly samples vertices on a graph to define neighbourhoods.
#' These are then refined by either computing the median profile for the neighbourhood
#' in reduced dimensional space and selecting the nearest vertex to this
#' position (refinement_scheme = "reduced_dim"), or by computing the vertex with the highest number of
#' triangles within the neighborhood (refinement_scheme = "graph").
#' Thus, multiple neighbourhoods may be collapsed down together to
#' prevent over-sampling the graph space.
#' @param x A \code{\linkS4class{Milo}} object with a non-empty \code{graph}
#' slot. Alternatively an \code{igraph} object on which neighbourhoods will
#' be defined.
#' @param prop A double scalar that defines what proportion of graph vertices
#' to randomly sample. Must be 0 < prop < 1.
#' @param k An integer scalar - the same k used to construct the input graph.
#' @param d The number of dimensions to use if the input is a matrix of cells
#' X reduced dimensions.
#' @param reduced_dims If x is an \code{\linkS4class{Milo}} object, a character indicating the name of the \code{reducedDim} slot in the
#' \code{\linkS4class{Milo}} object to use as (default: 'PCA'). If x is an \code{igraph} object, a
#' matrix of vertices X reduced dimensions with \code{rownames()} set to correspond to the cellIDs.
#' @param refined A logical scalar that determines the sampling behavior, default=TRUE implements a refined sampling scheme,
#' specified by the refinement_scheme argument.
#' @param refinement_scheme A character scalar that defines the sampling scheme, either "reduced_dim" or "graph".
#' Default is "reduced_dim".
#' @param BPPARAM A \code{BiocParallelParam} object for parallelizing KNN searches
#' in refined sampling. Defaults to \code{SerialParam()}.
#' @param BNPARAM A \code{BiocNeighborParam} object specifying the KNN algorithm.
#' Defaults to \code{KmknnParam()} (exact). Use \code{HnswParam()} or
#' \code{AnnoyParam()} for approximate methods on large datasets.
#'
#' @details
#' This function randomly samples graph vertices, then refines them to collapse
#' down the number of neighbourhoods to be tested. The refinement behaviour can
#' be turned off by setting \code{refine=FALSE}, however, we do not recommend
#' this as neighbourhoods will contain a lot of redundancy and lead to an
#' unnecessarily larger multiple-testing burden.
#'
#' @return A \code{\linkS4class{Milo}} object containing a list of vertices and
#' the indices of vertices that constitute the neighbourhoods in the
#' nhoods slot. If the input is a \code{igraph} object then the output
#' is a matrix containing a list of vertices and the indices of vertices that constitute the
#' neighbourhoods.
#'
#' @author
#' Emma Dann, Mike Morgan
#'
#' @examples
#'
#' require(igraph)
#' m <- matrix(rnorm(100000), ncol=100)
#' milo <- buildGraph(m, d=10)
#'
#' milo <- makeNhoods(milo, prop=0.1)
#' milo
#'
#' @export
#' @rdname makeNhoods
#' @importFrom BiocNeighbors findKNN queryKNN
#' @importFrom igraph neighbors neighborhood as_ids V
#' @importFrom Matrix sparseMatrix
#' @importFrom stats setNames
makeNhoods <- function(x, prop=0.1, k=21, d=30, refined=TRUE, reduced_dims="PCA",
                       refinement_scheme = "reduced_dim",
                       BPPARAM = BiocParallel::SerialParam(),
                       BNPARAM = BiocNeighbors::KmknnParam()) {
    if(is(x, "Milo")){
        message("Checking valid object")
        # check that a graph has been built
        if(!.valid_graph(miloR::graph(x))){
            stop("Not a valid Milo object - graph is missing. Please run buildGraph() first.")
        }
        x.graph <- miloR::graph(x)

        if(isTRUE(refined) & refinement_scheme == "reduced_dim"){
            X_reduced_dims  <- reducedDim(x, reduced_dims)
            if (d > ncol(X_reduced_dims)) {
                warning("Specified d is higher than the total number of dimensions in reducedDim(x, reduced_dims).
                        Falling back to using",ncol(X_reduced_dims),"dimensions\n")
                d <- ncol(X_reduced_dims)
            }
            X_reduced_dims  <- X_reduced_dims[,seq_len(d)]
            mat_cols <- ncol(x)
            match.ids <- all(rownames(X_reduced_dims) == colnames(x))
            if(!match.ids){
                stop("Rownames of reduced dimensions do not match cell IDs")
            }
        }

    } else if(is(x, "igraph")){

        if(isTRUE(refined) & refinement_scheme == "reduced_dim" & !is.matrix(reduced_dims)) {
            stop("No reduced dimensions matrix provided - required for refined sampling with refinement_scheme = reduced_dim.")
        }

        x.graph <- x

        if(isTRUE(refined) & refinement_scheme == "reduced_dim"){
            X_reduced_dims  <- reduced_dims
            mat_cols <- nrow(X_reduced_dims)
            if(is.null(rownames(X_reduced_dims))){
                stop("Reduced dim rownames are missing - required to assign cell IDs to neighbourhoods")
            }
        }

        if(isTRUE(refined) & refinement_scheme == "graph" & is.matrix(reduced_dims)){
            warning("Ignoring reduced dimensions matrix because refinement_scheme = graph was selected.")
        }

    } else{
        stop("Data format: ", class(x), " not recognised. Should be Milo or igraph.")
    }

    random_vertices <- .sample_vertices(x.graph, prop, return.vertices = TRUE)

    if (isFALSE(refined)) {
        sampled_vertices <- random_vertices
    } else if (isTRUE(refined)) {
        if(refinement_scheme == "reduced_dim"){
            sampled_vertices <- .refined_sampling(random_vertices, X_reduced_dims, k,
                                                    BPPARAM = BPPARAM, BNPARAM = BNPARAM)
        } else if (refinement_scheme == "graph") {
            sampled_vertices <- .graph_refined_sampling(random_vertices, x.graph)
        } else {
            stop("When refined == TRUE, refinement_scheme must be one of \"reduced_dim\" or \"graph\".")
        }
    } else {
        stop("refined must be TRUE or FALSE")
    }

    sampled_vertices <- unique(sampled_vertices)

    # Determine row count and names for the nhood matrix
    v.class <- V(x.graph)$name
    if(is(x, "Milo")){
        n_rows <- ncol(x)
        rnames <- colnames(x)
    } else if(is(x, "igraph")){
        n_rows <- length(V(x))
        if(is.null(v.class) & refinement_scheme == "reduced_dim"){
            rnames <- rownames(X_reduced_dims)
        } else if(!is.null(v.class)){
            rnames <- V(x.graph)$name
        } else {
            rnames <- NULL
        }
    }

    # Build nhood matrix from COO triplets (vectorized neighborhood call)
    nhood_lists <- neighborhood(x.graph, order = 1, nodes = sampled_vertices)
    col_indices <- rep(seq_len(length(sampled_vertices)), lengths(nhood_lists))
    row_indices <- unlist(nhood_lists)

    nh_mat <- sparseMatrix(
        i = row_indices, j = col_indices,
        x = rep(1, length(row_indices)),
        dims = c(n_rows, length(sampled_vertices)),
        dimnames = list(rnames, as.character(sampled_vertices))
    )
    if(is(x, "Milo")){
        nhoodIndex(x) <- as(sampled_vertices, "list")
        nhoods(x) <- nh_mat
        return(x)
    } else {
        return(nh_mat)
    }
}


#' @importFrom BiocNeighbors findKNN queryKNN
#' @importFrom matrixStats colMedians
.refined_sampling <- function(random_vertices, X_reduced_dims, k,
                              BPPARAM = BiocParallel::SerialParam(),
                              BNPARAM = BiocNeighbors::KmknnParam()){
    message("Running refined sampling with reduced_dim")
    vertex.knn <-
        findKNN(
            X = X_reduced_dims,
            k = k,
            subset = as.vector(random_vertices),
            get.index = TRUE,
            get.distance = FALSE,
            BPPARAM = BPPARAM,
            BNPARAM = BNPARAM
        )

    nh_reduced_dims <- t(apply(vertex.knn$index, 1, function(x) colMedians(X_reduced_dims[x,])))

    if(is.null(rownames(X_reduced_dims))){
        warning("Rownames not set on reducedDims - setting to row indices")
        rownames(X_reduced_dims) <- as.character(seq_len(nrow(X_reduced_dims)))
    }

    colnames(nh_reduced_dims) <- colnames(X_reduced_dims)

    ## Search nearest cell to each median profile using queryKNN
    nn_result <- queryKNN(
        X = X_reduced_dims,
        query = nh_reduced_dims,
        k = 1,
        BPPARAM = BPPARAM,
        BNPARAM = BNPARAM
    )
    sampled_vertices <- as.integer(nn_result$index[, 1])
    return(sampled_vertices)
}

#' @importFrom igraph is_igraph
.valid_graph <- function(x){
    # check for a valid graph
    if(isTRUE(is_igraph(x))){
        TRUE
    } else{
        FALSE
    }
}

#' @import igraph
.sample_vertices <- function(graph, prop, return.vertices=FALSE){
    # define a set of vertices and neihbourhood centers - extract the neihbourhoods of these cells
    random.vertices <- sample(V(graph), size=floor(prop*length(V(graph))))
    if(isTRUE(return.vertices)){
        return(random.vertices)
    } else{
        message("Finding neighbours of sampled vertices")
        vertex.list <- sapply(seq_len(length(random.vertices)), FUN=function(X) neighbors(graph, v=random.vertices[X]))
        return(list(random.vertices, vertex.list))
    }
}

#' @importFrom igraph count_triangles neighborhood set_vertex_attr induced_subgraph V
.graph_refined_sampling <- function(random_vertices, graph){
    message("Running refined sampling with graph")
    random_vertices <- as.vector(random_vertices)
    X_graph <- set_vertex_attr(graph, "name", value = 1:length(V(graph)))
    refined_vertices <- lapply(seq_along(random_vertices), function(i){
        target_vertices <- unlist(neighborhood(X_graph, order = 1, nodes = random_vertices[i])) #get neighborhood of random vertex
        target_vertices <- target_vertices[-1] #remove first entry which is the random vertex itself
        rv_induced_subgraph <- induced_subgraph(graph = X_graph, vids = target_vertices)
        triangles <- count_triangles(rv_induced_subgraph)
        max_triangles <- max(triangles)
        max_triangles_indices <- which(triangles == max_triangles)
        #note - take first max_ego_index in the next line of code
        resulting_vertices <- V(rv_induced_subgraph)[max_triangles_indices]$name[1]
        return(resulting_vertices)
    }) %>% unlist() %>% as.integer()
    return(refined_vertices)
}
