


# functions to find nearest neighbors and do the grouping operations


# naive nearest neighbor finder

findOrderedNN <- function( locs, m ){
     # find the m+1 nearest neighbors to locs[j,] in locs[1:j,]
     # by convention, this includes locs[j,], which is distance 0
     n <- dim(locs)[1]
     NNarray <- matrix(NA,n,m+1)
     for(j in 1:n ){
         distvec <- c(fields::rdist(locs[1:j,,drop=FALSE],locs[j,,drop=FALSE]) )
         NNarray[j,1:min(m+1,j)] <- order(distvec)[1:min(m+1,j)]
     }
     return(NNarray)
}




# this algorithm is a modification of the kdtree
# algorithm in the FNN package, adapted to the setting
# where the nearest neighbors must come from previous
# in the ordering
#' @export
#' @importFrom FNN get.knnx
findOrderedNN_kdtree <- function(locs,m,mult=2,printsearch=FALSE){

    # number of locations
    n <- nrow(locs)

    # to store the nearest neighbor indices
    NNarray <- matrix(NA,n,m+1)
    # to the first mult*m+1 by brutce force
    NNarray[1:(mult*m+1),] <- findOrderedNN(locs[1:(mult*m+1),],m)

    query_inds <- (mult*m+2):n
    data_inds <- 1:n

    msearch <- m

    while( length(query_inds) > 0 ){

        msearch <- min( max(query_inds), 2*msearch )
        data_inds <- 1:max(query_inds)
        NN <- FNN::get.knnx( locs[data_inds,,drop=FALSE], locs[query_inds,,drop=FALSE], msearch )$nn.index
        less_than_k <- t(sapply( 1:nrow(NN), function(k) NN[k,] <= query_inds[k]  ))
        sum_less_than_k <- apply(less_than_k,1,sum)
        ind_less_than_k <- which(sum_less_than_k >= m+1)
        NN_less_than_k <- NN[ind_less_than_k,]

        NN_m <- t(sapply(ind_less_than_k,function(k) NN[k,][less_than_k[k,]][1:(m+1)] ))

        NNarray[ query_inds[ind_less_than_k], ] <- NN_m

        query_inds <- query_inds[-ind_less_than_k]
        if( printsearch ) print(length(query_inds))

    }

    return(NNarray)
}
