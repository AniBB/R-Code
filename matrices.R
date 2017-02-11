EPSILON <- 1e-8

# Replaces values in a matrix or vector with 0 if the values are very close to 0 (within "tol" of 0).
# Helpful when dealing with roundoff errors from calculations.
round.to.zero <- function(x, tol = EPSILON) {

    if(class(x) == "matrix")
        y <- as.numeric(x)
    else
        y <- x
    
    for(i in 1:length(y))
        if(abs(y[i]) < tol)
            y[i] <- 0

    if(class(x) == "matrix")
        return(matrix(y, nrow(x), ncol(x)))

    y
}

# Returns the identity matrix of size n x n if "x" is a square matrix of dimensions "n". 
# If "x" is not square, this function creates a matrix of the same dimensions of "x" with all 0's except 1's on the 
# "diagonal" starting from the upper left element and going diagonally down as much as possible.
# Example: "x" is of size 3 x 4. The return value will be the following matrix: 
# [1 0 0 0 
#  0 1 0 0 
#  0 0 1 0]
# (diag(n) is a built-in function that creates an identity matrix of dimensions n x n , which is why this function is   named "diag.nonsq" to stress the fact that this does something similar for nonsquare matrices as well as square ones. 
diag.nonsq <- function(x) {

    if(nrow(x) == ncol(x))
        return(diag(nrow(x)))

    y <- matrix(0, nrow(x), ncol(x))

    for(i in 1:min.two(nrow(x), ncol(x)))
        y[i, i] <- 1

    return(y)
} 
#ARGUMENTS:
# x: a matrix 
# decomp: The decomposition to be performed: "none", "LU" or "LDV" are the options.
# upTri: TRUE, if "x" should be made upper triangular (using row operations); FALSE, if "x" should be made lower triangular (using row operations)
# special: TRUE if "x" should be made special upper or special lower triangular. (1's on the diagonal)
# diag: TRUE if the diagonal matrix produced from the special upper or lower triangular operations should be part of the output; FALSE otherwise.
# augment: TRUE if the matrix "x" should be augmented with the identity matrix (or something close if the matrix is not square); FALSE otherwise
# is.aug: TRUE to indicate the matrix "x" is augmented already; FALSE otherwise
# solvevec: A vector to be manipulated with the row operations performed on "x" (generally the "right hand side" of a system of linear equations); NULL if no such operations are desired.
# has.vec: TRUE if "x" already has a vector attached (augmented) to it; FALSE otherwise.

# For example, if "x" is the matrix
# [4  1  7 
#  5 -3  6 
# 19  0  8]
# and we want to attach the vector [0 1 18], this vector would be "solvevec", has.vec would be false, and is.aug would be FALSE.

#If decomp is "lu" or "ldv", other arguments are ignored. If decomp is "lu", or decomp is "none" and special is FALSE, diag is ignored.

mat.tri <- function (x, decomp = "none", upTri = TRUE, special = FALSE, diag = FALSE, augment = FALSE, is.aug = FALSE, 
                     solvevec = NULL, has.vec = FALSE) {
    #origX will be "x" before any modifications occur
    origX <- x
    decomp <- toupper(decomp)

    if(nrow(x) == 0) 
        stop("Matrix must contain elements")

    if(!is.null(solvevec) && length(solvevec) != nrow(x))
        stop("Right hand side must be same length as matrix columns")

    # "origcols" is the number of columns that the "original" matrix has. The original matrix is not "x"; it is what
    # "x" would be if it is not augmented AND has no vector attached.
    if(is.aug && has.vec)
        origcols <- (ncol(x) - 1) / 2
    else if(is.aug)
        origcols <- ncol(x) / 2
    else if(has.vec)
        origcols <- ncol(x) - 1
    else
        origcols <- ncol(x)
    
    if(!is.null(solvevec))
        x <- cbind(x, solvevec)

    if(augment)
        x <- cbind(x, diag.nonsq(x))

    # The following 3 "if" statements are mainly to correct for user errors in the inputs to this function.

    # There is no need for augmenting or vectors if doing an LU or LDV decomposition.
    # We need upper triangular for these decompositions; the lower triangular portion of the decomposition will
    # be computed afterward.
    if(decomp == "LU" || decomp == "LDV") {
        upTri <- TRUE
        augment <- FALSE
        solvevec <- NULL
    }

    # For LDV, the diagonal matrix is needed (the diagonal matrix is the "D" in "LDV").
    # Special upper triangular is also needed (this is the "V" in "LDV")
    if(decomp == "LDV") {
        diag <- TRUE
        special <- TRUE
    }

    # If special upper or lower is not desired, there will be no diagonal matrix 
    if(decomp == "LU" || (decomp != "LU" && decomp != "LDV" && special == FALSE))
       diag <- FALSE

    # This will be the permutation matrix in case an "LU" decomposition needs to be a permuted "LU" decomposition.
    p <- diag(nrow(x))

    # "tri" is the out-parameter that will contain the upper or lower triangular matrix, and "diagmat" will contain the
    # diagonal matrix if "diag" is TRUE. The matrix "x" is transposed because R stores matrices in column-major order 
    # whereas C stores arrays in row-major order. 
    out <- .C("trinonsq",
              rows = as.integer(nrow(x)),
              cols = as.integer(ncol(x)),
              origcols = as.integer(origcols),
              decp = as.integer(decomp == "LU" || decomp == "LDV"),
              upper = as.integer(upTri),
              special = as.integer(special),
              diagonal = as.integer(diag),
              tri = as.double(t(x)),
              pmat = as.double(p),
              diagmat = double(nrow(x)*nrow(x)),
             PACKAGE = "matrices")

    # TRUE needed as an argument to the matrix function because of the row vs. column-major order issue:
    # TRUE ensures that the elements in the vector (first argument) will be read in row-major order.
    ret <- matrix(out$tri, nrow(x), ncol(x), TRUE)
    ret <- round.to.zero(ret)
    
    diagonal <- matrix(out$diagmat, nrow(x), nrow(x), TRUE)
    
    if(decomp == "LU" || decomp == "LDV") {

        # permutation matrix
        perm <- matrix(out$pmat, nrow(x), nrow(x), TRUE)
        upper <- ret
        permX <- perm %*% origX

        # For "LU" decomposition, A = LU or more generally, PA = LU. This is why "upTri" is TRUE before the C
        # function call. This way, the L (lower triangular matrix) is easy to compute. It is just P*A*inverse of U.
        if(decomp == "LU")
            lower <- permX %*% solve(upper)

        else {

            # If the decomposition is "LDV", "upper" is special upper triangular: PA = LDV. But DV = U (the same U
            # in "LU" decomposition), so "realUpper" is this U (DV). 
            realUpper <- diagonal %*% upper
            lower <- permX %*% mat.inv(realUpper)
        }

        lower <- round.to.zero(lower)
        
        if(decomp == "LU") {
            ret <- list(permutation = perm, lower = lower, upper = upper)
            return(ret)
        }

        ret <- list(permutation = perm, lower = lower, diagonal = diagonal, upper = upper)
        return(ret)
    }

    if(diag == TRUE) {
        
        if(upTri)
            ret <- list(diagonal = round.to.zero(diagonal), upper = ret)
        else
            ret <- list(diagonal = round.to.zero(diagonal), lower = ret)
    }

    ret
}
