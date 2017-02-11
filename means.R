# Code for computing various types of averages.

EPSILON <- 1e-8

mean.norm <- function(vec, norm = 2) {

    if(is.infinite(norm) && norm > 0)
        return(max.or.min(vec))
        
    if(is.infinite(norm) && norm < 0)
        return(max.or.min(vec, FALSE))
    
    return(normp(vec, norm)/(length(vec)^(1/norm)))
}

#This is the straightforward geometric mean function that does not adjust
# for values of 0. 
# May produce undefined results for vectors containing negative values.
mean.geom <- function(vec) {

    p <- prod(vec)
    
    geom <- p ^(1/length(vec))
 
    geom
}

#This is the geometric mean where the vector is adjusted for values of 0.
# There are 4 methods of adjusting for 0: 
# "change": This replaces all 0's in the vector with 1's, computes the 
# geometric mean of this new vector, then subtracts 1 from the result.
# "add": This adds 1 to each element of the vector, computes the geometric mean
# of this new vector, then subtracts 1 from the result.
# "ignore": This removes 0's from the vector and computes the geometric mean
# of the zero-less vector.
# "comp" (default): This is an original, complicated way of replacing 0's with
# appropriate nearby values based on the nonzero elements of the vector.
# If an invalid method is passed in, "comp" will be used.
mean.geom.zadj <- function(vec, method = "comp") {
    
    method <- tolower(method)

    #Changes 0 to 1
    if(method == "change")  {
        for(i in 1:length(vec)) {
            if(vec[i] == 0)
                vec[i] = 1
        }

        return(mean.geom(vec) - 1)
            
    } else if(method == "add") {
          
        #Adds 1 to each value
        vec <- vec + 1
        return(mean.geom(vec) - 1)
            
    } else if(method == "ignore"){
        
        i <- 1
        
        #Removes 0 from "vec"
        while(i <= length(vec)) {
            if(vec[i] == 0) {
                vec <- vec[-1*i]
                i <- i - 1
            }

            i <- i + 1
        }
        
        return(mean.geom(vec))
    } else
        return(mean.zadj(vec, mean.geom))
}

#This computes the harmonic mean of the vector "vec". May produce undefined
# results for vectors containing negative values.
mean.harmonic <- function(vec) {
    sum <- 0

    for(i in 1:length(vec))
        sum <- sum + abs.val(1/vec[i])

    return(length(vec) / sum)
}

# This removes values of 0 from a vector and replaces the 0's with appropriate
# nearby values. The most obvious use is for computing the adjusted geometric mean of a vector
# containing 0's.
# "func" is the mean function to apply to the vector "vec". 
# "see.vec" should be TRUE if the modified (zero-adjusted) vector should be 
# part of the return value; FALSE if only the mean should be returned.
# The only time 0's are not removed is if "vec" contains only values of 0.
mean.zadj <- function(vec, func, see.vec = FALSE) {
    
    vec <- sort(vec)

    #Execute this loop until no more zeroes in "vec".
    while(TRUE) {

        # Indices in "vec" that contain 0's. 
        zind <- which(vec == 0)

        # Number of 0's.
        zeroes <- length(zind)
        
        # If there are no 0's, the loop is done.
        if(zeroes == 0)
            break

        # Lowest index in "vec" containing 0.
        lower <- min(which(vec == 0))

        # Greatest index in "vec" containing 0.
        upper <- max(which(vec == 0))

        # Subvector of "vec" containing negative values.
        lvec <- vec[vec < 0]

        # Subvector of "vec" containing positive values.
        uvec <- vec[vec > 0]

        # If the vector is ALL 0's, no replacement should occur.
        if(length(lvec) == 0 && length(uvec) == 0)
            break

        lmean <- -1*func(abs(lvec))
        umean = func(uvec)
        
        #If there is only one zero
        if(lower == upper && length(lvec) != 0 && length(uvec) != 0) {

            minU = min(uvec)
            maxL = max(lvec)

            # Trying to decide the best possible value to replace the 0 with.
            extremes <- minU + maxL
            lowMax.meanU <- maxL + umean
            lowMax.meanU.mult <- maxL + umean*minU
            highMin.meanL <- minU + lmean
            highMin.meanL.mult <- minU + lmean*maxL
            means <- lmean + umean
            all <- lmean*maxL + umean*minU

            possibles <- c(extremes, lowMax.meanU, lowMax.meanU.mult, highMin.meanL,
                           highMin.meanL.mult, means, all)
            
            # Pick the number of those numbers closest to 0.
            abs <- abs(possibles)
            index <- abs == min(abs)
            
            vec[lower] <- possibles[index]
            break
        } else {
            
            # If there is more than 1 zero.
            
            # How to handle negative numbers
            if(length(lvec) != 0) {

               # Ideal value to replace 0 with, unless this number is also 0.
               low <- lmean - vec[lower - 1]
            
               if(low < 0 && low > vec[lower - 1])
                   vec[lower] <- low
               else {

                   # If "low" is 0 or is too negative to replace the 0 with, use this
                   # formula which guarantees that the 0 will be replaced with a negative number
                   # greater than the max of the negative numbers in "vec".
                   vec[lower] <- vec[lower-1]*lmean / ((lmean-1)*(zeroes+length(lvec)))
               }
           }
        
            # How to handle positive numbers
            if(length(uvec) != 0) {

                # Ideal value to replace 0 with, unless this number is also 0.
                high <- umean - vec[upper + 1]
            
                if(high > 0 && high < vec[upper + 1])
                    vec[upper] <- high
                else {

                    # If "high" is 0 or is too positive to replace the 0 with, use this
                    # formula which guarantees that the 0 will be replaced with a positive number
                    # greater than the min of the positive numbers in "vec".
                    vec[upper] <- vec[upper+1]*umean / ((umean+1)*(zeroes+length(uvec)))
                }
            }
        }   
    }

    if(see.vec)
        list(newvector = vec, mean = func(vec))
    else
        func(vec)
}

# Returns TRUE if "num" is an actual number (i.e. not NaN, not Inf, not NA and not NULL)
is.good <- function(num) {

    for(i in 1:length(num))
        if(is.nan(num) || is.infinite(num) || is.na(num) || is.null(num))
            return(FALSE)

    return(TRUE)
}

#Removes NA, Inf, NaN and NULL from a vector "vec"
without.bad <- function(vec) {

    i <- 1
    
    while(i <= length(vec)) {
        
        if(!is.good(vec[i])) {
            
            vec <- vec[-1*i]
            i <- i - 1
        }

        i <- i + 1
    }

    vec
}

#func's 2 arguments will  be applied to
# 2 of the elements of vecs - which is a vector of numbers or vectors. func must return a single number. If vectors, pass in a list. "remove" indicates whether or not NA, NaN, Inf or NULL should be removed from a vector. The vectors need not be the same length. 
# Example use: to find the difference between every number in a vector and every other number in 
# the vector, pass in the vector as "vecs" and a 2-argument subtraction function as "func", e.g: 
# func = function(x,y) x - y
# This could also be used to compute a covariance matrix if "vecs" is a list of vectors.
cross.matrix <- function(vecs, type = "num", remove = FALSE, func) {
   type <- tolower(type)

   if(type != "vec")
       type = "num"

   if(remove && type == "num")
       vecs <- without.bad(vecs)

   else if(remove) {
       for(i in 1:length(vecs))
           vecs[[i]] <- without.bad(vecs[[i]])
   }
       
   i <- 1
   
   while(i <= length(vecs)) {

       if(length(vecs[i]) == 0) {

           vecs <- vecs[-1*i]
           i <- i - 1
       }
       
       i <- i + 1
   }

   sizes <- numeric(length(vecs))

   if(type == "vec") {
       for(i in 1:length(vecs))
           sizes[i] <- length(vecs[[i]])
   }
       
   if(type == "vec") {
    
       newlist <- vector("list", length(vecs))

       if(remove) {
           
           # If "vecs" is a list of vectors, and the vectors have different lengths,
           # the minimum of these lengths is considered so that "newlist" can be "square" (like 
           # a list version of a square matrix) without any NA values.
           newlen <- min(sizes)
           
           for(i in 1:length(newlist)) {
               for(j in 1:newlen)
                   newlist[[i]][[j]] <- vecs[[i]][[j]]
           }
           
       } else {

           newlen <- max.or.min(sizes)

           # NA may occur if vectors are not of the same length
           for(i in 1:length(newlist)) {
               for(j in 1:newlen) {
                   if(length(vecs[[i]]) <= newlen)
                      newlist[[i]][[j]] <- vecs[[i]]
                   else
                       newlist[[i]][[j]] <- NA
               }
           }
       }

       vecs <- newlist
   }

   cross <- matrix(0, length(vecs), length(vecs))

   # Fill matrix with result of function applications
   for(i in 1:nrow(cross)) {
       for(j in 1:ncol(cross)) {
           if(remove || is.good(vecs[[i]]) && is.good(vecs[[j]]))
              cross[i, j] <- func(vecs[[i]], vecs[[j]])
           else
               cross[i, j] <- NA
       }
    }
           
   cross
}

#Computes the "inner" mean of a set of means. Computing AGM is one obvious use, but this is more general
# because the "agm" function already exists.
# Arguments: 
# vec: vector to find mean of
# funcs: vector of functions that will operate on "vec"
# arg1: List of first arguments for funcs. (i.e., if "funcs" contains 2 functions and the first function
#       takes only the vector as an argument, but the second function takes both the vector and another 
#       parameter, arg1 should be the vector [NA arg] where arg is the argument for the second function
# arg2: List of second arguments for funcs. Put NA in this vector for functions that do not have a second
#      argument.
# tol: tolerance level; used to determine convergence.
# maxitt: maximum number of iterations
# To use this function to compute AGM (arithmetic-geometric mean), funcs should contain 
# the "mean" and "mean.geom" or "mean.geom.zadj" functions. 
# If funcs = c(mean, mean.geom.zadj) and a non-default method for adjusting for zero is chosen, arg1 should be c(NA, method).
mean.inner <- function(vec, funcs, arg1 = rep(NA, length(vec)), arg2 = rep(NA, length(vec)),
                       tol = EPSILON, maxitt = 40) {
        
    j <- 1
    start <- vec
    nextvec <- numeric(length(funcs))
    
    while(TRUE) {

        # Apply function to vector "start" with desired arguments.
        for(i in 1:length(funcs)) {

            if(is.na(arg1[i]) && is.na(arg2[i]))
                nextvec[i] <- funcs[[i]](start)

            else if(is.na(arg1[i]) && !is.na(arg2[i]))
                nextvec[i] <- funcs[[i]](start, arg2[i])

            else if(!is.na(arg1[i]) && is.na(arg2[i]))
                nextvec[i] <- funcs[[i]](start, arg1[i])

            else 
                nextvec[i] <- funcs[[i]](start, arg1[i], arg2[i])
        }
        
        # All possible differences
        cross.diffs <- cross.matrix(vecs = nextvec, remove = TRUE, func = function(x,y) x - y)

        # Convergence achieved
        if(max(cross.diffs) <= tol)
            break

        # There could be values of Inf. Should not happen with mean and geometric mean but 
        # for functions in general this could happen.
        if(length(cross.diffs) > length(without.bad(cross.diffs)) || j > maxitt)
            stop("Possible divergence")
        
        start <- nextvec
        j <- j + 1
    }

    return(nextvec[1])
}
