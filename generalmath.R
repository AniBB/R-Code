MAXBASE <- 36
MINBASE <- 2
EPSILON <- 1e-8

# Finds the p-norm of a vector "vec" where "norm" = p
norm.p <- function(vec, norm) {

    if(is.infinite(norm))
        return(max(sapply(vec, abs)))

    powers <- sapply(vec, function(x) x^norm)
    sum <- sum(powers)

    return(sum ^ (1/norm))
}

# Converts a nonnegative integer "num" from base "oldbase" into base
# base "newbase". Only works for bases between 2 and 36 because larger 
# bases may require the use of other characters, which may be confusing to 
# users. (10 digits + 26 letters = 36 possible characters; having lower- and 
# uppercase letters represent different digits may be confusing).
basechange.int <-  
function(num, oldbase = 10, newbase)
{
    if(num < 0)
        stop("Value must be nonnegative")
    
    if(oldbase < MINBASE || oldbase > MAXBASE ||
       newbase < MINBASE || newbase > MAXBASE) stop("Invalid base")

    if(!poss.base(num, oldbase))stop("Base is invalid for this
        number.")
    
    if(newbase == oldbase)
        return(as.character(num))

    if(num == 0)
        return("0")
    
    if(oldbase == 10) {

        currnum <- as.numeric(num)
        numvec <- character(0)
        
        while(currnum > 0)
        {
            remainder <- allchars[currnum %% newbase + 1]

            numvec <- c(remainder, numvec)
        
            currnum <- my.floor(currnum / newbase)
        }

        ans <- paste(numvec, sep = "", collapse = "")
        return(ans)
    }

    else {

        # Convert from other base ("oldbase") to base 10.
        currnum <- toupper(as.character(num))
        totdigits <- nchar(currnum)
        decimal <- 0
        
        for(i in 1:totdigits) {

            dig <- which(allchars == substr(currnum, i, i)) - 1
            power <- totdigits - i
            decimal <- decimal + dig*power.simp(oldbase, power)
        }

        # Convert from base 10 to desired "newbase"
        return(basechange.int(decimal, 10, newbase))
   }
}

# Same as basechange.int except works for floating point values as well.
# "acc" is the maximum number of digits to the right of the decimal point 
# to be included in the final answer. (Needed to prevent infinite loops in 
# the case of a repeating decimal.)
basechange.float <- function(num, old.base = 10, new.base, acc = 12) {

    if(old.base < MINBASE || old.base > MAXBASE ||
       new.base < MINBASE || new.base > MAXBASE) stop("Invalid base")
    
    if(!poss.base(num, old.base))stop("Base is invalid for this
        number.")
    
    if(new.base == old.base)
        return(as.character(num))

    str <- toupper(as.character(num))
    
    if(indexof(str, ".") != -1) {
        float <- substr(str, indexof(str, "."), nchar(str))
        whole <- substr(str, 1, indexof(str, "."))
    } else {

        #If there is no decimal point, just use integer version
        return(basechange.int(str, old.base, new.base))
    }

    #Convert integer portion
    new.whole <- basechange.int(whole, old.base, new.base)
    new.whole <- paste(new.whole, ".", sep = "", collapse = "")
    new.dec <- ""
    
    if(oldbase == 10) {

        float <- as.numeric(float)

        #Convert decimal portion; keep going until there is no more 
        # decimal portion OR we have reached "acc": 
        while(nchar(new.dec) < acc) {

            float <- float*new.base
            decnum <- floor(float)
            digit <- allchars[decnum + 1]
            new.dec <- paste(new.dec, digit, sep = "", collapse = "")

            if(is.int(float))
                break
            
            float <- float - my.floor(float)
        }

        out <- paste(new.whole, new.dec, sep = "", collapse = "")
        
    } else {

        tot.digits <- nchar(float)
        decimal <- 0
        
        for(i in 2:tot.digits) {

            dig <- which(allchars == substr(float, i, i)) - 1
            decimal <- decimal + dig*oldbase^(-1*(i-1))
        }

        decimal <- substr(as.character(decimal), 3, nchar(decimal))
        base.ten <- paste(new.whole, decimal, sep = "", collapse = "")
        base.ten <- substr(base.ten, 1, acc + nchar(new.whole))
        return(basechange.float(base.ten, 10, new.base, acc))
    }

    out
}

#Helper function for the base change functions
poss.base <- function(num, base) {

    string <- toupper(as.character(num))

    for(i in 1:nchar(string)) {

        if((which(allchars == substr(string, i, i)) < 1 ||
            which(allchars == substr(string, i, i)) > base) &&
            substr(string, i, i) != ".")
            return(FALSE)
    }

    return(TRUE)
}

# Finds the nearest integer power of "base" to "num":

nearest.int.power <- function(base, num) {

    # Edge cases
    if(num < 0 && base > 0)
        return(NaN)

    # 0^(any power except 0) = 0, so the answer would be almost anything. If num is not 0
    # 0^0 = 1, 0 or undefined. To avoid confusion, output NaN if base = 0.
    if(base == 0)
        return(NaN)
    
    if(num == 0 && abs(base) > 1)
        return(-Inf)

    if(num == 0 && abs(base) < 1)
        return(Inf)

    if(num == 0)
        return(NaN)
    
    if(base == 1 || num == 1)
        return(0)
    
    if(base == -1 && num < 0)
        return(1)
    
    if(base == -1 && num > 0)
        return(0)
    
    even <- NULL
    
    if(base < 0) {

        #negative bases with positive results have even exponents
        even <- TRUE
        base <- -1*base
    }

    # "num" will only be negative if "base" is negative
    if(num < 0) {

        #negative bases with negative results have odd exponents
        even <- FALSE
        num <- -1*num
    }

    i <- 0

    #If base is positive, powers can be even or odd, so iterate by one.
    if(is.null(even))
        nextI <- 1

    else {

        #If "base" is negative, powers are all even or all odd, so iterate
        # by 2.
        nextI <- 2

        #If powers are odd, start at 1 instead of 0.
        if(!even)
            i <- 1
    }

    # "sameside" is TRUE if "num" and "base" are on the same side of 1.
    if(num > 1 && base > 1) {

        sameside <- TRUE
        
       while(base^i < num)
           i <- i + nextI
    }
    
    else if(num < 1 && base > 1){

        sameside <- FALSE
        
        while(base^i > num) 
           i <- i - nextI
    }

    else if(num < 1 && base < 1) {

        sameside <- TRUE

        while(base^i > num)
            i <- i + nextI
    }

    else {

        sameside <- FALSE

        while(base^i < num)
            i <- i - nextI
    }
    
    # Get difference between "num" and nearest power above and below
    # Needed to decide which power is closest.
    if(sameside) {

        diffLower <- abs(num - base^(i-nextI))       
        diffUpper <- abs(num - base^i)
    }

    else {
       
       diffLower <- abs(num - base^i)
       diffUpper <- abs(num - base^(i+nextI))
    }
    
    if(diffLower < diffUpper && sameside)
        return(i - nextI)
    
    else if(diffUpper < diffLower && !sameside)
        return(i + nextI)

    else
        return(i)
}

# Generates the "n"th value in a recursive sequence.
# Arguments:
# start: first element in sequence
# n: nth element (will be a fixed number if sequence converges)
# seq.func: function determining sequence's recursive rule
# make.list: TRUE to generate all elements of sequence between start and n;
#            FALSE otherwise
# maxitt: max number of iterations.
# supp.warn: TRUE to suppress warnings of sequence divergence; FALSE otherwise
# tol: tolerance level; used to determine convergence
sequence.rec <- function(start, n = Inf, seq.func, make.list = FALSE, maxitt = 5000, supp.warn = FALSE, tol = EPSILON) {

    if(is.infinite(n)) 
        n = maxitt
    else 
        maxitt = n
    
    
    i <- 1
    
    dist <- numeric(2)
    new <- 0

    if(make.list) {
        list <- numeric(0)
        list[1] <- start
    }
    
    while(i < maxitt) {

        new <- seq.func(start)
        
        if(is.infinite(new))
            break

        # Attempt at eliminating unnecessary iterations if the sequence
        # appears to converge
        if(abs(new - start) < tol && length(list) > 3)
            break

        start <- new
        
        i <- i + 1

        if(make.list)
            list[i] <- new
    }

    len = length(list)

    # Tries to figure out if sequence is diverging.
    if(!supp.warn && (is.infinite(new) || (make.list &&
        abs(list[len] - list[len-1]) >= abs(list[len-1] - list[len-2]))))
        warning("Possible divergence")

    if(make.list)
        out <- list(ans = new, list = list)
    else
        out <- new

    out
}
