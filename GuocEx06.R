GaussJordan <- function(mat){                   # function that uses gauss jordan elimination to solve for a system of linear equations
  A <- mat                                      # the function accepts an augmented coefficient matrix
  n <- nrow(A)
  
  for (i in 1:n){
    if (i != n){
      pr <- which.max(abs(A[i:n, i])) + (i - 1)  # find pivot row
      if(A[pr,i] == 0){                          # if element is equal to 0, no unique solution exists
        return (NA)
      } else {
        A = swap(A, pr, i)                       # swap pivot row and ith row
      }
    }
    A[i,] <- A[i,] / A[i,i]                      # normalize row
    for (j in 1:n){              
      if (i == j){                               # will form diagonal
        next
      }
      A[j,] <- A[j,] - (A[j,i] * A[i,])
    }
  }
  # the last column of the resulting matrix will contain the solution
  return (list(augcoeff = A, coefficients = A[1:n, ncol(A)]))
}


PolynomialRegression <- function(n, mat){                   # function that takes an integer(determining the order of the polynomial) and a matrix of two columns containing x and y variables
  augcoeffmatrix <- matrix(nrow = n + 1, ncol = n + 2, dimnames = list(c(1:(n + 1)), c(1:(n + 1), 'RHS')))         # create empty matrix with n + 1 rows, n + 1 + rhs columns
  sumVec <- c()                                             # vector to store the sum of all Xi^0 - Xi^2n
  rhsVec <- c()                                             # vector to store the sum of (Xi^n)Yi for the values of the right hand side
  for(i in 0:(2 * n)){                        
    cSum <- 0
    rhsSum <- 0
    for(j in 1:nrow(mat)){                                  # for each value of x
      cSum <- cSum + ((mat[j, 1]) ^ i)                      # adds Xj^i to the sum
      rhsSum <- rhsSum + ((mat[j, 1]) ^ i) * mat [j, 2]     # same as cSum but is multiplied by Yj
    }
    rhsVec[i + 1]<- rhsSum                                  # adds value to vector
    sumVec[i + 1] <- cSum
  } 
  x = c()                                                   # vector to store x^n for column names of returned matrix
  for(i in 1:(n + 1)){                                      # place the values(sums) in the augmented coefficient matrix
    for(j in 1:(n + 1)){                                    
      augcoeffmatrix[i, j] <- sumVec[i + j - 1]             # [i + j - 1] increments index but also ensures the first number to be taken for each row is always shifted one spot to the right
      x[i] = paste(c('x^', i - 1), collapse = '')           
    }
  }
  for(i in 1:nrow(augcoeffmatrix)){                         # place the values for the right hand side 
    augcoeffmatrix[i, (n + 2)] <- rhsVec[i]                 # does not use all elements in the vector but fills up all rows
  }
  
  print(mat)
  # print(list(augcoeffmatrix = augcoeffmatrix, polynomial = polynom))
  polynom <- matrix(GaussJordan(augcoeffmatrix)$coefficients, nrow = 1, dimnames = list(NULL, x))        # matrix containing the results of the GaussJordan function to find the unknowns given the augmented coefficient matrix
  return(list(augcoeffmatrix = augcoeffmatrix, polynomial = polynom))                                    # return labeled list containing original augmented coefficient matrix, and matrix containing the coefficients of the polynomial
}
