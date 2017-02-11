#include <R.h>
#include <Rmath.h>

/*"tri" represents the matrix to be manipulated through row operations.
  This function is named "trinonsq" because it was originally intended only for square matrices,
  but was then modified to apply to nonsquare matrices too. */
void trinonsq(int *rows, int *cols, int *origcols, int *decp, int *upper, int *special,
	      int *diagonal, double *tri, double *pmat, double *diagmat) {

  int rows_, cols_, upper_,special_, rowA, opp, diagRow, diffRow, col, rowB, nonsq_shift;
  double diag, old, new, nume, coeff;
  char *str;
  
  rows_ = *rows;
  cols_ = *cols;
  upper_ = *upper;
  special_ = *special;

  /*This is necessary for nonsquare matrices because their "diagonals" look strange compared to
    square matrices. For square matrices, this variable = 0. */
  nonsq_shift = *origcols - rows_;

  /*This is the "grand loop" - the upper or lower triangular matrix, along with its diagonal
    matrix (if "special" is true) and permutation matrix EACH get modified together */
  for(rowA = 0; rowA < rows_; rowA++) {

    /*If "rowA" is the first row, then "opp" is the last row; if "rowA" is the second row,
    "opp" is the second-to-last row, and so on.*/
    opp = rows_ - 1 - rowA;

    //"diagRow" is the current row.
    diagRow = upper_ ? rowA : opp;

    //"diag" is the diagonal element in "diagRow"
    diag = upper_ ? tri[diagRow*cols_ + diagRow] : tri[diagRow*cols_ + diagRow + nonsq_shift];

    /* This rearranges rows in a matrix if a value on the diagonal = 0. (This forms the permutation
     * matrix of the "LU" and "LDV" decompositions.*/
    if(diag == 0) {

      diffRow = (upper_) ? (diagRow + 1) : (diagRow - 1);
      int broke = 0;

      //"diffRow" is the row that will be swapped with the row containing 0 as a diagonal element.
      while((upper_ && diffRow <= rows_ - 1) || (!upper_ && diffRow >= 0) {

	//Break once a row is found that does not have 0 on the diagonal.
	if(tri[diffRow*cols_ + diagRow] != 0) {
	  broke = 1;
	  break;
	}
	
	if(upper_) {
	  diffRow++;
	} else {
	  diffRow--;
	}
      }

      if(broke) {
	for(col = 0; col < cols_; col++) {
	    old = tri[diagRow*cols_ + col];
	    new = tri[diffRow*cols_ + col];
	    tri[diagRow*cols_ + col] = new;
	    tri[diffRow*cols_ + col] = old;     
	  }
      }
      
      /*Permutation matrix only modified if decomposition specified. Errors
      * may occur if this matrix is modified and the intended triangular matrix
      * is formed from an augmented matrix*/
      if(*decp && broke) {
	for(col = 0; col < rows_; col++) {
	  double oldP = pmat[diagRow*rows_ + col];
	  double newP = pmat[diffRow*rows_ + col];
	  pmat[diagRow*rows_ + col] = newP;
	  pmat[diffRow*rows_ + col] = oldP; 
	}
      }
    }

    //If special upper or lower specified, the diagonal matrix is modified.
    if(special_) {
      diag = upper_ ? tri[rowA*cols_ + rowA] : tri[opp*cols_ + opp + nonsq_shift];

      if(*diagonal)
	diagmat[diagRow*rows_ + diagRow] = diag;
      
      for(col = 0; col < cols_; col++) {
	      
	if(diag != 0)
	  tri[diagRow*cols_ + col] /= diag;
      }
    }

    //"rowB" is the row to have row operations performed on it with respect to the values in "rowA"
    rowB = upper_ ? rowA : opp;
    diag = upper_ ? tri[rowA*cols_ + rowA] : tri[opp*cols_ + opp + nonsq_shift];

    if(diag == 0)
      continue;
    
    while((upper_ && rowB < (rows_ - 1)) || (!upper_ && rowB > 0)) {	
      
      if(upper_)
	rowB++;
      else
	rowB--;

      //"nume" is the numerator of the coefficient 
      nume = upper_ ? tri[rowB*cols_ + rowA] : tri[rowB*cols_ + opp + nonsq_shift];

      if(nume == 0)
	continue;

      coeff = -1*nume/diag;

      //Row operation performed on each element of a given row ("rowB")
      for(col = 0; col < cols_; col++) {

	if(upper_) {
	  tri[rowB*cols_ + col] = coeff*tri[rowA*cols_ + col] 
	    + tri[rowB*cols_ + col];

	} else {
	  tri[rowB*cols_ + col] = tri[rowB*cols_ + col] 
	    + coeff*tri[opp*cols_ + col];

	} //end else
      } // end for "col" loop
    } // end while loop doing row operations with respect to given "rowA"
  } //end grand loop iterating over "rowA"     
} //end func

	
	
	  

	  
	  
	  
	  

	    
	    
	      
	      

	      
	

     
      
