#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

SEXP sample_beta(SEXP n, SEXP pL, SEXP XL, SEXP xL2, SEXP bL, SEXP e, SEXP varBj, SEXP varE, SEXP minAbsBeta)
{
        double *xj, *pXL, *pxL2, *pbL, *pe, *pvarBj;
	double rhs,c,sigma2e, smallBeta;
	
        int j,i, rows, cols;
        SEXP list;
	
  	GetRNGstate();
	
	rows=INTEGER_VALUE(n);
	cols=INTEGER_VALUE(pL);
	sigma2e=NUMERIC_VALUE(varE);
	smallBeta=NUMERIC_VALUE(minAbsBeta);
	
  	PROTECT(XL=AS_NUMERIC(XL));
        pXL=NUMERIC_POINTER(XL);

        PROTECT(xL2=AS_NUMERIC(xL2));
        pxL2=NUMERIC_POINTER(xL2);

        PROTECT(bL=AS_NUMERIC(bL));
        pbL=NUMERIC_POINTER(bL);

        PROTECT(e=AS_NUMERIC(e));
        pe=NUMERIC_POINTER(e);

        PROTECT(varBj=AS_NUMERIC(varBj));
        pvarBj=NUMERIC_POINTER(varBj);

        xj=(double *) R_alloc(rows,sizeof(double));

        for(j=0; j<cols;j++)
        {
	  rhs=0;
	  for(i=0; i<rows; i++)
	  {
	    xj[i]=pXL[i+j*rows];
	    pe[i] = pe[i] + pbL[j]*xj[i];
	    rhs+=xj[i]*pe[i];
	  }
	  rhs=rhs/sigma2e;
  	  c=pxL2[j]/sigma2e + 1.0/pvarBj[j];
	  pbL[j]=rhs/c + sqrt(1.0/c)*norm_rand();
	  
	  for(i=0; i<rows; i++)
	  {
	    pe[i] = pe[i] - pbL[j]*xj[i];
	  }
          if(fabs(pbL[j])<smallBeta)
          {
             pbL[j]=smallBeta;
          }
        }
        
        // Creating a list with 2 vector elements:
	PROTECT(list = allocVector(VECSXP, 2));
	// attaching bL vector to list:
	SET_VECTOR_ELT(list, 0, bL);
	// attaching e vector to list:
	SET_VECTOR_ELT(list, 1, e);

  	PutRNGstate();

  	UNPROTECT(6);

  	return(list);
}
