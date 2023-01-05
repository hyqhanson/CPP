#ifndef RESSPMV_H
#define RESSPMV_H

#include "mkl.h"


subroutine dgetrf	(	integer 	M,
integer 	N,
double precision, dimension( lda, * ) 	A,
integer 	LDA,
integer, dimension( * ) 	IPIV,
integer 	INFO 
)	

float *A;
int64  n;

dgetrf(m,n, A, lda, ipiv, &info )




template <typename INT, typename A_TYPE>
void dgetrf(INT m, INT n, A_TYPE *A, INT lda, INT *ipiv, INT *inf)

{ //  r = b- A*x
    for (INT row = 0; row < nrows; row++)
    {
        PREC sum_row = 0.0;
        // row * x
        for (INT j = row_ptr[row] - one_based; j < row_ptr[row + 1] - one_based; j++)
        {
            INT col = col_idx[j] - 1;
            sum_row += PREC(values[j]) * PREC(x[col]);
        }
        // printf("%g   %g \n", b[row], sum_row);
        residual[row] = RES_TYPE(((PREC)(b[row]) - sum_row));
    }
}




template <typename PREC, int one_based, typename INT, typename VAL_TYPE, typename X_TYPE, typename B_TYPE,
          typename RES_TYPE>
void ResidualCSR(INT nrows, const VAL_TYPE *values, const INT *row_ptr,
                 const INT *col_idx, const X_TYPE *x, const B_TYPE *b, RES_TYPE *residual)

{ //  r = b- A*x
    for (INT row = 0; row < nrows; row++)
    {
        PREC sum_row = 0.0;
        // row * x
        for (INT j = row_ptr[row] - one_based; j < row_ptr[row + 1] - one_based; j++)
        {
            INT col = col_idx[j] - 1;
            sum_row += PREC(values[j]) * PREC(x[col]);
        }
        // printf("%g   %g \n", b[row], sum_row);
        residual[row] = RES_TYPE(((PREC)(b[row]) - sum_row));
    }
}

template <typename INT, typename T, int one_based>
void SpMV(INT nrows, const T *values, const INT *row_ptr,
          const INT *col_idx, const T *x, T *b)

{ //  r = b- A*x
    for (INT row = 0; row < nrows; row++)
    {
        T sum_row = 0.0;
        // row * x
        for (INT j = row_ptr[row] - one_based; j < row_ptr[row + 1] - one_based; j++)
        {
            INT col = col_idx[j] - 1;
            sum_row += values[j] * x[col];
        }
        // printf("%g   %g \n", b[row], sum_row);
        b[row] = sum_row;
    }
}

double
relative_error(int n, double *x, double *xe)
{
    double *d = (double *)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++)
        d[i] = x[i] - xe[i];
    return cblas_dnrm2(n, d, 1) / cblas_dnrm2(n, xe, 1);
}

#define PRINT_TIME(str, time) \
    printf("TIMING: %-30s  %.2f s\n", str, time);
#define PRINT_TIME_ITER(str, time, iter) \
    printf("TIMING: %-30s  %.2f s   iterations %d\n", str, time, iter);

// Shahrooz
#define PRINT_TIME_BOTH_ITERS(str, time, pardiso_iter, ned_iter) \
    printf("TIMING: %-30s  %.2f s   Pardiso iters %d, Ned's iters %d\n", str, time, pardiso_iter, ned_iter);

#define BLACK "\033[30m"              /* Black */
#define RED "\033[31m"                /* Red */
#define GREEN "\033[32m"              /* Green */
#define YELLOW "\033[33m"             /* Yellow */
#define BLUE "\033[34m"               /* Blue */
#define MAGENTA "\033[35m"            /* Magenta */
#define CYAN "\033[36m"               /* Cyan */
#define WHITE "\033[37m"              /* White */
#define BOLDBLACK "\033[1m\033[30m"   /* Bold Black */
#define BOLDRED "\033[1m\033[31m"     /* Bold Red */
#define BOLDGREEN "\033[1m\033[32m"   /* Bold Green */
#define BOLDYELLOW "\033[1m\033[33m"  /* Bold Yellow */
#define BOLDBLUE "\033[1m\033[34m"    /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m" /* Bold Magenta */
#define BOLDCYAN "\033[1m\033[36m"    /* Bold Cyan */
#define BOLDWHITE "\033[1m\033[37m"   /* Bold White */
#define RESET "\x1B[0m"
#define BOLD "\x1B[1m" // Bold Text Formula...

#define PRINT_TOTAL_TIME(time) \
    printf("TIMING: total    %-30s  %.2f %ss\n", BOLDRED, time, RESET);

#define PRINT_DOUBLE(str, d) \
    printf("  %-30s %s %e %s \n", str, BOLDBLUE, d, RESET);
#endif