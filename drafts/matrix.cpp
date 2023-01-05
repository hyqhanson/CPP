
#include <iostream>
#include <vector>

template <template <class, class> class MATRIX, class INT, class FLOAT>
void LU(MATRIX<INT, FLOAT> &A, std::vector<INT> &pivot, int &exitflag)
{
  // https://github.com/cayek/TESS3/blob/master/tess3r/src/external/lapack/dgetrf.c
  A.print();
}

template <typename INT, typename FLOAT>
class DenseMatrix
{
public:
  DenseMatrix(){};

  DenseMatrix(INT m, INT n) : num_rows_(m), num_cols_(n)
  { // allocate memory
  }
  DenseMatrix(INT m, INT n, FLOAT *values) {}

  template <typename INT1, typename FLOAT1>
  DenseMatrix(const DenseMatrix<INT1, FLOAT1> &A){};

  INT get_num_rows() const { return num_rows_; };
  const FLOAT &get_ij(INT i, INT j) const { return values_(i * num_cols_ + j); }
  FLOAT &set_ij(INT i, INT j, FLOAT &v) const
  {
    return values_[i * num_cols_ + j] = v;
  }

  /*
    class Symmetric: public DenseMatrix{


    }
    */

  /*
    class SPD: public DenseMatrix{

  // Factorization is Cholesky
    }
    */

  // void LU(DenseMatrix & LU, std::vector<INT> & pivot,  int &exitflag){
  // LU(*this, pivot, exitflag);}

  // overload << for output
  // overload >> for input
  void print() const
  {
    std::cout << "Dense Matrix ";
    std::cout << "INT is " << typeid(INT).name() << "  FLOAT is "
              << typeid(FLOAT).name() << std::endl;
    // output elements
  }

private:
  INT num_rows_, num_cols_;
  FLOAT *values_;
};

// y := alpha*A*x + beta*y,
template <template <class, class> class MATRIX, class INT, class FLOAT,
          typename PREC>
void GEMV(const FLOAT &alpha, const MATRIX<INT, FLOAT> &A, const FLOAT &beta,
          std::vector<FLOAT> &y)
{
  // cblas_gemv
}

// r = b-A*X

template <template <class, class> class MATRIX, class INT, class FLOAT,
          typename PREC>
void Residual(const MATRIX<INT, FLOAT> &A, const FLOAT &x,
              const std::vector<FLOAT> &b, std::vector<FLOAT> &res)
{
  // y := alpha*A*x + beta*y,
  res = b;
  GEMV<PREC>(-1, A, x, 1, res);
}

template <typename INT, typename FLOAT>
class CSRMatrix
{
public:
  void print() const
  {
    std::cout << "CSR Matrix ";
    std::cout << "INT is " << typeid(INT).name() << "  FLOAT is "
              << typeid(FLOAT).name() << std::endl;
  }

private:
  INT *row_ptrs_;
  INT *col_indices;
};

// P*L*U = A
// A*x = bto

template <template <class, class> class MATRIX, class INT, class FLOAT>
class LU_Solver
{
public:
  template <typename INT1, typename FLOAT1>
  LU_Solver(const MATRIX<INT1, FLOAT1> &A){};
  void Factorize(){};
  void Solve(std::vector<FLOAT> &b){};

public:
  MATRIX<INT, FLOAT> A;
  MATRIX<INT, FLOAT> LU;
  std::vector<INT> pivot_;
};

template <template <class, class> class MATRIX, class INT, class FLOAT>
void LowerTriangularSolve(const MATRIX<INT, FLOAT> &LU,
                          const std::vector<FLOAT> &b, std::vector<FLOAT> &x)
{
  // Solve L*x = b, L lower triangula
}

template <template <class, class> class MATRIX, class INT, class FLOAT>
void UpperTriangularSolve(const MATRIX<INT, FLOAT> &LU,
                          const std::vector<FLOAT> &b, std::vector<FLOAT> &x)
{
  // Solve U*x = b, U upper triangular
}

template <template <class, class> class MATRIX, class INT, class FLOAT>
void SolveWithLU(const MATRIX<INT, FLOAT> &LU, const std::vector<FLOAT> &b,
                 std::vector<FLOAT> &x)
{
  // call the above two functions. / Solve U*x = b, U upper triangular
}

template <template <class, class> class MATRIX, class INT, class FLOAT>
void Solve(const MATRIX<INT, FLOAT> &A, const std::vector<FLOAT> &b,
           std::vector<FLOAT> &x)
{
  // Do LU factorization and call SolveWithLU
  // call the above two functions. / Solve U*x = b, U upper triangular
}

template <template <class, class> class MATRIX, class INT, class FLOAT,
          typename FACT_PREC>
void Solve_IR(const MATRIX<INT, FLOAT> &A, const std::vector<FLOAT> &b,
              std::vector<FLOAT> &x)
{
  MATRIX<INT, FACT_PREC> Af(A), LUf;
  std::vector<FACT_PREC> bf, xf, d, r;
  std::vector<INT> pivot;
  INT info;

  LU(Af, pivot, info);

  SolveWithLU(LUf, bf, x);

  // loop
  // r = b-A*x
  // LU*d = r
  // SolveWithLU(LUf, r, d);
  // x = x+r
}

int main()
{
  DenseMatrix<int, float> A, lu;

  LU_Solver<DenseMatrix, int, double> Solver(A);
  std::vector<double> b;
  Solver.Factorize();
  // Solve Ax= b
  Solver.Solve(b);

  DenseMatrix<int, double> B(A);

  std::vector<int> pivot;
  int info;
  LU(A, pivot, info);

  // A.LU(lu);

  DenseMatrix<int, double> Ad, lud;

  LU(Ad, pivot, info);

  CSRMatrix<int, double> Acsr, LUcsr;
  LU(Acsr, pivot, info);
}
