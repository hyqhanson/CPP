#include <iostream>
#include <stdexcept>
#include <vector>
#include <string>
#include <algorithm>
#include <random>
using namespace std;

template <typename INT, typename FLOAT>
class denseM
{

public:
    denseM(const INT &, const INT &);
    denseM(const INT &, const INT &, const vector<FLOAT> &);
    denseM(const vector<FLOAT> &);

    template <typename INT_alt, typename FLOAT_alt>
    denseM(const denseM<INT_alt, FLOAT_alt> &);

    const vector<FLOAT> &get_data() const
    {
        return matrix_;
    }

    INT get_num_rows() const
    {
        return rows_;
    }

    INT get_num_cols() const
    {
        return cols_;
    }

    bool is_symmetric() const;

    // check the value for a given row and col
    const FLOAT at(const INT &, const INT &) const;

    // check and modify the value in the 1d vector for a given index
    const FLOAT at(const INT &) const;

    // overloading [] for checking the value in the 1d vector for a given index, allows modification
    FLOAT &operator[](INT);

    const FLOAT &operator[](INT) const;

    /**
     * @brief Exception occurs when given matrix vector mismatch the given size
     * of rows and columns
     *
     */
    class invalid_size : public invalid_argument
    {
    public:
        invalid_size() : invalid_argument("The matrix size can only be positive."){};
    };

    /**
     * @brief Exception occurs when given matrix vector mismatch the given size
     * of rows and columns
     *
     */
    class size_mismatch : public invalid_argument
    {
    public:
        size_mismatch() : invalid_argument("The matrix vector does not have the same amount of elements as the given size."){};
    };

private:
    INT rows_ = 0;
    INT cols_ = 0;
    vector<FLOAT> matrix_;
};

// preallocate memory
template <typename INT, typename FLOAT>
denseM<INT, FLOAT>::denseM(const INT &rows, const INT &cols)
    : rows_(rows), cols_(cols)
{
    matrix_.resize(cols * rows);
}

template <typename INT, typename FLOAT>
denseM<INT, FLOAT>::denseM(const INT &rows, const INT &cols, const vector<FLOAT> &matrix)
    : rows_(rows), cols_(cols), matrix_(matrix) {}

template <typename INT, typename FLOAT>
denseM<INT, FLOAT>::denseM(const vector<FLOAT> &diagonal)
    : rows_((INT)(diagonal.size())), cols_((INT)(diagonal.size()))
{

    if (diagonal.size() == 0)
    {
        throw invalid_size();
    }
    matrix_.resize(diagonal.size() * diagonal.size(), 0);
    for (INT i = 0; i < diagonal.size(); i++)
    {
        matrix_[cols_ * i + i] = diagonal[i];
    }
}

template <typename INT, typename FLOAT>
template <typename INT_alt, typename FLOAT_alt>
denseM<INT, FLOAT>::denseM(const denseM<INT_alt, FLOAT_alt> &A)
    : rows_((INT)A.get_num_rows()), cols_((INT)A.get_num_cols())
{
    vector<FLOAT_alt> temp = A.get_data();
    for (FLOAT_alt &i : temp)
    {
        matrix_.push_back((FLOAT)i);
    }
}

template <typename INT, typename FLOAT>
bool denseM<INT, FLOAT>::is_symmetric() const
{
    // Matrix must be a square matrix
    if (rows_ != cols_)
    {
        return 0;
    }
    for (INT i = 0; i < rows_; i++)
    {
        for (INT j = 0; j < cols_; j++)
        {
            if (matrix_[i * cols_ + j] != matrix_[j * cols_ + i])
            {
                return 0;
            }
        }
    }
    return 1;
}

template <typename INT, typename FLOAT>
const FLOAT denseM<INT, FLOAT>::at(const INT &row, const INT &col) const
{
    if (row <= 0 || col <= 0)
    {
        throw invalid_size();
    }
    return matrix_[row * cols_ + col];
}

template <typename INT, typename FLOAT>
const FLOAT denseM<INT, FLOAT>::at(const INT &index) const
{
    if (index >= cols_ * rows_ || index <= 0)
    {
        throw invalid_size();
    }
    return matrix_[index];
}

template <typename INT, typename FLOAT>
FLOAT &denseM<INT, FLOAT>::operator[](INT index)
{
    if (index >= cols_ * rows_ || index < 0)
    {
        throw invalid_size();
    }
    return matrix_[index];
}

template <typename INT, typename FLOAT>
const FLOAT &denseM<INT, FLOAT>::operator[](INT index) const
{
    if (index >= cols_ * rows_ || index < 0)
    {
        throw invalid_size();
    }
    return matrix_[index];
}

template <typename INT, typename FLOAT>
ostream &operator<<(ostream &out, const denseM<INT, FLOAT> &M)
{
    vector<FLOAT> temp = M.get_data();
    for (uint64_t i = 0; i < temp.size(); i++)
    {
        out << temp[i] << " ";
        if ((i + 1) % M.get_num_cols() == 0)
        {
            out << "\n";
        }
    }
    return out;
}

template <typename INT, typename FLOAT>
denseM<INT, FLOAT> operator+(const denseM<INT, FLOAT> &M1, const denseM<INT, FLOAT> &M2)
{
    // Exception occurs when given matrices have different size
    if (M1.get_num_cols() != M2.get_num_cols() || M1.get_num_rows() != M2.get_num_rows())
    {
        throw invalid_argument("These two matrices are not able to add each other.");
    }

    INT Row = M1.get_num_rows();
    INT Col = M1.get_num_cols();
    denseM<INT, FLOAT> sum(Row, Col);

    for (INT i = 0; i < Row * Col; i++)
    {
        sum[i] += M1[i] + M2[i];
    }
    return sum;
}

template <typename INT, typename FLOAT>
denseM<INT, FLOAT> operator-(const denseM<INT, FLOAT> &M1, const denseM<INT, FLOAT> &M2)
{
    // Exception occurs when given matrices have different size
    if (M1.get_num_cols() != M2.get_num_cols() || M1.get_num_rows() != M2.get_num_rows())
    {
        throw invalid_argument("These two matrices are not able to subtract each other.");
    }

    INT Row = M1.get_num_rows();
    INT Col = M1.get_num_cols();

    denseM<INT, FLOAT> sub(Row, Col);

    for (INT i = 0; i < Row * Col; i++)
    {
        sub[i] += M1[i] - M2[i];
    }

    return sub;
}

template <typename INT, typename FLOAT>
denseM<INT, FLOAT> operator*(const denseM<INT, FLOAT> &M1, const denseM<INT, FLOAT> &M2)
{
    // Exception occurs when given matrices are not able to multiply each other
    if (M1.get_num_cols() != M2.get_num_rows())
    {
        throw invalid_argument("These two matrices are not able to multiply each other.");
    }

    INT row1 = M1.get_num_rows();
    INT row2 = M2.get_num_rows();
    INT col1 = M1.get_num_cols();
    INT col2 = M2.get_num_cols();

    // size of the multiplication result is rows of matrix1 with cols of matrix2
    denseM<INT, FLOAT> result(row1, col2);

    for (uint64_t i = 0; i < row1; i++)
    {
        for (uint64_t j = 0; j < col2; j++)
        {
            for (uint64_t k = 0; k < row2; k++)
            {
                result[col2 * i + j] += M1[col1 * i + k] * M2[col2 * k + j];
            }
        }
    }

    return result;
}

/*
// infinity norm of denseM
template <typename INT, typename FLOAT>
FLOAT norm(const denseM<INT, FLOAT> &_M)
{
    INT Row = _M.get_row();
    INT Col = _M.get_col();
    vector<FLOAT> M = _M.get_vec();
    FLOAT maxRowSum = 0.0;

    for (INT i = 0; i < Row; i++)
    {
        FLOAT RowSum = 0.0;
        for (INT j = 0; j < Col; j++)
        {
            RowSum += fabs(M[i * Col + j]);
        }
        maxRowSum = max(maxRowSum, RowSum);
    }

    return maxRowSum;
}

template <typename INT, typename FLOAT>
void LU_withPivot(const denseM<INT, FLOAT> &_A,
                  denseM<INT, FLOAT> &_P, denseM<INT, FLOAT> &_L, denseM<INT, FLOAT> &_U)
{
    // Exception occurs when given matrix is not a square matrix
    if (_A.get_col() != _A.get_row())
    {
        throw invalid_argument("This matrix is not a square matrix.");
    }

    INT size = _A.get_row();
    vector<FLOAT> A = _A.get_vec();

    // Create an identity matrix for initializing permutation matrix
    vector<FLOAT> P(size * size);
    for (INT i = 0; i < size; i++)
    {
        P[size * i + i] = 1;
    }

    for (INT i = 0; i < size; i++)
    {
        INT current_pivot = i;

        // compare the pivot in the current row with the number in
        // the following rows
        for (INT j = i + 1; j < size; j++)
        {
            // find the biggest absolute value for each column, choose it
            // as the new current pivot
            if (fabs(A[j * size + i]) > fabs(A[current_pivot * size + i]))
            {
                current_pivot = j;
            }
        }

        // If pivot changed, swap rows
        if (current_pivot != i)
        {
            for (INT n = 0; n < size; n++)
            {
                swap(A[i * size + n], A[current_pivot * size + n]);
                swap(P[i * size + n], P[current_pivot * size + n]);
            }
        }

        for (INT j = i + 1; j < size; j++)
        {
            A[j * size + i] /= A[i * size + i];
            for (INT k = i + 1; k < size; k++)
            {
                A[j * size + k] -= A[i * size + k] * A[j * size + i];
            }
        }
    }

    // Pass the permutation matrix to _P
    denseM<INT, FLOAT> tempP(size, size, P);
    _P = tempP;

    // Separate lower-triangular matrix out
    // Create an identity matrix since lower-triangular matrix has diagonal with 1s
    vector<FLOAT> lower(size * size);
    for (INT i = 0; i < size; i++)
    {
        lower[size * i + i] = 1;
    }

    for (INT i = 0; i < size; i++)
    {

        for (INT j = 0; j < size; j++)
        {
            if (j < i)
            {
                lower[i * size + j] = A[i * size + j];
            }
        }
    }
    denseM<INT, FLOAT> Lower(size, size, lower);
    _L = Lower;

    // Separate upper-triangular matrix out
    vector<FLOAT> upper;
    for (INT i = 0; i < size; i++)
    {
        for (INT zero = 0; zero < i; zero++)
        {
            upper.push_back(0);
        }
        for (INT j = i; j < size; j++)
        {
            upper.push_back(A[i * size + j]);
        }
    }
    denseM<INT, FLOAT> Upper(size, size, upper);
    _U = Upper;
}

template <typename INT, typename FLOAT>
denseM<INT, FLOAT> forwardSub(const denseM<INT, FLOAT> &_L, const denseM<INT, FLOAT> &_b)
{
    // Exception occurs when given matrix L is not a square matrix,
    // or given matrix b does not have the same rows as L
    if (_b.get_row() != _L.get_col())
    {
        throw invalid_argument("b should have the same row number as L.");
    }
    else if (_L.get_col() != _L.get_row())
    {
        throw invalid_argument("L is not a square matrix.");
    }

    INT size = _L.get_col();
    vector<FLOAT> L = _L.get_vec();

    vector<FLOAT> b = _b.get_vec();
    // initializing the solution x
    vector<FLOAT> x(size, 0);

    for (INT i = 0; i < size; i++)
    {
        FLOAT temp = 0;
        // temp will store the summation of all known xi value times corresponding
        // L value in each row
        for (INT j = 0; j <= i; j++)
        {
            temp += x[j] * L[i * size + j];
        }
        // b subtracts temp will left out the value of the current xi times
        // corresponding L value, divide this L number, we can get xi
        x[i] = (b[i] - temp) / L[i * size + i];
    }
    denseM<INT, FLOAT> _x(size, 1, x);
    return _x;
}

template <typename INT, typename FLOAT>
denseM<INT, FLOAT> backwardSub(const denseM<INT, FLOAT> &_U, const denseM<INT, FLOAT> &_b)
{
    // Exception occurs when given matrix U is not a square matrix,
    // or given matrix b does not have the same rows as U
    if (_b.get_row() != _U.get_col())
    {
        throw invalid_argument("b should have the same row number as U.");
    }
    else if (_U.get_col() != _U.get_row())
    {
        throw invalid_argument("U is not a square matrix.");
    }

    INT size = _U.get_col();
    vector<FLOAT> U = _U.get_vec();
    vector<FLOAT> b = _b.get_vec();
    // initializing the solution x
    vector<FLOAT> x(size, 0);

    // backward substituion starts from the last row
    for (INT i = size; i >= 1; i--)
    {
        FLOAT temp = 0;
        // temp will store the summation of all known xi value times corresponding
        // U value in each row
        for (INT j = i - 1; j < size; j++)
        {
            temp += x[j] * U[(i - 1) * size + j];
        }
        // b subtracts temp will left out the value of the current xi times
        // corresponding U value, divide this U number, we can get xi
        x[i - 1] = (b[i - 1] - temp) / U[(i - 1) * size + (i - 1)];
    }
    denseM<INT, FLOAT> _x(size, 1, x);
    return _x;
}

template <typename INT, typename FLOAT>
denseM<INT, FLOAT> LU_solver(const denseM<INT, FLOAT> &_L, const denseM<INT, FLOAT> &_U, const denseM<INT, FLOAT> &_P, const denseM<INT, FLOAT> &_b)
{
    // Starts with LUx = Pb
    // Let y = Ux, then Ly = Pb, solve by forward substituion
    denseM<INT, FLOAT> y = forwardSub(_L, _P * _b);
    // Solve Ux = y by backward substituion
    denseM<INT, FLOAT> x = backwardSub(_U, y);

    return x;
}

template <typename FACT, typename RES, typename INT, typename FLOAT>
denseM<INT, FLOAT> IR(const denseM<INT, FLOAT> &_A, const denseM<INT, FLOAT> &_b)
{
    denseM<INT, FACT> Af(A);
    // permutation matrix from LU-decomposition
    denseM<INT, FLOAT> P;
    // lower triangular matrix from LU-decomposition
    denseM<INT, FLOAT> L;
    // upper triangular matrix from LU-decomposition
    denseM<INT, FLOAT> U;
    // Get P, L, U
    LU_withPivot(Af, P, L, U);

    // L,U must be in FACT precision
    // LU_solver returns x in FACT precision
    denseM<INT, FLOAT> x = LU_solver(L, U, P, _b);

    INT max_iter = 1000;
    INT iter = 0;

    // residual
    denseM<INT, RES> r;
    // infinity norm of r
    FLOAT residual = 1;
    // tolerance for stopping iteration
    FLOAT tol = 1e-10;

    // correction
    denseM<INT, FLOAT> c;
    while (iter != max_iter && residual > tol && residual != 0)
    {
        r = _b - (_A * denseM<INT,RES>(x);
        residual = norm<INT, FLOAT>(r);

        // for now, using LU get correction
        c = LU_solver(L, U, P, r);
        x = x + c;
        iter++;
    }

    return x;
}

denseM<int, double> IR<float, double>(_A, _b)
*/