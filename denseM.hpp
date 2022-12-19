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
private:
    INT rows = 0;
    INT cols = 0;
    vector<FLOAT> matrix;

public:
    denseM()
    {
        matrix = {0};
    }

    denseM(const INT &_rows, const INT &_cols)
    {
        rows = _rows;
        cols = _cols;

        if (rows <= 0 || cols <= 0)
        {
            throw invalid_size();
        }

        random_device rd;
        mt19937 mt(rd());
        uniform_real_distribution<FLOAT> uid(-50, 50);
        for (INT i = 0; i < _rows * _cols; i++)
        {
            matrix.push_back(uid(rd));
        }
    }

    denseM(const INT &_rows, const INT &_cols, const vector<FLOAT> &_matrix)
    {
        rows = _rows;
        cols = _cols;

        if (rows <= 0 || cols <= 0)
        {
            throw invalid_size();
        }

        if (_matrix.size() != rows * cols)
        {
            throw size_mismatch();
        }
        matrix = _matrix;
    }

    denseM(const vector<FLOAT> &diagonal)
    {
        rows = (INT)(diagonal.size());
        cols = (INT)(diagonal.size());

        if (diagonal.size() == 0)
        {
            throw invalid_size();
        }
        matrix.resize(diagonal.size() * diagonal.size(), 0);
        for (INT i = 0; i < diagonal.size(); i++)
        {
            matrix[cols * i + i] = diagonal[i];
        }
    }

    // change precision
    template <typename INT_alt, typename FLOAT_alt>
    denseM(const denseM<INT_alt, FLOAT_alt> &A)
    {
        rows = (INT)A.get_row();
        cols = (INT)A.get_col();

        vector<FLOAT_alt> temp = A.get_vec();
        for (FLOAT_alt &i : temp)
        {
            matrix.push_back((FLOAT)i);
        }
    }

    vector<FLOAT> get_vec() const
    {
        return matrix;
    }

    INT get_row() const
    {
        return rows;
    }

    INT get_col() const
    {
        return cols;
    }

    FLOAT at(const INT &_row, const INT &_col) const
    {
        if (_row <= 0 || _col <= 0)
        {
            throw invalid_size();
        }
        return matrix[(_row - 1) * cols + (_col - 1)];
    }

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
};

template <typename INT, typename FLOAT>
ostream &operator<<(ostream &out, const denseM<INT, FLOAT> &_M)
{
    vector<FLOAT> temp = _M.get_vec();
    for (uint64_t i = 0; i < temp.size(); i++)
    {
        out << temp[i] << " ";
        if ((i + 1) % _M.get_col() == 0)
        {
            out << "\n";
        }
    }
    return out;
}

template <typename INT, typename FLOAT>
denseM<INT, FLOAT> operator+(const denseM<INT, FLOAT> &_M1, const denseM<INT, FLOAT> &_M2)
{
    // Exception occurs when given matrices have different size
    if (_M1.get_col() != _M2.get_col() || _M1.get_row() != _M2.get_row())
    {
        throw invalid_argument("These two matrices are not able to add each other.");
    }

    INT Row = _M1.get_row();
    INT Col = _M1.get_col();
    vector<FLOAT> M1 = _M1.get_vec();
    vector<FLOAT> M2 = _M2.get_vec();
    vector<FLOAT> sum(Row * Col, 0.0);

    for (INT i = 0; i < Row; i++)
    {
        for (INT j = 0; j < Col; j++)
        {
            sum[i * Col + j] += M1[i * Col + j] + M2[i * Col + j];
        }
    }

    denseM<INT, FLOAT> _sum(Row, Col, sum);
    return _sum;
}

template <typename INT, typename FLOAT>
denseM<INT, FLOAT> operator-(const denseM<INT, FLOAT> &_M1, const denseM<INT, FLOAT> &_M2)
{
    // Exception occurs when given matrices have different size
    if (_M1.get_col() != _M2.get_col() || _M1.get_row() != _M2.get_row())
    {
        throw invalid_argument("These two matrices are not able to subtract each other.");
    }

    INT Row = _M1.get_row();
    INT Col = _M1.get_col();
    vector<FLOAT> M1 = _M1.get_vec();
    vector<FLOAT> M2 = _M2.get_vec();
    vector<FLOAT> sub(Row * Col, 0.0);

    for (INT i = 0; i < Row; i++)
    {
        for (INT j = 0; j < Col; j++)
        {
            sub[i * Col + j] += M1[i * Col + j] - M2[i * Col + j];
        }
    }

    denseM<INT, FLOAT> _sub(Row, Col, sub);
    return _sub;
}

template <typename INT, typename FLOAT>
denseM<INT, FLOAT> operator*(const denseM<INT, FLOAT> &_M1, const denseM<INT, FLOAT> &_M2)
{
    // Exception occurs when given matrices are not able to multiply each other
    if (_M1.get_col() != _M2.get_row())
    {
        throw invalid_argument("These two matrices are not able to multiply each other.");
    }

    INT row1 = _M1.get_row();
    INT row2 = _M2.get_row();
    INT col1 = _M1.get_col();
    INT col2 = _M2.get_col();

    // size of the multiplication result is rows of matrix1 with cols of matrix2
    vector<FLOAT> result(row1 * col2);
    vector<FLOAT> M1 = _M1.get_vec();
    vector<FLOAT> M2 = _M2.get_vec();

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
    denseM<INT, FLOAT> mult(row1, col2, result);
    return mult;
}

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
void LU_withPivot(const denseM<INT, FLOAT> &_A, denseM<INT, FLOAT> &_P, denseM<INT, FLOAT> &_L, denseM<INT, FLOAT> &_U)
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

template <typename INT, typename FLOAT>
denseM<INT, FLOAT> IR(const denseM<INT, FLOAT> &_A, const denseM<INT, FLOAT> &_b)
{
    // permutation matrix from LU-decomposition
    denseM<INT, FLOAT> P;
    // lower triangular matrix from LU-decomposition
    denseM<INT, FLOAT> L;
    // upper triangular matrix from LU-decomposition
    denseM<INT, FLOAT> U;
    // Get P, L, U
    LU_withPivot(_A, P, L, U);

    denseM<INT, FLOAT> x = LU_solver(L, U, P, _b);

    INT max_iter = 1000;
    INT iter = 0;

    // residual
    denseM<INT, FLOAT> r;
    // infinity norm of r
    FLOAT residual = 1;
    // tolerance for stopping iteration
    FLOAT tol = 1e-10;

    // correction
    denseM<INT, FLOAT> c;
    while (iter != max_iter && residual > tol && residual != 0)
    {
        r = _b - (_A * x);
        residual = norm<INT, FLOAT>(r);

        // for now, using LU get correction
        c = LU_solver(L, U, P, r);
        x = x + c;
        iter++;
    }

    return x;
}