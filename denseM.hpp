/**
 * @file denseM.hpp
 * @author Yiqi Huang (huany132@mcmaster.ca)
 * @brief Construct a class that allow to store a dense matrix in any floating-point
 * number type vector, and overload several basic matrix arithmetics. The main goal is
 * using the class and the related functions, to implement LU-decomposition and Iterative
 * Refinement which supports three precisions.
 * @version 0.1
 * @date 2022-12-15
 *
 * @copyright Copyright (c) 2022
 *
 */
#include <iostream>
#include <chrono>
#include <stdexcept>
#include <vector>
#include <string>
#include <algorithm>
#include <random>
using namespace std;

/**
 * @brief Construct the class of dense matrix
 *
 * @tparam INT Any integer precision (unsigned or signed are both available)
 * @tparam FLOAT Any floating-point number precision
 */
template <typename INT, typename FLOAT>
class denseM
{
public:
    /**
     * @brief Construct a new denseM object with zeros, given the sizes
     *
     */
    denseM(const INT &, const INT &);

    /**
     * @brief Construct a new denseM object with specific elements, given the
     * sizes and any floating-point number type vector contains all the elements
     *
     */
    denseM(const INT &, const INT &, const vector<FLOAT> &);

    /**
     * @brief Construct a new denseM object accepts positive-definite matrix,
     * by given the size, the matrix itself, and a is-positive-definite checker
     *
     */
    denseM(const INT &, const vector<FLOAT> &, const bool &);

    /**
     * @brief Construct a new denseM object by given a floating-point number
     * vector. It will construct a diagonal matrix
     *
     */
    denseM(const vector<FLOAT> &);

    /**
     * @brief Construct a new denseM object with another precision
     *
     * @tparam INT_alt Another integer precision (It could remain the same)
     * @tparam FLOAT_alt Another floating-point number precision
     */
    template <typename INT_alt, typename FLOAT_alt>
    denseM(const denseM<INT_alt, FLOAT_alt> &);

    /**
     * @brief Get the vector which stores all the elements of matrix
     *
     * @return const vector<FLOAT>& matrix_ the reference of the vector stores matrix numbers
     */
    const vector<FLOAT> &get_data() const
    {
        return matrix_;
    }

    /**
     * @brief Get the number of rows
     *
     * @return INT rows_ is the number of rows of the current matrix
     */
    INT get_num_rows() const
    {
        return rows_;
    }

    /**
     * @brief Get the number of columns
     *
     * @return INT cols_ is the number of columns of the current matrix
     */
    INT get_num_cols() const
    {
        return cols_;
    }

    /**
     * @brief Get the matrix's is_pos_def information
     *
     * @return true if the matrix is positive definite
     * @return false if the matrix is not positive definite
     */
    bool get_is_pos_def() const
    {
        return is_pos_def;
    }

    /**
     * @brief Check if the matrix is symmetric
     *
     * @return true if the matrix is symmetric
     * @return false if the matrix is not symmetric
     */
    bool is_symmetric() const;

    /**
     * @brief Check the value for a given row and col
     *
     * @return const FLOAT
     */
    const FLOAT at(const INT &, const INT &) const;

    /**
     * @brief overloading operator[] for checking the value in the 1d vector
     * for a given index, allows modification
     *
     * @return FLOAT& the reference of the elements found by index
     */
    FLOAT &operator[](INT);

    /**
     * @brief overloading operator[] for checking the value in the 1d vector
     * for a given index, does not allow modification
     *
     * @return const FLOAT& the reference of the elements found by index
     */
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
    bool is_pos_def = 0;
};

/**
 * @brief Construct a new denseM object with zeros, given the sizes.
 * Preallocating the memory.
 *
 * @tparam INT Any integer precision (unsigned or signed are both available)
 * @tparam FLOAT Any floating-point number precision
 * @param rows number of rows
 * @param cols number of columns
 */
template <typename INT, typename FLOAT>
denseM<INT, FLOAT>::denseM(const INT &rows, const INT &cols)
    : rows_(rows), cols_(cols)
{
    if (rows <= 0 || cols <= 0)
    {
        throw invalid_size();
    }
    matrix_.resize(cols * rows);
}

/**
 * @brief Construct a new denseM object with specific elements, given the
 * sizes and any floating-point number type vector contains all the elements
 *
 * @tparam INT Any integer precision (unsigned or signed are both available)
 * @tparam FLOAT Any floating-point number precision
 * @param rows Number of rows
 * @param cols Number of columns
 * @param matrix A floating-point number vector contains all elements of matrix
 */
template <typename INT, typename FLOAT>
denseM<INT, FLOAT>::denseM(const INT &rows, const INT &cols, const vector<FLOAT> &matrix)
    : rows_(rows), cols_(cols), matrix_(matrix) {}

/**
 * @brief Check if the matrix is symmetric
 *
 * @tparam INT Any integer precision (unsigned or signed are both available)
 * @tparam FLOAT Any floating-point number precision
 * @return true if the matrix is symmetric
 * @return false if the matrix is not symmetric
 */
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

/**
 * @brief Construct a new denseM object with the size, the matrix vector, and a boolean
 * is-positive-definite checker. The user should know if the matrix they create is or
 * not a positive definite. Input 1 in the third argument if the user know the matrix is
 * positive definite, input 0 if it's not and it will create a regular square matrix
 *
 * @tparam INT Any integer precision (unsigned or signed are both available)
 * @tparam FLOAT Any floating-point number precision
 * @param size The size of matrix
 * @param matrix A floating-point number vector contains all elements of matrix
 * @param pos_def Boolean variable, input 1 if it is positive definite, otherwise input 0
 */
template <typename INT, typename FLOAT>
denseM<INT, FLOAT>::denseM(const INT &size, const vector<FLOAT> &matrix, const bool &pos_def)
    : rows_(size), cols_(size), matrix_(matrix)
{
    if (pos_def == 1 && this->is_symmetric() == 1)
    {
        is_pos_def = 1;
    }
}

/**
 * @brief Construct a new denseM object with another precision
 *
 * @tparam INT Any integer precision (unsigned or signed are both available)
 * @tparam FLOAT Any floating-point number precision
 * @tparam INT_alt Another integer precision (It could remain the same)
 * @tparam FLOAT_alt Another floating-point number precision
 * @param M matrix needs to change precision
 */
template <typename INT, typename FLOAT>
template <typename INT_alt, typename FLOAT_alt>
denseM<INT, FLOAT>::denseM(const denseM<INT_alt, FLOAT_alt> &M)
    : rows_((INT)M.get_num_rows()), cols_((INT)M.get_num_cols())
{
    vector<FLOAT_alt> temp = M.get_data();
    for (FLOAT_alt &i : temp)
    {
        matrix_.push_back((FLOAT)i);
    }
}

/**
 * @brief Show the value of matrix by given row number and column number.
 * It will use the mathematical index of matrix. Example M.at(3,2) will show
 * the value at 3rd row and 2nd column
 *
 * @tparam INT Any integer precision (unsigned or signed are both available)
 * @tparam FLOAT Any floating-point number precision
 * @param row Row number
 * @param col Column number
 * @return const FLOAT The value at the corresponding location in matrix
 */
template <typename INT, typename FLOAT>
const FLOAT denseM<INT, FLOAT>::at(const INT &row, const INT &col) const
{
    if (row < 0 || col < 0 || (row - 1) > rows_ || (col - 1) > cols_)
    {
        throw invalid_size();
    }
    return matrix_[(row - 1) * cols_ + (col - 1)];
}

/**
 * @brief Overload [] operator for finding specific elements given the index.
 * Allow modification
 *
 * @tparam INT Any integer precision (unsigned or signed are both available)
 * @tparam FLOAT Any floating-point number precision
 * @param index The element's index
 * @return FLOAT& The element
 */
template <typename INT, typename FLOAT>
FLOAT &denseM<INT, FLOAT>::operator[](INT index)
{
    if (index >= cols_ * rows_ || index < 0)
    {
        throw invalid_size();
    }
    return matrix_[index];
}

/**
 * @brief Overload [] operator for finding specific elements given the index.
 * Does not allow modification
 *
 * @tparam INT Any integer precision (unsigned or signed are both available)
 * @tparam FLOAT Any floating-point number precision
 * @param index The element's index
 * @return FLOAT& The element
 */
template <typename INT, typename FLOAT>
const FLOAT &denseM<INT, FLOAT>::operator[](INT index) const
{
    if (index >= cols_ * rows_ || index < 0)
    {
        throw invalid_size();
    }
    return matrix_[index];
}

/**
 * @brief Overloaded binary operator << to print the matrix
 * stored in denseM object. Each element separated by white space
 *
 * @tparam INT Any integer precision (unsigned or signed are both available)
 * @tparam FLOAT Any floating-point number precision
 * @param out An ostream object
 * @param M A denseM object
 * @return ostream& Returns a reference to an ostream object.
 */
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

/**
 * @brief Overloaded binary operator + to calculate the summation of the
 * first denseM object and the second denseM object
 *
 * @tparam INT Any integer precision (unsigned or signed are both available)
 * @tparam FLOAT Any floating-point number precision
 * @param M1 First dense matrix
 * @param M2 Second dense matrix
 * @return denseM<INT, FLOAT> A dense matrix as the result of addition.
 */
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

    // elements are adding each other correspondingly
    for (INT i = 0; i < Row * Col; i++)
    {
        sum[i] += M1[i] + M2[i];
    }
    return sum;
}

/**
 * @brief Overloaded binary operator - to calculate the subtraction of the
 * first denseM object and the second denseM object
 *
 * @tparam INT Any integer precision (unsigned or signed are both available)
 * @tparam FLOAT Any floating-point number precision
 * @param M1 First dense matrix
 * @param M2 Second dense matrix
 * @return denseM<INT, FLOAT> A dense matrix as the result of subtraction.
 */
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

    // elements are subtracting each other correspondingly
    for (INT i = 0; i < Row * Col; i++)
    {
        sub[i] += M1[i] - M2[i];
    }

    return sub;
}

/**
 * @brief Overloaded binary operator * to calculate the multiplication of the
 * first denseM object and the second denseM object
 *
 * @tparam INT Any integer precision (unsigned or signed are both available)
 * @tparam FLOAT Any floating-point number precision
 * @param M1 First dense matrix
 * @param M2 Second dense matrix
 * @return denseM<INT, FLOAT> A dense matrix as the result of multiplication.
 */
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

    for (INT i = 0; i < row1; i++)
    {
        for (INT j = 0; j < col2; j++)
        {
            for (INT k = 0; k < row2; k++)
            {
                result[col2 * i + j] += M1[col1 * i + k] * M2[col2 * k + j];
            }
        }
    }

    return result;
}

/**
 * @brief infinity norm of a denseM
 *
 * @tparam INT Any integer precision (unsigned or signed are both available)
 * @tparam FLOAT Any floating-point number precision
 * @param M Dense matrix
 * @return FLOAT The infinity norm of M
 */
template <typename INT, typename FLOAT>
FLOAT norm(const denseM<INT, FLOAT> &M)
{
    INT Row = M.get_num_rows();
    INT Col = M.get_num_cols();
    FLOAT maxRowSum = 0.0;

    // find the max row sum
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

/**
 * @brief LU-decomposition with partial pivoting, it will change the denseM object A to the
 * combination of L and U for saving memory, and modify the integer vector P to record row
 * swapping in permutation matrix
 *
 * @tparam INT Any integer precision (unsigned or signed are both available)
 * @tparam FLOAT Any floating-point number precision
 * @param A Dense matrix, it will be modified to the combination of L and U
 * @param P A integer vector P = {0,1,2,...,n-1} for nxn matrix A
 * @return INT exit-flag of LU-decomposition
 * (return 0 means success, return >0 means completed but U is singular)
 */
template <typename INT, typename FLOAT>
INT LU_withPivot(denseM<INT, FLOAT> &A, vector<INT> &P)
{
    // Since we want to solve linear system, only consider LU decomposition in square matrices
    // Exception occurs when given matrix is not a square matrix
    if (A.get_num_cols() != A.get_num_rows())
    {
        throw invalid_argument("This matrix is not a square matrix.");
    }

    INT size = A.get_num_rows();

    sort(P.begin(), P.end());

    for (INT i = 0; i < size; i++)
    {
        INT current_pivot = i;
        // INT max = 0;

        // compare the pivot in the current row with the number in
        // the following rows
        for (INT j = i + 1; j < size; j++)
        {
            // find the biggest absolute value for each column, choose it
            // as the new current pivot
            if (fabs(A[j * size + i]) > fabs(A[current_pivot * size + i]))
            {
                current_pivot = j;

                // max = fabs(A[j * size + i]);
            }
        }
        /*
        if (max == 0){
            throw invalid_argument("There is a column only contains 0, cannot compute LU")
        }*/

        // If pivot changed, swap rows
        if (current_pivot != i)
        {
            for (INT n = 0; n < size; n++)
            {
                FLOAT temp = A[current_pivot * size + n];
                A[current_pivot * size + n] = A[i * size + n];
                A[i * size + n] = temp;
            }
            swap(P[i], P[current_pivot]);
        }

        // By Gaussian elimination, A will be modified to the combination of L and U
        for (INT j = i + 1; j < size; j++)
        {
            // Preventing division by 0
            if (A[i * size + i] != 0)
            {
                A[j * size + i] /= A[i * size + i];
            }
            for (INT k = i + 1; k < size; k++)
            {
                A[j * size + k] -= A[i * size + k] * A[j * size + i];
            }
        }
    }

    // Count 0s on diagonal, if there is, means U is singular
    INT count = size;
    for (INT i = 0; i < size; i++)
    {
        if (A[i * size + i] == 0)
        {
            count = i;
        }
    }

    // print out messages with exit flag
    // if LU_withPivot = i > 0, U(i,i) is exactly zero
    if (count != size)
    {
        cout << "LU-decomposition is complete, the upper-triangular matrix's diagonal number in row"
             << count + 1 << " is 0, it is singular."
             << "\n";
        return count;
    }

    // successful exit, LU_withPivot = 0
    // cout << "LU-decomposition is success" << "\n";
    return 0;
}

/**
 * @brief Use the updated A and P from LU_withPivot(A,P), solve linear system Ax = y
 * by LU-decomposition
 *
 * @tparam INT Any integer precision (unsigned or signed are both available)
 * @tparam FLOAT Any floating-point number precision
 * @param A nxn dense matrix
 * @param b nx1 dense matrix
 * @param P A integer vector P = {0,1,2,...,n-1} for nxn matrix A
 * @return denseM<INT, FLOAT> nx1 dense matrix, it is the solution of the linear system
 */
template <typename INT, typename FLOAT>
denseM<INT, FLOAT> LU_solver(denseM<INT, FLOAT> &A, const denseM<INT, FLOAT> &b, vector<INT> &P)
{
    // Receive updated A and P
    // A contains L and U key parameters
    // P contains the row swapping information
    INT exit_flag = LU_withPivot(A, P);
    if (exit_flag != 0)
    {
        throw invalid_argument("This matrix cannot solve by LU");
    }

    INT size = A.get_num_cols();

    // Starts with LUx = Pb
    // Calculate Pb
    denseM<INT, FLOAT> Pb(size, 1);
    for (INT i = 0; i < size; i++)
    {
        Pb[i] = b[(P[i])];
    }

    // Let y = Ux, then Ly = Pb, solve y by forward substituion
    denseM<INT, FLOAT> y(size, 1);
    for (INT i = 0; i < size; i++)
    {
        FLOAT temp = 0;
        // temp will store the summation of all known xi value times corresponding
        // L value in each row
        for (INT j = 0; j <= i; j++)
        {
            // The number on the diagonal should be 1
            if (j == i)
            {
                temp += y[j] * 1;
            }
            // Lower-triangular matrix elements other than the diagonal numbers
            else if (j < i)
            {
                temp += y[j] * A[i * size + j];
            }
        }
        // b subtracts temp will left out the value of the current yi times
        // corresponding L value, divide this L number which is on diagonal,
        // we know it is 1. Then we can get yi
        y[i] = (Pb[i] - temp) / 1;
    }

    // Solve x from Ux = y by backward substituion
    denseM<INT, FLOAT> x(size, 1);
    for (INT i = size; i >= 1; i--)
    {
        FLOAT temp = 0;
        // temp will store the summation of all known xi value times corresponding
        // U value in each row
        for (INT j = i - 1; j < size; j++)
        {
            temp += x[j] * A[(i - 1) * size + j];
        }
        // y subtracts temp will left out the value of the current xi times
        // corresponding U value, divide this U number, we can get xi
        x[i - 1] = (y[i - 1] - temp) / A[(i - 1) * size + (i - 1)];
    }
    return x;
}

/**
 * @brief
 *
 * @tparam INT
 * @tparam FLOAT
 * @param A
 * @return denseM<INT, FLOAT>
 */
template <typename INT, typename FLOAT>
denseM<INT, FLOAT> cholesky(const denseM<INT, FLOAT> &A)
{
    INT size = A.get_num_cols();
    // INT exit_flag = 0;
    denseM<INT, FLOAT> A_chole(size, size);
    for (INT i = 0; i < size; i++)
    {
        for (INT j = 0; j < (i + 1); j++)
        {
            FLOAT sum = 0;
            for (INT k = 0; k < j; k++)
            {
                sum += A_chole[i * size + k] * A_chole[j * size + k];
            }
            // on diagonal
            if (i == j)
            {
                /*if (A[i * size + i] < sum)
                {
                    exit_flag = i + 1;
                }*/
                A_chole[i * size + j] = sqrt(A[i * size + i] - sum);
            }
            // not on diagonal
            else
            {
                A_chole[i * size + j] = (1.0 / A_chole[j * size + j] * (A[i * size + j] - sum));
                // It is symmetric, A_chole[j,i] = A_chole[i,j]
                A_chole[j * size + i] = A_chole[i * size + j];
            }
        }
    }

    /*
    if (exit_flag != 0){
       cout << "Cholesky-decomposition is complete, the L's diagonal number in row"
             << exit_flag << " is a complex number."
             << "\n";
        return exit_flag;
    }
    return exit_flag;*/
    return A_chole;
}

template <typename INT, typename FLOAT>
denseM<INT, FLOAT> cholesky_solver(const denseM<INT, FLOAT> &A, const denseM<INT, FLOAT> &b)
{
    // Get A_chole (which contains both the lower triangular L and upper
    // triangular L's conjugate transpose L_T) from cholesky decomposition
    denseM<INT, FLOAT> A_chole = cholesky(A);

    INT size = A.get_num_cols();

    // Start with L(L_T)x = b;
    // Let y = (L_T)x, then Ly = b, solve y by forward substituion
    denseM<INT, FLOAT> y(size, 1);
    for (INT i = 0; i < size; i++)
    {
        FLOAT temp = 0;
        // temp will store the summation of all known xi value times corresponding
        // L value in each row
        for (INT j = 0; j <= i; j++)
        {

            temp += y[j] * A_chole[i * size + j];
        }
        // b subtracts temp will left out the value of the current yi times
        // corresponding L value. Divide this L number, which is the number on diagonal
        y[i] = (b[i] - temp) / A_chole[i * size + i];
    }

    // Solve x from (L_T)x = y by backward substituion
    denseM<INT, FLOAT> x(size, 1);
    for (INT i = size; i >= 1; i--)
    {
        FLOAT temp = 0;
        // temp will store the summation of all known xi value times corresponding
        // U value in each row
        for (INT j = i - 1; j < size; j++)
        {
            temp += x[j] * A_chole[(i - 1) * size + j];
        }
        // y subtracts temp will left out the value of the current xi times
        // corresponding U value, divide this U number, we can get xi
        x[i - 1] = (y[i - 1] - temp) / A_chole[(i - 1) * size + (i - 1)];
    }
    return x;
}
/**
 * @brief
 *
 * @tparam INT Any integer precision (unsigned or signed are both available)
 * @tparam FACT floating-point number precision used for LU-factorization
 * (need to be a low-accurate precision)
 * @tparam FLOAT Any floating-point number precision
 * (the accuracy need to be between FACT and RES)
 * @tparam RES floating-point number precision used for calculating residual
 * (need to be a accurate precision)
 * @param A nxn dense matrix in precision FLOAT
 * @param b nx1 dense matrix in precision FLOAT
 * @return denseM<INT, FLOAT> nx1 dense matrix, it is the solution of the linear system
 * solved by iterative refinement
 */
template <typename INT, typename FACT, typename FLOAT, typename RES>
denseM<INT, FLOAT> IR(denseM<INT, FLOAT> &A, const denseM<INT, FLOAT> &b)
{
    INT size = A.get_num_cols();

    // change to FACT precision for factorization
    denseM<INT, FACT> Af(A);
    denseM<INT, FACT> bf(b);

    // Creating original permutation matrix in vector form
    // P = {0,1,2,...,n-1} for nxn matrix A
    vector<INT> P(size);
    for (INT i = 0; i < size; i++)
    {
        P[i] = i;
    }

    // max iteration
    INT max_iter = 1000;
    // iteration counter
    INT iter = 0;

    // residual
    denseM<INT, RES> r(size, size);
    // infinity norm of r
    RES residual = 1;
    // tolerance for stopping iteration
    FLOAT tol = 1e-16;

    // correction
    denseM<INT, FLOAT> c(size, 1);

    // timing the process
    chrono::time_point start_time = chrono::steady_clock::now();

    // Calculate x0 for starting the iteration
    // L,U in Af must be in FACT precision
    // LU_solver returns x in FACT precision
    denseM<INT, FLOAT> x = LU_solver(Af, bf, P);
    Af = denseM<INT, FACT>(A);
    while (iter != max_iter && residual > tol && residual != 0)
    {
        // residual must be in RES precision
        r = b - (A * denseM<INT, RES>(x));
        residual = norm<INT, RES>(r);

        // Using LU again in FACT precision to get correction
        c = LU_solver(Af, denseM<INT, FACT>(r), P);
        Af = denseM<INT, FACT>(A);
        x = x + c;
        iter++;
    }
    chrono::time_point end_time = chrono::steady_clock::now();
    chrono::duration<double> elapsed_time_seconds = end_time - start_time;
    cout << "Elapsed time: " << elapsed_time_seconds.count() << " seconds\n";
    cout << "The total iteration is " << iter << "\n";
    cout << "The error in the last iteration is " << residual << "\n";
    return x;
}
