/**
 * @file denseM.hpp
 * @author Yiqi Huang (huany132@mcmaster.ca)
 * @brief Construct a class that allow to store a dense matrix in any floating-point
 * number type vector, and overload several basic matrix arithmetics. The main goal is
 * using the class and the related functions, to implement LU-decomposition/Cholesky
 * decomposition and Iterative Refinement which uses three floating-point number precisions.
 * @version 0.1
 * @date 2022-12-14
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
#include <cmath>
#include <fstream>
using namespace std;

/**
 * @brief Construct the class of dense matrix
 *
 * @tparam FLOAT Any floating-point number precision
 */
template <typename FLOAT>
class denseM
{
public:
    /**
     * @brief Construct a new denseM object with zeros, given the sizes
     *
     */
    denseM(const uint64_t &, const uint64_t &);

    /**
     * @brief Construct a new denseM object with specific elements, given the
     * sizes and any floating-point number type vector contains all the elements
     *
     */
    denseM(const uint64_t &, const uint64_t &, const vector<FLOAT> &);

    /**
     * @brief Construct a new denseM object accepts positive-definite matrix,
     * by given the size, the matrix itself, and a is-positive-definite checker
     *
     */
    denseM(const uint64_t &, const vector<FLOAT> &, const bool &);

    /**
     * @brief Construct a new denseM object by given a integer
     *  It will construct a nxn identity matrix
     *
     */
    denseM(const uint64_t &);

    /**
     * @brief Construct a new denseM object given the size and file name
     *
     */
    denseM(const uint64_t &, const uint64_t &, const string &);

    /**
     * @brief Construct a new denseM object with another precision
     *
     * @tparam INT_alt Another integer precision (It could remain the same)
     * @tparam FLOAT_alt Another floating-point number precision
     */
    template <typename FLOAT_alt>
    denseM(const denseM<FLOAT_alt> &);

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
     * @return uint64_t Rows_ is the number of rows of the current matrix
     */
    uint64_t get_num_rows() const
    {
        return rows_;
    }

    /**
     * @brief Get the number of columns
     *
     * @return uint64_t Cols_ is the number of columns of the current matrix
     */
    uint64_t get_num_cols() const
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
    const FLOAT at(const uint64_t &, const uint64_t &) const;

    /**
     * @brief overloading operator[] for checking the value in the 1d vector
     * for a given index, allows modification
     *
     * @return FLOAT& the reference of the elements found by index
     */
    FLOAT &operator[](uint64_t);

    /**
     * @brief overloading operator[] for checking the value in the 1d vector
     * for a given index, does not allow modification
     *
     * @return const FLOAT& the reference of the elements found by index
     */
    const FLOAT &operator[](uint64_t) const;

    /**
     * @brief Output denseM object into a file
     *
     */
    void output(const string &) const;

    /**
     * @brief Exception occurs when given matrix vector has invalid number of rows or columns
     *
     */
    class invalid_size : public invalid_argument
    {
    public:
        invalid_size() : invalid_argument("The matrix size can only be positive."){};
    };

    /**
     * @brief Exception occurs when given matrix vector has invalid number of rows or columns
     *
     */
    class index_overflow : public invalid_argument
    {
    public:
        index_overflow() : invalid_argument("The input index is beyond the size of the matrix."){};
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
    uint64_t rows_ = 0;
    uint64_t cols_ = 0;
    vector<FLOAT> matrix_;
    bool is_pos_def = 0;
};

/**
 * @brief Construct a new denseM object with zeros in the matrix, given the sizes.
 * It is for preallocating the memory.
 *
 * @tparam FLOAT Any floating-point number precision
 * @param rows number of rows
 * @param cols number of columns
 */
template <typename FLOAT>
denseM<FLOAT>::denseM(const uint64_t &rows, const uint64_t &cols)
    : rows_(rows), cols_(cols)
{
    if (rows == 0 || cols == 0)
    {
        throw invalid_size();
    }
    matrix_.resize(cols * rows);
}

/**
 * @brief Construct a new denseM object with specific elements, given the
 * sizes and a floating-point number type vector contains all the elements
 * of the matrix
 *
 * @tparam FLOAT Any floating-point number precision
 * @param rows Number of rows
 * @param cols Number of columns
 * @param matrix A floating-point number vector contains all elements of matrix
 */
template <typename FLOAT>
denseM<FLOAT>::denseM(const uint64_t &rows, const uint64_t &cols, const vector<FLOAT> &matrix)
    : rows_(rows), cols_(cols), matrix_(matrix)
{
    if (rows == 0 || cols == 0)
    {
        throw invalid_size();
    }
    else if (rows * cols != matrix.size())
    {
        throw size_mismatch();
    }
}

/**
 * @brief Check if the matrix is symmetric
 *
 * @tparam FLOAT Any floating-point number precision
 * @return true if the matrix is symmetric
 * @return false if the matrix is not symmetric
 */
template <typename FLOAT>
bool denseM<FLOAT>::is_symmetric() const
{
    // Matrix must be a square matrix
    if (rows_ != cols_)
    {
        return 0;
    }
    for (uint64_t i = 0; i < rows_; i++)
    {
        for (uint64_t j = 0; j < cols_; j++)
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
 * not positive definite. Input 1 in the third argument if the user know the matrix is
 * positive definite, input 0 if it's not and it will create a regular square matrix
 *
 * @tparam FLOAT Any floating-point number precision
 * @param size The size of matrix
 * @param matrix A floating-point number vector contains all elements of matrix
 * @param pos_def Boolean variable, input 1 if it is positive definite, otherwise input 0
 */
template <typename FLOAT>
denseM<FLOAT>::denseM(const uint64_t &size, const vector<FLOAT> &matrix, const bool &pos_def)
    : rows_(size), cols_(size), matrix_(matrix)
{
    if (pos_def == 1 && this->is_symmetric() == 1)
    {
        is_pos_def = 1;
    }

    if (size == 0)
    {
        throw invalid_size();
    }
    else if (size * size != matrix.size())
    {
        throw size_mismatch();
    }
}

/**
 * @brief Construct a new denseM with an integer, create a nxn identity matrix
 *
 * @tparam FLOAT Any floating-point number precision
 * @param size Size of the identity matrix
 */
template <typename FLOAT>
denseM<FLOAT>::denseM(const uint64_t &size)
    : rows_(size), cols_(size)
{
    matrix_.resize(size * size);
    for (uint64_t i = 0; i < size; i++)
    {
        matrix_[i * size + i] = 1;
    }
}

/**
 * @brief Construct a new denseM object given size and a file name which contains matrix.
 * The matrix number in the file can be separated by comma or white space.
 *
 * @tparam FLOAT Any floating-point number precision
 * @param rows Number of rows
 * @param cols Number of columns
 * @param name File's name
 */
template <typename FLOAT>
denseM<FLOAT>::denseM(const uint64_t &rows, const uint64_t &cols, const string &name)
    : rows_(rows), cols_(cols)
{
    if (rows == 0 || cols == 0)
    {
        throw invalid_size();
    }
    ifstream input(name);
    if (!input.is_open())
    {
        cout << "There is a problem when opening the file " << name;
        exit(EXIT_FAILURE);
    }

    // Read floating-point number from file into variable
    FLOAT value;
    while (input >> value)
    {
        // Push the value into the vector
        matrix_.push_back(value);
        // ignore commas
        if (input.peek() == ',')
        {
            input.ignore();
        }
    }
    input.close();

    if (rows * cols != matrix_.size())
    {
        throw size_mismatch();
    }
}

/**
 * @brief Construct a new denseM object based on a existed denseM
 * object with another precision.
 *
 * @tparam FLOAT Any floating-point number precision
 * @tparam FLOAT_alt Another floating-point number precision
 * @param M matrix needs to change precision
 */
template <typename FLOAT>
template <typename FLOAT_alt>
denseM<FLOAT>::denseM(const denseM<FLOAT_alt> &M)
    : rows_(M.get_num_rows()), cols_(M.get_num_cols())
{
    vector<FLOAT_alt> temp = M.get_data();
    for (FLOAT_alt &i : temp)
    {
        matrix_.push_back((FLOAT)i);
    }
}

/**
 * @brief Show the value of matrix by given row number and column number.
 * It will use the mathematical index of matrix. For example, M.at(3,2) will show
 * the value at 3rd row and 2nd column in matrix M
 *
 * @tparam FLOAT Any floating-point number precision
 * @param row Row number
 * @param col Column number
 * @return const FLOAT The value at the corresponding location in matrix
 */
template <typename FLOAT>
const FLOAT denseM<FLOAT>::at(const uint64_t &row, const uint64_t &col) const
{
    if (row == 0 || col == 0)
    {
        throw invalid_size();
    }
    else if ((row - 1) > rows_ || (col - 1) > cols_)
    {
        throw index_overflow();
    }

    return matrix_[(row - 1) * cols_ + (col - 1)];
}

/**
 * @brief Overload [] operator for finding specific elements given the index.
 * Allow modification
 *
 * @tparam FLOAT Any floating-point number precision
 * @param index The element's index
 * @return FLOAT& The element
 */
template <typename FLOAT>
FLOAT &denseM<FLOAT>::operator[](uint64_t index)
{
    if (index >= cols_ * rows_)
    {
        throw index_overflow();
    }
    return matrix_[index];
}

/**
 * @brief Overload [] operator for finding specific elements given the index.
 * Does not allow modification
 *
 * @tparam FLOAT Any floating-point number precision
 * @param index The element's index
 * @return FLOAT& The element
 */
template <typename FLOAT>
const FLOAT &denseM<FLOAT>::operator[](uint64_t index) const
{
    if (index >= cols_ * rows_)
    {
        throw index_overflow();
    }
    return matrix_[index];
}

template <typename FLOAT>
void denseM<FLOAT>::output(const string &filename) const
{
    ofstream out(filename);
    if (!out.is_open())
    {
        cout << "There is a problem when opening the file " << filename;
        exit(EXIT_FAILURE);
    }

    // Find the extension of the filename
    size_t dot_pos = filename.find_last_of(".");
    string extension = filename.substr(dot_pos);
    for (uint64_t i = 0; i < rows_; i++)
    {
        for (uint64_t j = 0; j < cols_; j++)
        {
            // If this is an excel file, use comma to separate numbers
            if (extension == ".csv")
            {
                out << matrix_[i * cols_ + j];
                if (j + 1 != cols_)
                {
                    out << ",";
                }
            }
            else
            {
                out << matrix_[i * cols_ + j] << " ";
            }
        }
        out << "\n";
    }

    out.close();
}

/**
 * @brief Overloaded binary operator << to print the matrix
 * stored in denseM object. Each element separated by white space
 *
 * @tparam FLOAT Any floating-point number precision
 * @param out An ostream object
 * @param M A denseM object
 * @return ostream& Returns a reference to an ostream object.
 */
template <typename FLOAT>
ostream &operator<<(ostream &out, const denseM<FLOAT> &M)
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
 * @tparam FLOAT Any floating-point number precision
 * @param M1 First dense matrix
 * @param M2 Second dense matrix
 * @return denseM<FLOAT> A dense matrix as the result of addition.
 */
template <typename FLOAT>
denseM<FLOAT> operator+(const denseM<FLOAT> &M1, const denseM<FLOAT> &M2)
{
    // Exception occurs when given matrices have different size
    if (M1.get_num_cols() != M2.get_num_cols() || M1.get_num_rows() != M2.get_num_rows())
    {
        throw invalid_argument("These two matrices are not able to add each other.");
    }

    uint64_t Row = M1.get_num_rows();
    uint64_t Col = M1.get_num_cols();
    denseM<FLOAT> sum(Row, Col);

    // elements are adding each other correspondingly
    for (uint64_t i = 0; i < Row * Col; i++)
    {
        sum[i] += M1[i] + M2[i];
    }
    return sum;
}

/**
 * @brief Overloaded binary operator - to calculate the subtraction of the
 * first denseM object and the second denseM object
 *
 * @tparam FLOAT Any floating-point number precision
 * @param M1 First dense matrix
 * @param M2 Second dense matrix
 * @return denseM<FLOAT> A dense matrix as the result of subtraction.
 */
template <typename FLOAT>
denseM<FLOAT> operator-(const denseM<FLOAT> &M1, const denseM<FLOAT> &M2)
{
    // Exception occurs when given matrices have different size
    if (M1.get_num_cols() != M2.get_num_cols() || M1.get_num_rows() != M2.get_num_rows())
    {
        throw invalid_argument("These two matrices are not able to subtract each other.");
    }

    uint64_t Row = M1.get_num_rows();
    uint64_t Col = M1.get_num_cols();

    denseM<FLOAT> sub(Row, Col);

    // elements are subtracting each other correspondingly
    for (uint64_t i = 0; i < Row * Col; i++)
    {
        sub[i] += M1[i] - M2[i];
    }

    return sub;
}

/**
 * @brief Overloaded binary operator * to calculate the multiplication of the
 * first denseM object and the second denseM object
 *
 * @tparam FLOAT Any floating-point number precision
 * @param M1 First dense matrix
 * @param M2 Second dense matrix
 * @return denseM<FLOAT> A dense matrix as the result of multiplication.
 */
template <typename FLOAT>
denseM<FLOAT> operator*(const denseM<FLOAT> &M1, const denseM<FLOAT> &M2)
{
    // Exception occurs when given matrices are not able to multiply each other
    if (M1.get_num_cols() != M2.get_num_rows())
    {
        throw invalid_argument("These two matrices are not able to multiply each other.");
    }

    uint64_t row1 = M1.get_num_rows();
    uint64_t row2 = M2.get_num_rows();
    uint64_t col1 = M1.get_num_cols();
    uint64_t col2 = M2.get_num_cols();

    // size of the multiplication result is rows of matrix1 with cols of matrix2
    denseM<FLOAT> result(row1, col2);

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

/**
 * @brief infinity norm of a denseM. Finding the largest row sum of the matrix.
 *
 * @tparam FLOAT Any floating-point number precision
 * @param M Dense matrix
 * @return FLOAT The infinity norm of M
 */
template <typename FLOAT>
FLOAT norm(const denseM<FLOAT> &M)
{
    uint64_t Row = M.get_num_rows();
    uint64_t Col = M.get_num_cols();
    FLOAT maxRowSum = 0.0;

    // find the max row sum
    for (uint64_t i = 0; i < Row; i++)
    {
        FLOAT RowSum = 0.0;
        for (uint64_t j = 0; j < Col; j++)
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
 * @tparam FLOAT Any floating-point number precision
 * @param A Dense matrix, it will be modified to the combination of L and U
 * @param P A integer vector P = {0,1,2,...,n-1} for nxn matrix A
 * @return INT exit-flag of LU-decomposition
 * (return 0 means success, return >0 means completed but U is singular)
 */
template <typename FLOAT>
uint64_t LU_withPivot(denseM<FLOAT> &A, vector<uint64_t> &P)
{
    // Since we want to solve linear system, only consider LU decomposition in square matrices
    // Exception occurs when given matrix is not a square matrix
    if (A.get_num_cols() != A.get_num_rows())
    {
        throw invalid_argument("This matrix is not a square matrix.");
    }

    uint64_t size = A.get_num_rows();

    sort(P.begin(), P.end());

    for (uint64_t i = 0; i < size; i++)
    {
        uint64_t current_pivot = i;

        // compare the pivot in the current row with the number in
        // the following rows
        for (uint64_t j = i + 1; j < size; j++)
        {
            // find the biggest absolute value for each column, choose it
            // as the new current pivot
            if (fabs(A[j * size + i]) > fabs(A[current_pivot * size + i]))
            {
                current_pivot = j;

                // max = fabs(A[j * size + i]);
            }
        }

        // If pivot changed, swap rows
        if (current_pivot != i)
        {
            for (uint64_t n = 0; n < size; n++)
            {
                FLOAT temp = A[current_pivot * size + n];
                A[current_pivot * size + n] = A[i * size + n];
                A[i * size + n] = temp;
            }
            swap(P[i], P[current_pivot]);
        }

        // By Gaussian elimination, A will be modified to the combination of L and U
        for (uint64_t j = i + 1; j < size; j++)
        {
            // Preventing division by 0
            if (A[i * size + i] != 0)
            {
                A[j * size + i] /= A[i * size + i];
            }
            for (uint64_t k = i + 1; k < size; k++)
            {
                A[j * size + k] -= A[i * size + k] * A[j * size + i];
            }
        }
    }

    // Count 0s on diagonal, if there is, means U is singular
    uint64_t count = size;
    for (uint64_t i = 0; i < size; i++)
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

template <typename FLOAT>
denseM<FLOAT> Lower_triangular(const denseM<FLOAT> &A)
{
    uint64_t size = A.get_num_cols();

    // Separate lower-triangular matrix out
    // Create an identity matrix since lower-triangular matrix has diagonal with 1s
    denseM<FLOAT> lower(size, size);
    for (uint64_t i = 0; i < size; i++)
    {
        lower[size * i + i] = 1;
    }

    for (uint64_t i = 0; i < size; i++)
    {

        for (uint64_t j = 0; j < size; j++)
        {
            if (j < i)
            {
                lower[i * size + j] = A[i * size + j];
            }
        }
    }
    return lower;
}

template <typename FLOAT>
denseM<FLOAT> Upper_triangular(const denseM<FLOAT> &A)
{

    uint64_t size = A.get_num_cols();

    // Separate upper-triangular matrix out
    // The upper part of A is just the upper-triangular matrix U
    denseM<FLOAT> upper(size, size);
    for (uint64_t i = 0; i < size; i++)
    {
        for (uint64_t zero = 0; zero < i; zero++)
        {
            upper.push_back(0);
        }
        for (uint64_t j = i; j < size; j++)
        {
            upper.push_back(A[i * size + j]);
        }
    }

    return upper;
}

/**
 * @brief Use the updated A and P from LU_withPivot(A,P), solve linear system LUx = Pb
 * by forward and backward substitution
 *
 * @tparam FLOAT Any floating-point number precision
 * @param A nxn dense matrix
 * @param b nx1 dense matrix
 * @param P A integer vector P = {0,1,2,...,n-1} for nxn matrix A
 * @return denseM<FLOAT> nx1 dense matrix, it is the solution of the linear system
 */
template <typename FLOAT>
denseM<FLOAT> LU_solver(denseM<FLOAT> &A, const denseM<FLOAT> &b, vector<uint64_t> &P)
{
    if (A.get_num_rows() != b.get_num_rows())
    {
        throw invalid_argument("b should have the same column number as A.");
    }

    // Receive updated A and P
    // A contains L and U key parameters
    // P contains the row swapping information
    uint64_t exit_flag = LU_withPivot(A, P);
    if (exit_flag != 0)
    {
        throw invalid_argument("This matrix cannot solve by LU");
    }

    uint64_t size = A.get_num_cols();

    // Starts with LUx = Pb
    // Calculate Pb
    denseM<FLOAT> Pb(size, 1);
    for (uint64_t i = 0; i < size; i++)
    {
        Pb[i] = b[(P[i])];
    }

    // Let y = Ux, then Ly = Pb, solve y by forward substituion
    denseM<FLOAT> y(size, 1);
    for (uint64_t i = 0; i < size; i++)
    {
        FLOAT temp = 0;
        // temp will store the summation of all known xi value times corresponding
        // L value in each row
        for (uint64_t j = 0; j <= i; j++)
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
    denseM<FLOAT> x(size, 1);
    for (uint64_t i = size; i >= 1; i--)
    {
        FLOAT temp = 0;
        // temp will store the summation of all known xi value times corresponding
        // U value in each row
        for (uint64_t j = i - 1; j < size; j++)
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
 * @brief Cholesky decomposition will decompose a positive definite matrix A into a
 * lower-triangular matrix L and its conjugate transpose L_T, such that L*L_T = A
 *
 * @tparam FLOAT Any floating-point number precision
 * @param A Dense matrix, it will not be modified
 * @return denseM<FLOAT> A_chole, contains the combination of L and L_T
 */
template <typename FLOAT>
denseM<FLOAT> cholesky(const denseM<FLOAT> &A)
{
    uint64_t size = A.get_num_cols();
    denseM<FLOAT> A_chole(size, size);

    // Modified from
    // https://en.wikipedia.org/wiki/Cholesky_decomposition#The_Cholesky_algorithm
    for (uint64_t i = 0; i < size; i++)
    {
        for (uint64_t j = 0; j < (i + 1); j++)
        {
            FLOAT sum = 0;
            for (uint64_t k = 0; k < j; k++)
            {
                sum += A_chole[i * size + k] * A_chole[j * size + k];
            }
            // on diagonal
            if (i == j)
            {
                A_chole[i * size + j] = sqrt(A[i * size + i] - sum);
            }
            // not on diagonal
            else
            {
                A_chole[i * size + j] = ((FLOAT)(1.0) / A_chole[j * size + j] * (A[i * size + j] - sum));
                // It is symmetric, A_chole[j,i] = A_chole[i,j]
                A_chole[j * size + i] = A_chole[i * size + j];
            }
        }
    }
    // Modification ends here

    return A_chole;
}

/**
 * @brief Get the combination matrix A_chole of L and L_T from cholesky(A,P),
 * solve linear system LL_Tx=b by forward and backward substitution
 *
 * @tparam FLOAT Any floating-point number precision
 * @param A nxn dense matrix
 * @param b nx1 dense matrix
 * @return denseM<FLOAT> nx1 dense matrix, it is the solution of the linear system
 */
template <typename FLOAT>
denseM<FLOAT> cholesky_solver(const denseM<FLOAT> &A, const denseM<FLOAT> &b)
{
    if (A.get_num_rows() != b.get_num_rows())
    {
        throw invalid_argument("b should have the same column number as A.");
    }

    // Get A_chole (which contains both the lower triangular L and upper
    // triangular L's conjugate transpose L_T) from cholesky decomposition
    denseM<FLOAT> A_chole = cholesky(A);

    uint64_t size = A.get_num_cols();

    // Start with L(L_T)x = b;
    // Let y = (L_T)x, then Ly = b, solve y by forward substituion
    denseM<FLOAT> y(size, 1);
    for (uint64_t i = 0; i < size; i++)
    {
        FLOAT temp = 0;
        // temp will store the summation of all known xi value times corresponding
        // L value in each row
        for (uint64_t j = 0; j <= i; j++)
        {

            temp += y[j] * A_chole[i * size + j];
        }
        // b subtracts temp will left out the value of the current yi times
        // corresponding L value. Divide this L number, which is the number on diagonal
        y[i] = (b[i] - temp) / A_chole[i * size + i];
    }

    // Solve x from (L_T)x = y by backward substituion
    denseM<FLOAT> x(size, 1);
    for (uint64_t i = size; i >= 1; i--)
    {
        FLOAT temp = 0;
        // temp will store the summation of all known xi value times corresponding
        // U value in each row
        for (uint64_t j = i - 1; j < size; j++)
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
 * @brief Using iterative refinement to solve linear system for less round off error.
 * Including LU-decomposition and Cholesky decomposition to solve linear system inside
 * iterations. Using three floating-point number precisions to accelerate decomposition
 * by low-accuracy precision, and getting a more precise residual for updating the
 * solution in each iteration by using a high-accuracy precision
 *
 * @tparam FACT floating-point number precision used for LU-factorization
 * (need to be a low-accurate precision)
 * @tparam FLOAT Any floating-point number precision
 * (the accuracy need to be between FACT and RES)
 * @tparam RES floating-point number precision used for calculating residual
 * (need to be a accurate precision)
 * @param A nxn dense matrix in precision FLOAT
 * @param b nx1 dense matrix in precision FLOAT
 * @return denseM<FLOAT> nx1 dense matrix, it is the solution of the linear system
 * solved by iterative refinement
 */
template <typename FACT, typename FLOAT, typename RES>
denseM<FLOAT> IR(denseM<FLOAT> &A, const denseM<FLOAT> &b)
{
    uint64_t size = A.get_num_cols();

    // change to FACT precision for factorization
    denseM<FACT> Af(A);
    denseM<FACT> bf(b);

    // Creating original permutation matrix in vector form
    // P = {0,1,2,...,n-1} for nxn matrix A
    vector<uint64_t> P(size);
    for (uint64_t i = 0; i < size; i++)
    {
        P[i] = i;
    }

    // max iteration
    uint64_t max_iter = 10000;
    // iteration counter
    uint64_t iter = 0;

    // residual
    denseM<RES> r(size, size);
    // infinity norm of r
    RES residual = 1;
    // tolerance for stopping iteration
    FLOAT tol = 1e-16;
    FLOAT tol2 = 1e-14;
    // correction
    denseM<FLOAT> c(size, 1);

    cout << "Starting iterative refinement: "
         << "\n";
    // timing the process
    chrono::time_point start_time = chrono::steady_clock::now();

    // If A is a positive definite matrix, use Cholesky-decomposition
    if (A.get_is_pos_def() == 1)
    {
        denseM<FLOAT> x = cholesky_solver(Af, bf);
        while (iter != max_iter && residual > tol2 && residual != 0)
        {
            // residual must be in RES precision
            r = b - (A * denseM<RES>(x));
            residual = norm<RES>(r);

            // Using Cholesky again in FACT precision to get correction
            c = cholesky_solver(Af, denseM<FACT>(r));
            x = x + c;
            iter++;
        }

        chrono::time_point end_time = chrono::steady_clock::now();
        chrono::duration<double> elapsed_time_seconds = end_time - start_time;
        cout << "Elapsed time: " << elapsed_time_seconds.count() << " seconds\n";
        cout << "The total iteration is " << iter << "\n";
        cout << "The error in the last iteration is " << residual << "\n";
        cout << "Iterative refinement succeeded!"
             << "\n"
             << "x = ";
        return x;
    }

    // Calculate x0 for starting the iteration
    // L,U in Af must be in FACT precision
    // LU_solver returns x in FACT precision
    denseM<FLOAT> x = LU_solver(Af, bf, P);
    Af = denseM<FACT>(A);
    while (iter != max_iter && residual > tol && residual != 0)
    {
        // residual must be in RES precision
        r = b - (A * denseM<RES>(x));
        residual = norm<RES>(r);

        // Using LU again in FACT precision to get correction
        c = LU_solver(Af, denseM<FACT>(r), P);
        Af = denseM<FACT>(A);
        x = x + c;
        iter++;
    }
    chrono::time_point end_time = chrono::steady_clock::now();
    chrono::duration<double> elapsed_time_seconds = end_time - start_time;
    cout << "Elapsed time: " << elapsed_time_seconds.count() << " seconds\n";
    cout << "The total iteration is " << iter << "\n";
    cout << "The error in the last iteration is " << residual << "\n";
    cout << "Iterative refinement succeeded!"
         << "\n"
         << "x = ";
    return x;
}
