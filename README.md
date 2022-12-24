# Iterative-refinement-with-three-precisions
## Main Goal
The C++ program *main.cpp* are focusing on demonstrating the class **denseM** which is storing dense matrix in a 1d floating-point number vector. Using the class and its member functions to implement LU-factorization and Cholesky-factorization including the linear system solver with these factorization. Eventually, using these factorizations to implement iterative refinement in three precisions, for getting a more accurate result than regular linear system solver and less time cost than fixed-precision iterative refinement. All declarations and implementations are in *denseM.hpp* header file. 
<br/>

The main implementations includes:
1. [Constructor by two any precision integers](#constructor-for-preallocation)
2. [Constructor by two any precision integers, and one any floating-point number type vector](#constructor-for-specific-matrix)
3. [Constructor by one any precision integer, one any floating-point number type vector, and one boolean variable (0 or 1)](#constructor-for-positive-definite-matrix)
4. [Constructor by a existed denseM object, change its integer/floating-point number precision](#constructor-for-change-precision)
5. [Member function "at"](#member-function-at)
6. [Member function overloading operator "[]"](#member-function-overloading-operator)
7. [Matrix addition](#matrix-addition)
8. [Matrix subtraction](#matrix-subtraction)
9. [Matrix multiplication](#matrix-multiplication)
10. [Infinity norm of matrix](#infinity-norm-of-matrix)
11. [LU-decomposition and solver](#lu-decomposition-and-solver)
12. [Cholesky-decomposition and solver](#cholesky-decomposition-and-solver)
13. [Iterative refinement]

<br/>

## Environment
This program contains several C++ header files: 

`<iostream>`, `<chrono>`, `<stdexcept>`, `<vector>`, `<string>`, `<algorithm>`, `<cmath>`

The integer number can use any signed or unsigned integer type, the width of the integer is depending on the matrix size the user inputs. <br/>
The floating-point number can use `float` and `double` by default. 
<br/>

## class denseM
Class *denseM* has private variables:
1. `INT rows_ = 0`. Recording the number of rows, default number is 0.
2. `INT cols_ = 0`. Recording the number of columns, default number is 0.
3. `vector<FLOAT> matrix_`. Recording the elements of the matrix in a floating-point number type vector.
4. `bool is_pos_def = 0`. Status of the matrix if or not it's positive definite.
<br/>
Every initialized denseM object will have the specific number of rows and columns in a integer type `INT`, and a 1d vector of the matrix in a floating-point number type `FLOAT`.
### Constructors:
1. ### Constructor for preallocation
Construct a new denseM object with zeros in the matrix, given the sizes. It is for preallocating the memory.
```cpp
// Preallocate 3x3 matrix
denseM<uint32_t, double> m0(3, 3);
cout << "The undefined matrix is: "
     << "\n"
     << m0 << "\n";
```
```
The undefined matrix is: 
0 0 0
0 0 0
0 0 0
```
2. ### Constructor for specific matrix
Construct a new denseM object with specific elements, given the sizes and a floating-point number type vector contains all the elements of the matrix.
```cpp
// Create 2x3 matrix with given vector
vector<double> v1 = {1.1234, 2.5, 3.3333, 4.19, 5, 6.2};
denseM<uint32_t, double> m1(2, 3, v1);
cout << "2x3 matrix m1 is: "
     << "\n"
     << m1 << "\n";
// Create 3x2 matrix with same vector
denseM<uint32_t, double> m2(3, 2, v1);
cout << "3x2 matrix m2 is: "
     << "\n"
     << m2 << "\n";
```
```
2x3 matrix m1 is:
1.1234 2.5 3.3333
4.19 5 6.2

3x2 matrix m2 is:
1.1234 2.5
3.3333 4.19
5 6.2
```
3. ### Constructor for positive definite matrix
Construct a new denseM object with the size, the matrix vector, and a boolean is-positive-definite checker. The user should know if the matrix they create is or not positive definite. Input 1 in the third argument if the user know the matrix is
positive definite, input 0 if it's not and it will create a regular square matrix.
```cpp
// Create a 3x3 matrix which is a positive definite matrix
vector<double> v2 = {4, 12, -16, 12, 37, -43, -16, -43, 98};
denseM<uint32_t, double> m3(3, v2, 1);
cout << "3x3 matrix m3 is: "
     << "\n"
     << m3 << "\n";
```
```
3x3 matrix m3 is:
4 12 -16
12 37 -43
-16 -43 98
```
4. ### Constructor for change precision
Construct a new denseM object based on a existed denseM object with another precision.
```cpp
// Change precision of m1 from double to float
vector<double> v1 = {1.1234, 2.5, 3.3333, 4.19, 5, 6.2};
denseM<uint32_t, double> m1(2, 3, v1);
denseM<uint32_t, float> m1_F(m1);
cout << "m1 is double precision is:"
     << "\n"
     << m1 << "\n";
cout << "m1 in float precision is: "
     << "\n"
     << m1_F << "\n";
```
```
m1 in double precision is:
1.1234 2.5 3.3333
4.19 5 6.2

m1 in float precision is:
1.123399972915649 2.5 3.333300113677979
4.190000057220459 5 6.199999809265137
```
5. ### Member function "at"
Show the value of matrix by given row number and column number. It will use the mathematical index of matrix. <br/>
For example, M.at(3,2) will show the value at 3rd row and 2nd column in matrix M.
```cpp
// Check the value in 1st row, 2nd column in m1 by "at"
vector<double> v1 = {1.1234, 2.5, 3.3333, 4.19, 5, 6.2};
denseM<uint32_t, double> m1(2, 3, v1);
cout << "2x3 matrix m1 is: "
     << "\n"
     << m1 << "\n";
cout << "m1's value in 1st row, 2nd column is: " << m1.at(1, 2) << "\n";
```
```
2x3 matrix m1 is:
1.1234 2.5 3.3333
4.19 5 6.2

m1's value in 1st row, 2nd column is: 2.5
```
6. ### Member function overloading operator "[]"
Overload [] operator for finding specific elements given the index. It has two versions, one is the value at the index is modifiable, another one is non-modifiable.
```cpp
vector<double> v1 = {1.1234, 2.5, 3.3333, 4.19, 5, 6.2};
denseM<uint32_t, double> m1(2, 3, v1);
// Check the value on index 4 in m1 by "[]"
cout << "m1's value on index 4 is: " << m1[4] << "\n";
// Change the value on index 4 in m1 by "[]"
m1[4] = 5.1234;
cout << "m1's changed value on index 4 is: " << m1[4] << "\n\n";
```
```
m1's value at index 4 is: 5
m1's changed value on index 4 is: 5.1234
```

7. ### Matrix addition
Overloaded binary operator + to calculate the summation of the first denseM object and the second denseM object. It uses the overloaded operator [] for direct-accessing denseM object and change the value for saving the memories.
```cpp
// Addition
          vector<double> a1 = {1, 2, 3, 4.5};
          vector<double> a2 = {-1, 2.5, -6, -0.4};
          denseM<uint32_t, double> A1(2, 2, a1);
          denseM<uint32_t, double> A2(2, 2, a2);
          cout << "matrix A1: "
               << "\n"
               << A1 << "matrix A2: "
               << "\n"
               << A2 << "\n";
          cout << "A1 + A2 = "
               << "\n"
               << A1 + A2 << "\n";
```
```
matrix A1:
1 2
3 4.5
matrix A2:
-1 2.5
-6 -0.4

A1 + A2 =
0 4.5
-3 4.1
```

8. ### Matrix subtraction
Overloaded binary operator - to calculate the subtraction of the first denseM object and the second denseM object. It uses the overloaded operator [] for direct-accessing denseM object and change the value for saving the memories.
```cpp
// Subtraction
          vector<double> a1 = {1, 2, 3, 4.5};
          vector<double> a2 = {-1, 2.5, -6, -0.4};
          denseM<uint32_t, double> A1(2, 2, a1);
          denseM<uint32_t, double> A2(2, 2, a2);
          cout << "matrix A1: "
               << "\n"
               << A1 << "matrix A2: "
               << "\n"
               << A2 << "\n";
          cout << "A1 - A2 = "
               << "\n"
               << A1 - A2 << "\n";
```
```
matrix A1:
1 2
3 4.5
matrix A2:
-1 2.5
-6 -0.4

A1 - A2 =
2 -0.5
9 4.9
```

9. ### Matrix multiplication
Overloaded binary operator * to calculate the multiplication of the first denseM object and the second denseM object. It uses the overloaded operator [] for direct-accessing denseM object and change the value for saving the memories.
```cpp
// Multiplication
          vector<double> a1 = {1, 2, 3, 4.5};
          vector<double> a2 = {-1, 2.5, -6, -0.4};
          denseM<uint32_t, double> A1(2, 2, a1);
          denseM<uint32_t, double> A2(2, 2, a2);
          cout << "matrix A1: "
               << "\n"
               << A1 << "matrix A2: "
               << "\n"
               << A2 << "\n";
          cout << "A1 * A2 = "
               << "\n"
               << A1 * A2 << "\n";
```
```
matrix A1:
1 2
3 4.5
matrix A2:
-1 2.5
-6 -0.4

A1 * A2 =
-13 1.7
-30 5.7
```
10. ### Infinity norm of matrix
Finding the largest row sum of the matrix.
```cpp
vector<double> a1 = {1, 2, 3, 4.5};
denseM<uint32_t, double> A1(2, 2, a1);
cout << "matrix A1: " << "\n" << A1 << "\n";
cout << "Infinity norm of A1: " << norm(A1) << "\n";
```
```
matrix A1:
1 2
3 4.5

Infinity norm of A1: 7.5
```
11. ### LU-decomposition and solver
Using LU-decomposition with partial pivoting. By the idea of Gaussian elimination, it will modify the denseM object A to the combination of L and U for saving memory (U will be the top-half of the modified A including the diagonal, L will be the bottom-half of the modified A excluding the diagonal, since the diagonal of L is always 1). Modify the integer vector P to record row swapping as in permutation matrix. It will return exit flag as integer to tell user the status of the result. Return 0 means success, return >0 means decomposition is completed but U is singular.
```cpp
// LU-decomposition
          vector<double> d1 = {5.23, 2.11, 3.15, 0, 1.67, 4.57, 10.111, 6.223, 0};
          denseM<uint32_t, double> D1(3, 3, d1);
          cout << "The original matrix before LU-decomposition is : "
               << "\n"
               << D1 << "\n";
          vector<uint32_t> P1 = {0, 1, 2};
          uint32_t exit_flag = LU_withPivot(D1, P1);
          cout << "After LU-decomposition, combined both L and U: "
               << "\n"
               << D1 << "\n"
               << "The exit_flag is: " << exit_flag
               << "\n";
```
```
The original matrix before LU-decomposition is :
5.23 2.11 3.15
0 1.67 4.57
10.111 6.223 0

After LU-decomposition, combined both L and U:
10.111 6.223 0
0 1.67 4.57
0.5172584314113342 -0.6640115081872652 6.184532592415803

The exit_flag is: 0
```
LU-decomposition solver uses the modified A and P from LU_withPivot(A,P), solve linear system $LUx = Pb$ by forward and backward substitution.
```cpp
// LU-solver
          // Need to reset D1
          denseM<uint32_t, double> D1_origin(3, 3, d1);
          vector<double> b1 = {1, 2, 3};
          denseM<uint32_t, double> B1(3, 1, b1);
          denseM<uint32_t, double> x1 = LU_solver(D1_origin, B1, P1);
          cout << "The result from LU solver is: "
               << "\n"
               << x1;
          denseM<uint32_t, double> D1_origin2(3, 3, d1);
          cout << "Verify the result by norm(b - Ax): "
               << "\n"
               << norm(B1 - (D1_origin2 * x1)) << "\n\n";
```
```
The result from LU solver is:
-0.2289842022255823
0.8541313303395248
0.1255143716264756
Verify the result by norm(b - Ax):
1.110223024625157e-16
```
12. ### Cholesky-decomposition and solver
Using Cholesky decomposition will decompose a positive definite matrix A into a lower-triangular matrix L and its conjugate transpose $L^T$, such that $L*L^T = A$.
$L$ can be found by formula: <br/>
$L_{i,j} = \sqrt{A_{i,j}-\sum_{k=1}^{j-1}L_{j,k}L^{*}_{j,k}}$ (for i = j) <br/>
$L_{i,j} = \frac{1}{L_{j,j}}(A_{i,j}-\sum_{k=1}^{j-1}L_{i,k}L^{*}_{j,k})$ (for i > j) <br/>
$L^T$ is the conjugate transpose of $L$, so they are stored as one single symmetric matrix as the return value of this function. It is generally faster when decompose a positive definite matrix than using LU-decomposition.
```cpp
// Cholesky decomposition
          vector<double> d2 = {4, 12, -16, 12, 37, -43, -16, -43, 98};
          denseM<uint32_t, double> D2(3, d2, 1);
          cout << "The original matrix before Cholesky-decomposition is : "
               << "\n"
               << D2 << "\n";
          denseM<uint32_t, double> D2_chole = cholesky(D2);
          cout << "After Cholesky-decomposition, combined both L and L_T: "
               << "\n"
               << D2_chole << "\n";
```
```
The original matrix before Cholesky-decomposition is :
4 12 -16
12 37 -43
-16 -43 98

After Cholesky-decomposition, combined both L and L_T:
2 6 -8
6 1 5
-8 5 3
```
Get the combination matrix A_chole of $L$ and $L^T$ from `cholesky(A,P)`, solve linear system $LL^Tx=b$ by forward and backward substitution.
```cpp
// Cholesky solver
          vector<double> b2 = {-22.2, -47.74, 213.12};
          denseM<uint32_t, double> B2(3, 1, b2);
          denseM<uint32_t, double> x2 = cholesky_solver(D2, B2);
          cout << "The result from Cholesky solver is: "
               << "\n"
               << x2;
          cout << "Verify the result by norm(b - Ax): "
               << "\n"
               << norm(B2 - (D2 * x2)) << "\n\n";
```
```
The result from Cholesky solver is:
1.245555555555659
2.182222222222194
3.33555555555556
Verify the result by norm(b - Ax):
2.842170943040401e-14
```
13. ### Iterative refinement
Using iterative refinement to solve linear system for less round off error. Including LU-decomposition and Cholesky decomposition to solve linear system inside iterations. Using three floating-point number precisions to accelerate decomposition by low-accuracy precision, and getting a more precise residual for updating the solution in each iteration by using a high-accuracy precision. Most variable are store in the middle precision. <br/>
To use this function, user needs to put one integer type and all three floating-point number types in the bracket with order from low accuracy to high accuracy. For example: `IR<uint32_t, float, double, double>(IR_A1, IR_B1)` <br/><br/>
If the matrix is regular dense matrix, it will use LU-decomposition. The tolerance is 1e-16.
```cpp
// IR for regular dense matrix
          vector<double> ir_a1 = {5.23, 2.11, 3.15, 0, 1.67, 4.57, 10.111, 6.223, 0};
          denseM<uint32_t, double> IR_A1(3, 3, ir_a1);
          vector<double> ir_b1 = {1, 2, 3};
          denseM<uint32_t, double> IR_B1(3, 1, ir_b1);
          // Three precisions are float, double and double
          denseM<uint32_t, double> result1 = IR<uint32_t, float, double, double>(IR_A1, IR_B1);
          cout << result1 << "\n";
```
```
Starting iterative refinement:
Elapsed time: 0.0002631 seconds
The total iteration is 9
The error in the last iteration is 0
Iterative refinement succeeded!
x = -0.2289842022255823
0.8541313303395248
0.1255143716264757
```
If the matrix is positive definite matrix, it will use Cholesky-decomposition. The tolerance is 1e-14.
```cpp
// IR for positive definite matrix
          vector<double> ir_a2 = {4, 12, -16, 12, 37, -43, -16, -43, 98};
          denseM<uint32_t, double> IR_A2(3, ir_a2, 1);
          vector<double> ir_b2 = {-22.2, -47.74, 213.12};
          denseM<uint32_t, double> IR_B2(3, 1, ir_b2);
          // Three precisions are float, double and double
          denseM<uint32_t, double> result2 = IR<uint32_t, float, double, double>(IR_A2, IR_B2);
          cout << result2 << "\n";
```
```
Starting iterative refinement:
Elapsed time: 0.000208 seconds
The total iteration is 9
The error in the last iteration is 7.105427357601002e-15
Iterative refinement succeeded!
x = 1.245555555555548
2.182222222222224
3.335555555555555
```
