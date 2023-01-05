## Main Goal
The C++ program *main.cpp* are focusing on demonstrating the class **denseM** which is storing dense matrix in a 1d floating-point number vector. Using the class and its member functions to implement LU-factorization and Cholesky-factorization including the linear system solver with these factorization. Eventually, using these factorizations to implement iterative refinement in three precisions, for getting a more accurate result than regular linear system solver and less time cost than fixed-precision iterative refinement. In addition, it also includes Generalized minimal residual method(GMRES) for GMRES-iterative refinement in mixed precisions. GMRES-IR can handle situations where the matrix has high condition numbers and cannot get a good approximation in the regular iterative refinement. All declarations and implementations are in *denseM.hpp* header file. 
<br/>

The main implementations includes:
1. [Constructor by two unsigned 64bits integers](#constructor-for-preallocation)
2. [Constructor by two unsigned 64bits integers, and one any floating-point number type vector](#constructor-for-specific-matrix)
3. [Constructor by one any precision integer, one any floating-point number type vector, and one boolean variable (0 or 1)](#constructor-for-positive-definite-matrix)
4. [Constructor by one unsigned 64bits integers](#constructor-for-identity-matrix)
5. [Constructor by two unsigned 64bits integers, and an existed file contains matrix](#constructor-for-building-matrix-by-reading-file)
6. [Constructor by a existed denseM object, change its integer/floating-point number precision](#constructor-for-change-precision)
7. [LU-decomposition and solver](#lu-decomposition-and-solver)
8. [Cholesky-decomposition and solver](#cholesky-decomposition-and-solver)
9. [Generalized minimal residual method(GMRES)](#gmres)
10. [Iterative refinement](#mixed-precision-iterative-refinement)
11. [GMRES-iterative refinement](#mixed-precision-gmres-ir)

<br/>

## Environment
This program contains several C++ header files: 

`<iostream>`, `<chrono>`, `<stdexcept>`, `<vector>`, `<string>`, `<algorithm>`, `<cmath>`, `<fstream>`

The floating-point number type can use `float` and `double` by default. 
<br/>

## class `denseM`
Templated Class `denseM<FLOAT>` has private variables:
1. `uint64_t rows_ = 0`. Recording the number of rows, default number is 0.
2. `uint64_t cols_ = 0`. Recording the number of columns, default number is 0.
3. `vector<FLOAT> matrix_`. Recording the elements of the matrix in a floating-point number type vector.
4. `bool is_pos_def = 0`. Status of the matrix if or not it's positive definite.
<br/>
Every initialized denseM object will have the specific number of rows and columns in the integer type `uint64_t`, and a templated 1d vector of the matrix, with the floating-point number type `FLOAT`.

## Constructors:
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

4. ### Constructor for identity matrix
Construct a new denseM with an integer, create a nxn identity matrix.
```cpp
// Create a 3x3 matrix which is a positive definite matrix
          denseM<double> m4(3);
          cout << "3x3 identity matrix m4 is: "
               << "\n"
               << m4 << "\n";
```
```
3x3 identity matrix m4 is:
1 0 0
0 1 0 
0 0 1
```
5. ### Constructor for building matrix by reading file
Construct a new denseM object given size and a file name which contains matrix. The matrix number in the file can be separated by comma or white space. (If the file has extension "csv", using comma to separate numbers can have a better format inside the file) </br>
User can also use function `output` to write a file.
```cpp
// "matrix.csv"
// 1,2,3
// 4,5,6
// 7,8,9
// Create a 3x3 matrix given a file
          denseM<double> f1(3, 3, "matrix.csv");
          cout << "3x3 matrix f1 is: "
               << "\n"
               << f1 << "\n";
          // Output and write the matrix into a file
          f1.output("m2.csv");
```
```
3x3 matrix f1 is: 
1 2 3
4 5 6
7 8 9

m2.csv generated successfully.
```

6. ### Constructor for change precision
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

## Decompositions and solvers
1. ### LU-decomposition and solver
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
2. ### Cholesky-decomposition and solver
Using Cholesky decomposition will decompose a positive definite matrix A into a lower-triangular matrix L and its conjugate transpose $L^T$, such that $L*L^T = A$. 
$L$ can be found by formula: 
1. $L_{i,j} = \sqrt{A_{i,j}-\sum\nolimits_{k=1}^{j-1} L_{j,k}L^{*}_{j,k}}$ (for i = j) 
2. $L_{i,j} = \frac{1}{L_{j,j}}(A_{i,j}-\sum\nolimits_{k=1}^{j-1}L_{i,k}L^{*}_{j,k})$ (for i > j) <br>
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
3. ### GMRES
Generalized minimal residual method (GMRES) is an iterative method for solving linear system. There are three arguments, the linear system A and b, and an initial guess $x_{0}$. The method approximates solution of the linear system by using Krylov subspace vectors $K_{n}(A,r_{0})= span \,\{r_{0},Ar_{0},A^{2}r_{0},\ldots ,A^{n-1}r_{0}\}$, where $r_{0} = b-Ax_{0}$. </br>
Starting iteration on the first vector. We need to use the Arnoldi iteration is used to find each vector in its orthonormal form in $q_1,\ldots,q_n$ stores in Q and produces a Hessenberg matrix H. Then we need to solve the least square problem $||Hy-||r||e_{1}||$ , by finding the y that can minimize it. Using Given rotation matrix to help converting H into a upper triangular matrix R, and updating $||r_{0}||e_{1}$ by given rotation matrix's elements to get a vector $g_{n}$. Solve $Ry=g_{n}$ to get y. </br>
Finally, update the initial guess $x_0$ by $x = x_{0} + Qy$. If the new x doesn't meet the criteria, it will recall the function and using x as the new initial guess.
```cpp
          vector<float> gmres_a1 = {5.23, 2.11, 3.15, 0, 1.67, 4.57, 10.111, 6.223, 0};
          denseM<float> GMRES_A1(3, 3, gmres_a1);
          vector<float> gmres_b1 = {1, 2, 3};
          denseM<float> GMRES_B1(3, 1, gmres_b1);
          vector<float> guess = {1.24, 2.22, 3};
          denseM<float> Guess(3, 1, guess);
          denseM<float> G_result1 = GMRES<float>(GMRES_A1, GMRES_B1, Guess);
          cout << "GMRES result is: "
               << "\n"
               << G_result1 << "\n";
```
```
GMRES result is:
-0.2289841175079346
0.8541311025619507
0.1255145072937012
```

## Different Iterative refinements:
1. ### Mixed-precision iterative refinement
Using iterative refinement to solve linear system for less round off error. In the beginning of the iteration, using LU-decomposition or Cholesky decomposition to solve linear system. Using three floating-point number precisions to accelerate the decomposition by low-accuracy precision, and getting a more precise residual for updating the solution in each iteration by using a high-accuracy precision. Most variable are store in the middle precision. <br/>
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

2. ### Mixed-precision GMRES-IR
Preconditioned GMRES-based Iterative Refinement is for solving linear system even with some ill-conditioned matrices. In the beginning of the iteration, using LU-decomposition or Cholesky decomposition to solve linear system. Using three floating-point number precisions to accelerate the decomposition by low-accuracy precision, and getting a more precise residual for updating the solution in each iteration by using a high-accuracy precision. Inside the iterative refinement, the correction c will be updated in medium-accuracy precision by GMRES after preconditioned by $\hat{U}^{-1}\hat{L}^{-1}$, where $\hat{U}$ and $\hat{L}$ are got from LU-decomposition in low-accuracy precision.
```cpp
          vector<float> gmres_a1 = {5.23, 2.11, 3.15, 0, 1.67, 4.57, 10.111, 6.223, 0};
          denseM<float> GMRES_A1(3, 3, gmres_a1);
          vector<float> gmres_b1 = {1, 2, 3};
          denseM<float> GMRES_B1(3, 1, gmres_b1);
          // GMRES-IR
          denseM<double> GMRES_IR_result = GMRES_IR<float, double, double>(GMRES_A1, GMRES_B1);
          cout << GMRES_IR_result << "\n";
```
```
Starting iterative refinement:
Elapsed time: 0.0005964 seconds
The total iteration is 2
The error in the last iteration is 0
Iterative refinement succeeded!
x = -0.2289841815397392
0.854131292168907
0.1255143888812452
```
