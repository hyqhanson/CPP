## Main Goal
The C++ program `main.cpp` are focusing on demonstrating the class `denseM` which is storing dense matrix in a 1d floating-point number vector. Using the class and its member functions to implement LU-factorization and Cholesky-factorization including the linear system solver with these factorization. Eventually, using these factorizations to implement iterative refinement in three precisions, for getting a generally more accurate result than regular linear system solver and less time-cost than fixed-precision iterative refinement. In addition, it also includes Generalized minimal residual method(GMRES) for preconditioned GMRES-iterative refinement in mixed precisions. GMRES-IR can generally handle situations where the matrix has high condition numbers with a better solution or less time-cost. All declarations and implementations are in *denseM.hpp* header file. 
<br/>

### class `denseM`
Templated Class `denseM<FLOAT>` has private variables:
1. `uint64_t rows_ = 0`. Recording the number of rows, default number is 0.
2. `uint64_t cols_ = 0`. Recording the number of columns, default number is 0.
3. `vector<FLOAT> matrix_`. Recording the elements of the matrix in a floating-point number type vector.
4. `bool is_pos_def = 0`. Status of the matrix if or not it's positive definite.
<br/>
Every initialized denseM object will have the specific number of rows and columns in the integer type `uint64_t`, and a templated 1d vector of the matrix, with the floating-point number type `FLOAT`.
<br>

### Environment
This program contains several C++ header files: 

`<iostream>`, `<chrono>`, `<stdexcept>`, `<vector>`, `<string>`, `<algorithm>`, `<cmath>`, `<cstring>`, `<fstream>`

The floating-point number type can use `float` and `double` by default. 
<br/>
<br>

## Menu
1. [Matrix file format](#matrix-file-format)
2. [How to use this program](#how-to-use-this-program)

The main functions inside this program includes:
1. [LU-decomposition and solver](#lu-decomposition-and-solver)
2. [Cholesky-decomposition and solver](#cholesky-decomposition-and-solver)
3. [Generalized minimal residual method(GMRES)](#gmres)
4. [Iterative refinement](#mixed-precision-iterative-refinement)
5. [GMRES-iterative refinement](#mixed-precision-gmres-ir)

<br/> 



## Matrix file format:
This program requires the user to input matrices by files (.txt or .csv). Each element in the matrix should be a floating-point number. If there is non-digit elements inside, it will be ignored. The user also needs to know the size of the matrix as one of the input argument later around. <br>
All the input and output matrices can be in .txt or .csv file format. Whitespace or comma is required for separating each elements. <br><br>
Input matrices in files can be formatted in the regular way (have a new line for each row):
```
Ex: A 4x4 matrix in Ex1.txt file:
24.2 2 5.3 3
1.56 6 6 7.8
1 23 45.3 3
3 3.4 7 2.222
```
Or just write them as a long vector in one row:
```
Ex: A 2x2 matrix in Ex2.txt file:
2 4.2 7 1.2
```
<br>

Output matrices in files will be formatted in the regular way. 
```
Ex: A 3x1 matrix in LU_result.csv file. Generated by `main LU_solver double 3 LU1.csv LU2.csv LU_result.csv`:
-0.228984
0.854131
0.125514
```
 
<br>

## How to use this program:
The user needs to run `main.cpp` at least once to compile the program. It will show a message:
```
You need to put some arguments in command line.
``` 
Then the user needs to use commend line in terminal for running this program.  Make sure the current path is in the directory where the source files and testing files are. <br>
For **powershell:** type “.\” before the program name(it is `main` in here), followed by other commend line arguments separated by whitespace. The program name doesn’t need its extension, all the other files' name needs its extension (.txt or .csv). <br>
```
Ex:  .\main LU_solver double 3 LU1.csv LU2.csv LU_result.csv
```
For **cwd:** starting with the program name directly(it is `main` in here), followed by other commend line arguments separated by whitespace. The program name doesn’t need its extension, all the other files' name needs its extension (.txt or .csv).
```
Ex:  main LU_solver double 3 LU1.csv LU2.csv LU_result.csv
```
The details of arguments will be introduced in each function section.
<br>
<br>

## Main functions:
### Decompositions and solvers
1. ### LU-decomposition and solver
The LU-decomposition with partial pivoting is implemented in `LU_withPivot`. By the idea of Gaussian elimination, it will modify the denseM object A to the combination of L and U for saving memory (U will be the top-half of the modified A including the diagonal, L will be the bottom-half of the modified A excluding the diagonal, since the diagonal of L is always 1). Modify a integer vector P such that P = $\{0,1,2,...,n-1\}$ for nxn matrix A, to record the row-swapping in the permutation matrix. `LU_withPivot` will return exit flag in integer to tell the status of the result. Return 0 means LU-decomposition succeed, return >0 means decomposition is completed but U is singular and cannot be used to solve a linear system. <br> 
If exit flag is 0, LU-decomposition solver `LU_solver` can use the modified A and P from `LU_withPivot(A,P)`, to solve linear system $LUx = Pb$ by forward and backward substitution. <br>
For solving a linear system by LU factorization using `LU_solver`, user needs to input 7 arguments in command line. For example:
```
main LU_solver double 3 LU1.csv LU2.csv LU_result.csv
```
<br>

### Arguments explanation:
**1st argument `main`:** the program's name <br>
**2nd argument `LU_solver`:** the function's name, using LU factorization to solve linear system <br>
**3rd argument:** the floating-point number's precision, it can be `single` or `double`<br>
**4th argument:** the size n of the nxn matrix A and the nx1 matrix b, it needs to be an integer <br>
**5th argument:** the filename of input nxn matrix A<br>
**6th argument:** the filename of input nx1 matrix b<br>
**7th argument:** the filename of output nx1 matrix x. It can be replaced by `N` if the user does not want an output file, it will print out the result directly. <br>

An example when the 7th argument is `N`, with the corresponding result
```
main LU_solver double 3 LU1.csv LU2.csv N

The solution of the linear system by LU decomposition in double precision is: 
-0.228984 
0.854131 
0.125514 

Verify the result by norm(b - Ax): 
1.11022e-16
```
<br>

### Handle Exceptions:
- If the number of total arguments is not 7, it will print the message
```
The number of arguments is incorrect. You need 7 arguments in total. 
Example: main LU_solver double 3 LU1.csv LU2.csv LU_result.csv
```
- If a wrong function name is entered in the 2nd argument, it will print the message
```
The second argument should be one of the function below:
LU_solver
cholesky_solver
GMRES
IR
GMRES_IR
```
- If the precision in the 3rd argument is not `single` nor `double`, it will print the message
```
Please use either single or double precision in the 3rd argument.
```
- If the size in the 4th argument contains non-numerical characters, or it is a negative integer, it will print the message
```
The size of the matrix needs be an positive integer!
```
- If the size in the 4th argument is 0, it will print the message
```
The matrix size can not be zero!
```
- When reading the file in the 5th or 6th argument has problems, it will print the message
```
There is a problem when opening the file **filename**
```
- If any matrix in the file mismatches the size given in the 4th argument, it will print the message
```
File: **filename** has issues. 
The matrix vector does not have the same amount of elements as the given size.
```
- If the output file name in the 7th argument is invalid, such as a file name without an extension '.csv' or '.txt'(except `N`), or a file name has no character in front of the extension, it will print the message
```
The extension of the output file should be '.csv' or '.txt', and the name of the file must have some character in front of the extension.
```
- It will print the message if the matrix A from the 5th argument is singular, try with:
```
main LU_solver double 2 singularM.csv Err2.csv N
```
```
LU-decomposition is complete, the upper-triangular matrix's diagonal number in row2 is 0, it is singular.
Error: The LU-factorization exit-flag is not 0. This linear system cannot be solved.
```
<br>

2. ### Cholesky-decomposition and solver
Using Cholesky decomposition will decompose a positive definite matrix A into a lower-triangular matrix L and its conjugate transpose $L^T$, such that $L*L^T = A$.
$L$ can be found by formula: <br/>
$L_{i,j} = \sqrt{A_{i,j}-\sum\nolimits_{k=1}^{j-1}L_{j,k}L^{*}_{j,k}}$ (for i = j) <br/>
$L_{i,j} = \frac{1}{L_{j,j}}(A_{i,j}-\sum\nolimits_{k=1}^{j-1}L_{i,k}L^{*}_{j,k})$ (for i > j) <br/>
$L^T$ is the conjugate transpose of $L$, so they are stored as one single symmetric matrix as the return value of this function. It is generally faster when decompose a positive definite matrix than using LU-decomposition. <br>
Get the combination matrix A_chole of $L$ and $L^T$ from `cholesky(A,P)`, solve linear system $LL^Tx=b$ by forward and backward substitution in `cholesky_solver`. <br>
For solving a linear system by Cholesky factorization using `cholesky_solver`, user needs to input 7 arguments in command line. For example:
```
main cholesky_solver double 3 C1.csv C2.csv C_result.csv
```
<br>

### Arguments explanation:
**1st argument `main`:** the program's name <br>
**2nd argument `cholesky_solver`:** the function's name, using Cholesky factorization to solve linear system <br>
**3rd argument:** the floating-point number's precision, it can be `single` or `double`<br>
**4th argument:** the size n of the nxn matrix A and the nx1 matrix b, it needs to be an integer <br>
**5th argument:** the filename of input nxn matrix A<br>
**6th argument:** the filename of input nx1 matrix b<br>
**7th argument:** the filename of output nx1 matrix x. It can be replaced by `N` if the user does not want an output file, it will print out the result directly. <br>

An example when the 7th argument is `N`, with the corresponding result
```
main cholesky_solver double 3 C1.csv C2.csv N

The solution of the linear system by Cholesky decomposition in double precision is: 
1.24556 
2.18222 
3.33556 

Verify the result by norm(b - Ax): 
2.84217e-14
```
<br>

### Handle Exceptions:
- If the number of total arguments is not 7, it will print the message
```
The number of arguments is incorrect. You need 7 arguments in total. 
Example: main cholesky_solver double 3 C1.csv C2.csv C_result.csv
```
- If a wrong function name is entered in the 2nd argument, it will print the message
```
The second argument should be one of the function below:
LU_solver
cholesky_solver
GMRES
IR
GMRES_IR
```
- If the precision in the 3rd argument is not `single` nor `double`, it will print the message
```
Please use either single or double precision in the 3rd argument.
```
- If the size in the 4th argument contains non-numerical characters, or it is a negative integer, it will print the message
```
The size of the matrix needs be an positive integer!
```
- If the size in the 4th argument is 0, it will print the message
```
The matrix size can not be zero!
```
- When reading the file in the 5th or 6th argument has problems, it will print the message
```
There is a problem when opening the file **filename**
```
- If any matrix in the file mismatches the size given in the 4th argument, it will print the message
```
File: **filename** has issues. 
The matrix vector does not have the same amount of elements as the given size.
```
- If the output file name in the 7th argument is invalid, such as a file name without an extension '.csv' or '.txt'(except `N`), or a file name has no character in front of the extension, it will print the message
```
The extension of the output file should be '.csv' or '.txt', and the name of the file must have some character in front of the extension.
```
- It will print the message if the matrix A from 5th argument is non-symmetric:
```
Error: This matrix is not symmetric, cannot be decomposed by Cholesky factorization
```
If the matrix A is symmetric but not positive definite, try with: 
```
main cholesky_solver double 2 Err1.csv Err2.csv N
```
```
Error: This matrix is not positive definite, cannot be decomposed by Cholesky factorization
```
<br>

3. ### GMRES
Generalized minimal residual method (GMRES) is an iterative method for solving linear system. There are three arguments, the linear system A and b, and an initial guess $x_{0}$. The method approximates solution of the linear system by using Krylov subspace vectors $K_{n}(A,r_{0})= span \,\{r_{0},Ar_{0},A^{2}r_{0},\ldots ,A^{n-1}r_{0}\}$, where $r_{0} = b-Ax_{0}$. </br>
Starting iteration on the first vector. We need to use the Arnoldi iteration is used to find each vector in its orthonormal form in $q_{0},\ldots,q_{n-1}$ stores in Q and produces a Hessenberg matrix H. The detailed steps of Arnoldi iteration is:
> $r_{0} = b - Ax_{0}$
>
> $q_{0} = \frac{r_{0}}{||r_{0}||}$
>
> for n = (1,2,3,...,max_iter)
>
>> $v = Aq_{0}$
>>
>> for m = (0,1,2,...,n-1)
>>
>>> $h_{m+1,n} = (q_{m})'*v$
>>>
>>> $v = v - h_{m+1,n}*q_{m}$
>>>
>> End for
>>
>> $h_{n+1,n} = ||v||$
>> 
>> $q_{n} = \frac{v}{h_{n+1,n}}$
>>
> End for

Then we need to solve the least square problem $||Hy-||r||e_{1}||$ , by finding the y that can minimize it. Using Given rotation matrix to help converting H into a upper triangular matrix R, and updating $||r_{0}||e_{1}$ by given rotation matrix's elements to get a vector $g_{n}$. Solve $Ry=g_{n}$ to get y. </br>
Finally, update the initial guess $x_0$ by $x = x_{0} + Qy$. If the new x doesn't meet the criteria, it will recall the function and using x as the new initial guess. <br>
For solving a linear system by GMRES using `GMRES`, user needs to input 8 arguments in command line. For example:
```
main GMRES double 3 G1.csv G2.csv guess.csv G_result.csv
```
<br>

### Arguments explanation:
**1st argument `main`:** the program's name <br>
**2nd argument `GMRES`:** the function's name, using GMRES to solve linear system <br>
**3rd argument:** the floating-point number's precision, it can be `single` or `double`<br>
**4th argument:** the size n of the nxn matrix A and the nx1 matrix b, it needs to be an integer <br>
**5th argument:** the filename of input nxn matrix A<br>
**6th argument:** the filename of input nx1 matrix b<br>
**7th argument:** the filename of input nx1 matrix of initial guess. It can be replaced by `N` if the user does not want to provide an initial guess, it will use a default nx1 zero matrix as the initial guess instead. <br>
**8th argument:** the filename of output nx1 matrix x. It can be replaced by `N` if the user does not want an output file, it will print out the result directly. <br>

An example when the 7th and 8th argument are both `N`, with the corresponding result
```
main GMRES double 3 G1.csv G2.csv N N

The solution of the linear system by GMRES in double precision is:
158.374 
-5.62334
3.90327

Verify the result by norm(b - Ax):
2.16005e-12
```
Another example when the 7th argument is a given initial guess
```
main GMRES double 3 G1.csv G2.csv guess.csv N

The solution of the linear system by GMRES in double precision is: 
158.374 
-5.62334 
3.90327 

Verify the result by norm(b - Ax): 
1.13687e-13
```
*A proper initial guess may result a more accurate solution.*

<br>

### Handle Exceptions:
- If the number of total arguments is not 8, it will print the message
```
The number of arguments is incorrect. You need 8 arguments in total.
Example: main GMRES double 3 G1.csv G2.csv guess.csv G_result.csv
```
- If a wrong function name is entered in the 2nd argument, it will print the message
```
The second argument should be one of the function below:
LU_solver
cholesky_solver
GMRES
IR
GMRES_IR
```
- If the precision in the 3rd argument is not `single` nor `double`, it will print the message
```
Please use either single or double precision in the 3rd argument.
```
- If the size in the 4th argument contains non-numerical characters, or it is a negative integer, it will print the message
```
The size of the matrix needs be an positive integer!
```
- If the size in the 4th argument is 0, it will print the message
```
The matrix size can not be zero!
```
- When reading the file in the 5th, 6th or 7th argument(if the 7th argument is not `N`) has problems, it will print the message
```
There is a problem when opening the file **filename**
```
- If any matrix in the file mismatches the size given in the 4th argument, it will print the message
```
File: **filename** has issues. 
The matrix vector does not have the same amount of elements as the given size.
```
- If the output file name in the 8 th argument is invalid, such as a file name without an extension '.csv' or '.txt'(except `N`), or a file name has no character in front of the extension, it will print the message
```
The extension of the output file should be '.csv' or '.txt', and the name of the file must have some character in front of the extension.
```
<br>

### Different Iterative refinements:
1. ### Mixed-precision iterative refinement
Using iterative refinement to solve linear system for less round off error. Before the iteration, using LU-decomposition or Cholesky decomposition to get the first solution $x_{0}$ of this linear system. Then starting the iteration by calculating the residual $r=||b-Ax_{0}||$. Find the correction $c$ by solving another linear system $Ac=r$, and update x by $x=x+c$. This loop stops either reaches the maximum iteration number(10000) or x reaches the desired tolerance(1e-14). There are three positions for putting precisions. Accelerating the decomposition by using the 1st precision (normally with a low-accuracy precision), and getting a more precise residual for updating the solution in each loop by using the 3rd precision (normally with a high-accuracy precision). Most variable are store in the 2nd precision(medium-accuracy precision). Users can choose precisions(single or double) in each position freely.  <br/>
To use the function `IR`, user needs to input 10 arguments in command line. For example:
```
main IR single double double 3 LU1.csv LU2.csv IR_result.csv LU
```

<br>

### Arguments explanation:
**1st argument `main`:** the program's name <br>
**2nd argument `IR`:** the function's name, using mixed precision IR to solve linear system <br>
**3rd argument:** the floating-point number's precision used for LU/Cholesky-factorization, it can be `single` or `double`<br>
**4th argument:** the floating-point number's precision used for the solution x and other trivial matrices or numbers, it can be `single` or `double`<br>
**5th argument:** the floating-point number's precision used for residual, it can be `single` or `double`<br>
**6th argument:** the size n of the nxn matrix A and the nx1 matrix b, it needs to be an integer <br>
**7th argument:** the filename of input nxn matrix A<br>
**8th argument:** the filename of input nx1 matrix b<br>
**9th argument:** the filename of output nx1 matrix x. It can be replaced by `N` if the user does not want an output file, it will print out the result directly. <br>
**10th argument:** the choice of using LU-factorization or Cholesky-factorization inside the iterative refinement, it can be `LU` or `C` <br>

If the matrix A is a regular dense matrix, it will use LU-decomposition. 
An example with 10th argument is `LU`:
```
main IR single double double 4 IR_LU1.txt IR_LU2.txt N LU 

Starting iterative refinement: 
Using LU decomposition
Elapsed time: 0.0766016 seconds
The total iteration is 10000
Iterative refinement succeeded!
The solution of the linear system by IR is:      
-0.697691
4.7591
8.28142
-3.74779

Verify the result by norm(b - Ax):
1.77636e-14
```
If the matrix A is positive definite matrix, it will use Cholesky-decomposition for a possibly less time-costing. 
```
main IR single double double 3 C1.csv C2.csv N C 

Starting iterative refinement: 
Using Cholesky decomposition
Elapsed time: 0.0007016 seconds
The total iteration is 9
Iterative refinement succeeded!
The solution of the linear system by IR is:      
1.24556
2.18222
3.33556

Verify the result by norm(b - Ax):
7.10543e-15
```
<br>

### Handle Exceptions:
- If the number of total arguments is not 10, it will print the message
```
The number of arguments is incorrect. You need 10 arguments in total.
Example: main IR single double double 4 IR_LU1.txt IR_LU2.txt N LU 
```
- If a wrong function name is entered in the 2nd argument, it will print the message
```
The second argument should be one of the function below:
LU_solver
cholesky_solver
GMRES
IR
GMRES_IR
```
- If the precision in the 3rd/4th/5th argument is not `single` nor `double`, it will print the message
```
Please use either single or double precision in the 3rd/4th/5th argument.
```
- If the size in the 6th argument contains non-numerical characters, or it is a negative integer, it will print the message
```
The size of the matrix needs be an positive integer!
```
- If the size in the 6th argument is 0, it will print the message
```
The matrix size can not be zero!
```
- When reading the file in the 7th/8th argument has problems, it will print the message
```
There is a problem when opening the file **filename**
```
- If any matrix in the file mismatches the size given in the 6th argument, it will print the message
```
File: **filename** has issues. 
The matrix vector does not have the same amount of elements as the given size.
```
- If the output file name in the 9th argument is invalid, such as a file name without an extension '.csv' or '.txt'(except `N`), or a file name has no character in front of the extension, it will print the message
```
The extension of the output file should be '.csv' or '.txt', and the name of the file must have some character in front of the extension.
```
- If the 10th argument is not `LU` nor `C`, it will print the message
```
You need to choose either 'LU' or 'C' in the last argument
```
- If the 10th argument is `LU`, but the matrix A is singular. It will print the message
```
Error: The LU-factorization exit-flag is not 0. This linear system cannot be solved.
```
- If the 10th argument is `C`, but the matrix A is not positive definite. It will print the message if the matrix is non-symmetric:
```
Starting iterative refinement: 
Using Cholesky decomposition
Error: This matrix is not symmetric, cannot be decomposed by Cholesky factorization
```
If the matrix is symmetric but not positive definite, try with: 
```
main IR single double double 2 Err1.csv Err2.csv N C
```
```
Starting iterative refinement: 
Using Cholesky decomposition
Error: This matrix is not positive definite, cannot be decomposed by Cholesky factorization
```
<br>

2. ### Mixed-precision GMRES-IR
Preconditioned GMRES-based Iterative Refinement is for solving linear system even with some ill-conditioned matrices. Similar stricture as the regular mixed-precision IR. using LU-decomposition to get the first solution $x_{0}$ of this linear system. Then starting the iteration by calculating the residual $r=||b-Ax_{0}||$. The linear system of getting correction $c$ will firstly modified by preconditioned system $\hat{U}^{-1}\hat{L}^{-1}P$, where $\hat{U}$, $\hat{L}$ and $P$ is calculated by partial pivoting LU-decomposition on A in a lower-accuracy precision. Solving the new linear system $\hat{U}^{-1}\hat{L}^{-1}PAc=\hat{U}^{-1}\hat{L}^{-1}Pr$, and update x by $x=x+c$. The purpose of preconditioned system is for reducing the condition number. This loop stops either reaches the maximum iteration number(10000) or x reaches the desired tolerance(1e-14). There are three positions for putting precisions. Accelerating the decomposition by using the 1st precision (normally with a low-accuracy precision), and getting a more precise residual for updating the solution in each loop and calculating $\hat{U}^{-1}\hat{L}^{-1}P$ by using the 3rd precision (normally with a high-accuracy precision). Most variable are store in the 2nd precision(medium-accuracy precision). Users can choose precisions(single or double) in each position freely.  <br/>
To use the function `GMRES_IR`, user needs to input 9 arguments in command line. For example:
```
main GMRES_IR single double double 3 G_IR1.csv G_IR2.csv G_IR_result.csv
```

<br>

### Arguments explanation:
**1st argument `main`:** the program's name <br>
**2nd argument `GMRES_IR`:** the function's name, using mixed precision IR to solve linear system <br>
**3rd argument:** the floating-point number's precision used for LU-factorization, it can be `single` or `double`<br>
**4th argument:** the floating-point number's precision used for the solution x and other trivial matrices or numbers, it can be `single` or `double`<br>
**5th argument:** the floating-point number's precision used for residual, it can be `single` or `double`<br>
**6th argument:** the size n of the nxn matrix A and the nx1 matrix b, it needs to be an integer <br>
**7th argument:** the filename of input nxn matrix A<br>
**8th argument:** the filename of input nx1 matrix b<br>
**9th argument:** the filename of output nx1 matrix x. It can be replaced by `N` if the user does not want an output file, it will print out the result directly. <br>

If the matrix A is an ill-conditioned dense matrix, try `GMRES_IR` without output a file
```
main GMRES_IR single double double 3 G_IR1.csv G_IR2.csv N 

Starting GMRES-Iterative refinement: 
Elapsed time: 0.0038251 seconds
The total iteration is 14
GMRES-Iterative refinement succeeded!
The solution of the linear system by GMRES_IR is:
-5.9988e+24
-6e+24
3e+24

Verify the result by norm(b - Ax):
0.01
```

<br>

### Handle Exceptions:
- If the number of total arguments is not 9, it will print the message
```
The number of arguments is incorrect. You need 9 arguments in total.
Example: main GMRES_IR single double double 3 G_IR1.csv G_IR2.csv G_IR_Result.csv
```
- If a wrong function name is entered in the 2nd argument, it will print the message
```
The second argument should be one of the function below:
LU_solver
cholesky_solver
GMRES
IR
GMRES_IR
```
- If the precision in the 3rd/4th/5th argument is not `single` nor `double`, it will print the message
```
Please use either single or double precision in the 3rd/4th/5th argument.
```
- If the size in the 6th argument contains non-numerical characters, or it is a negative integer, it will print the message
```
The size of the matrix needs be an positive integer!
```
- If the size in the 6th argument is 0, it will print the message
```
The matrix size can not be zero!
```
- When reading the file in the 7th/8th argument has problems, it will print the message
```
There is a problem when opening the file **filename**
```
- If any matrix in the file mismatches the size given in the 6th argument, it will print the message
```
File: **filename** has issues. 
The matrix vector does not have the same amount of elements as the given size.
```
- If the matrix A from file in the 7th argument is singular, try with:
```
main GMRES_IR single double double 2 singularM.csv Err2.csv N 
```
It will print the message
```
Starting GMRES-Iterative refinement: 
LU-decomposition is complete, the upper-triangular matrix's diagonal number in row2 is 0, it is singular.
Error: The LU-factorization exit-flag is not 0. This linear system cannot be solved.
```
- If the output file name in the 9th argument is invalid, such as a file name without an extension '.csv' or '.txt'(except `N`), or a file name has no character in front of the extension, it will print the message
```
The extension of the output file should be '.csv' or '.txt', and the name of the file must have some character in front of the extension.
```



