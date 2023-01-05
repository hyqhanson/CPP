#include <iostream>
#include <chrono>
#include <stdexcept>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include "denseM.hpp"
using namespace std;

int main()
{
     try
     {
          cout.precision(16);
          // 1. Constructors
          // Preallocate 3x3 matrix
          denseM<double> m0(3, 3);
          cout << "The undefined matrix is: "
               << "\n"
               << m0 << "\n";
          // Exception occurs when sizes are non-positive
          // Example:
          // denseM<double> err1(3, 0);

          // Create 2x3 matrix with given vector
          vector<double> v1 = {1.1234, 2.5, 3.3333, 4.19, 5, 6.2};
          denseM<double> m1(2, 3, v1);
          cout << "2x3 matrix m1 is: "
               << "\n"
               << m1 << "\n";
          // Create 3x2 matrix with same vector
          denseM<double> m2(3, 2, v1);
          cout << "3x2 matrix m2 is: "
               << "\n"
               << m2 << "\n";
          // Create a 3x3 matrix which is a positive definite matrix
          vector<double> v2 = {4, 12, -16, 12, 37, -43, -16, -43, 98};
          denseM<double> m3(3, v2, 1);
          cout << "3x3 matrix m3 is: "
               << "\n"
               << m3 << "\n";

          // Create a 3x3 identity matrix
          denseM<double> m4(3);
          cout << "3x3 identity matrix m4 is: "
               << "\n"
               << m4 << "\n";

          // Create a 3x3 matrix given a file
          denseM<double> f1(3, 3, "matrix.csv");
          cout << "3x3 matrix f1 is: "
               << "\n"
               << f1 << "\n";
          // Output and write the matrix into a file
          f1.output("m2.csv");

          // Change precision of m1 from double to float
          denseM<float> m1_F(m1);
          cout << "m1 in double precision is:"
               << "\n"
               << m1 << "\n";
          cout << "m1 in float precision is: "
               << "\n"
               << m1_F << "\n";

          // 2. Important Member functions
          // Check the value in 1st row, 2nd column in m1 by "at"
          cout << "2x3 matrix m1 is: "
               << "\n"
               << m1 << "\n";
          cout << "m1's value in 1st row, 2nd column is: " << m1.at(1, 2) << "\n";

          // Check the value on index 4 in m1 by "[]"
          cout << "m1's value at index 4 is: " << m1[4] << "\n";
          // Change the value on index 4 in m1 by "[]"
          m1[4] = 5.1234;
          cout << "m1's changed value on index 4 is: " << m1[4] << "\n\n";

          // Restore the matrix
          m1[4] = 5;

          // Transpose
          cout << "The transpose of a 2x3 matrix m1 is: "
               << "\n"
               << transpose(m1) << "\n";

          // Resize
          denseM<double> M1(3, 3, "matrix.csv");
          M1.resize(2, 5);
          cout << "Resize the 3x3 matrix m1 to 2x5 matrix: "
               << "\n"
               << M1 << "\n";

          // 3. Matrix operation
          vector<double> a1 = {1, 2, 3, 4.5};
          vector<double> a2 = {-1, 2.5, -6, -0.4};
          denseM<double> A1(2, 2, a1);
          denseM<double> A2(2, 2, a2);
          cout << "matrix A1: "
               << "\n"
               << A1 << "matrix A2: "
               << "\n"
               << A2 << "\n";

          // Infinity-norm
          cout << "Infinity norm of A1: " << norm_inf(A1) << "\n";
          cout << "Infinity norm of A2: " << norm_inf(A2) << "\n";

          // Inverse
          vector<double> a3 = {2, 4, 6, 3, 5, 1, 6, -2, 2};
          denseM<double> A3(3, 3, a3);
          cout << "Inverse of A3: \n"
               << inverse(A3) << "\n";
          cout << "Test inverse of A3 by A3 * A3^-1 = I: \n"
               << A3 * inverse(A3) << "\n";

          // 4. Decompositions and solvers
          // LU-decomposition
          vector<double> d1 = {5.23, 2.11, 3.15, 0, 1.67, 4.57, 10.111, 6.223, 0};
          denseM<double> D1(3, 3, d1);
          cout << "The original matrix before LU-decomposition is : "
               << "\n"
               << D1 << "\n";
          vector<uint64_t> P1 = {0, 1, 2};
          uint64_t exit_flag = LU_withPivot(D1, P1);
          cout << "After LU-decomposition, combined both L and U: "
               << "\n"
               << D1 << "\n"
               << "The exit_flag is: " << exit_flag
               << "\n\n";

          // LU-solver
          // Need to reset D1
          denseM<double> D1_origin(3, 3, d1);
          vector<double> b1 = {1, 2, 3};
          denseM<double> B1(3, 1, b1);
          denseM<double> x1 = LU_solver(D1_origin, B1, P1);
          cout << "The result from LU solver is: "
               << "\n"
               << x1;

          cout << "Verify the result by norm(b - Ax): "
               << "\n"
               << norm_inf(B1 - (D1_origin * x1)) << "\n\n";

          // Cholesky decomposition
          vector<double> d2 = {4, 12, -16, 12, 37, -43, -16, -43, 98};
          denseM<double> D2(3, d2, 1);
          cout << "The original matrix before Cholesky-decomposition is : "
               << "\n"
               << D2 << "\n";
          denseM<double> D2_chole = cholesky(D2);
          cout << "After Cholesky-decomposition, combined both L and L_T: "
               << "\n"
               << D2_chole << "\n";
          // Cholesky solver
          vector<double> b2 = {-22.2, -47.74, 213.12};
          denseM<double> B2(3, 1, b2);
          denseM<double> x2 = cholesky_solver(D2, B2);
          cout << "The result from Cholesky solver is: "
               << "\n"
               << x2;
          cout << "Verify the result by norm(b - Ax): "
               << "\n"
               << norm_inf(B2 - (D2 * x2)) << "\n\n";

          // 5. Iterative Refinement
          // IR for regular dense matrix
          vector<double> ir_a1 = {5.23, 2.11, 3.15, 0, 1.67, 4.57, 10.111, 6.223, 0};
          denseM<double> IR_A1(3, 3, ir_a1);
          vector<double> ir_b1 = {1, 2, 3};
          denseM<double> IR_B1(3, 1, ir_b1);
          // Three precisions are float, double and double
          denseM<double> result1 = IR<float, double, double>(IR_A1, IR_B1);
          cout << result1 << "\n";

          // IR for positive definite matrix
          vector<double> ir_a2 = {4, 12, -16, 12, 37, -43, -16, -43, 98};
          denseM<double> IR_A2(3, ir_a2, 1);
          vector<double> ir_b2 = {-22.2, -47.74, 213.12};
          denseM<double> IR_B2(3, 1, ir_b2);
          // Three precisions are float, double and double
          denseM<double> result2 = IR<float, double, double>(IR_A2, IR_B2);
          cout << result2 << "\n";

          // GMRES
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

          // GMRES-IR
          denseM<double> GMRES_IR_result = GMRES_IR<float, double, double>(GMRES_A1, GMRES_B1);
          cout << GMRES_IR_result << "\n";
     }

     catch (const denseM<double>::size_mismatch &e)
     {
          cout << "The matrix vector does not have the same amount of elements as the given size.\n";
     }
     catch (const denseM<double>::invalid_size &e)
     {
          cout << "The matrix size can only be positive.\n";
     }
     catch (const denseM<double>::index_overflow &e)
     {
          cout << "The input index is beyond the size of the matrix.\n";
     }
     catch (const invalid_argument &e)
     {
          cout << "Error: " << e.what() << '\n';
     }
}
