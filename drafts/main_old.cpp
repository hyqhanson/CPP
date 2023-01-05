#include <iostream>
#include <stdexcept>
#include <vector>
#include <string>
#include <algorithm>
#include <random>
#include "denseM_old.hpp"
using namespace std;

int main()
{
     try
     {
          cout.precision(10);
          // Constructors
          // Create random 2x3 matrix
          denseM<uint64_t, double> d1(2, 3);
          cout << d1 << "\n";
          // Create random 3x2 matrix
          denseM<uint64_t, double> d2(3, 2);
          cout << d2 << "\n";
          cout << d1 * d2 << "\n";
          // Change precision of d2 from double to float
          denseM<uint32_t, float> d3(d2);
          cout << d3 << "\n";
          // create 3x3 matrix with diagonal {1,2,3}
          vector<double> v1 = {1, 2, 3};
          denseM<uint32_t, double> d4(v1);
          cout << d4 << "\n";

          // Given the vector and size, create specific matrix
          vector<double> m = {2.231, 5.42, 1.45, -2.2111, 3.33, -2.3, 0, 7.2, 9};
          denseM<uint32_t, double> A(3, 3, m);
          cout << "Inf norm is: " << norm<uint32_t, double>(A) << "\n";
          cout << A << "\n";
          // LU decomposition
          denseM<uint32_t, double> P;
          denseM<uint32_t, double> U;
          denseM<uint32_t, double> L;
          LU_withPivot(A, P, L, U);
          cout << A << "\n";
          cout << "Permutation matrix: "
               << "\n"
               << P << "\n";
          cout << "Lower-triangular matrix: "
               << "\n"
               << L << "\n";
          cout << "Upper-triangular matrix: "
               << "\n"
               << U << "\n";
          // cout << P * A << "\n";
          // cout << L * U << "\n";
          vector<double> v2 = {1, 2, 3};
          denseM<uint32_t, double> b(3, 1, v2);

          // Solve by LU-decomposition
          denseM<uint32_t, double> x = LU_solver(L, U, P, b);
          cout << x << "\n";

          // Solve by IR
          denseM<uint32_t, double> x2 = IR(A, b);
          cout << "x2 is :" << x2 << "\n";
     }
     catch (const invalid_argument &e)
     {
          cout << "Error: " << e.what() << '\n';
     }
}