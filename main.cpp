#include <iostream>
#include <chrono>
#include <stdexcept>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <cstring>
#include "denseM.hpp"
using namespace std;

int main(int argc, char *argv[])
{
     try
     {
          // If running the program without putting any commend line argument,
          // print error message.
          if (argc == 1)
          {
               cout << "You need to put some arguments in command line.";
               return 1;
          }

          // Using LU-factorization to solve linear system
          // Example: main LU_solver double 3 LU1.csv LU2.csv LU_result.csv
          // Example meaning: solve linear system by LU in double precision, using a 3x3 matrix
          // from LU1.csv and a 3x1 matrix from LU2.csv, store the result in LU_result.csv
          if (strcmp(argv[1], "LU_solver") == 0)
          {
               // Check if the number of command line argument is correct.
               if (argc != 7)
               {
                    cout << "The number of arguments is incorrect. "
                         << "You need 7 arguments in total. \n"
                         << "Example: main LU_solver double 3 LU1.csv LU2.csv LU_result.csv";
                    return 1;
               }

               // Check if the input size a valid positive integer
               // it can also prevent negative integers
               for (uint64_t i = 0; i < strlen(argv[3]); i++)
               {
                    if (isdigit(argv[3][i]) == 0)
                    {
                         cout << "The size of the matrix needs be an positive integer!";
                         return 1;
                    }
               }

               uint64_t row = (uint64_t)(stoi(argv[3]));
               uint64_t col = (uint64_t)(stoi(argv[3]));
               // If user want the precision be float
               if (strcmp(argv[2], "single") == 0)
               {
                    // A is a nxn matrix
                    denseM<float> A(row, col, argv[4]);
                    // b is a nx1 matrix
                    denseM<float> b(row, 1, argv[5]);
                    denseM<float> x = LU_solver(A, b);
                    // If the last command line argument is "N", don't output a file but print
                    // the result on the terminal
                    if (strcmp(argv[6], "N") == 0)
                    {
                         cout << "The solution of the linear system by LU decomposition in single precision is: \n"
                              << x << "\n";
                         cout << "Verify the result by norm(b - Ax): "
                              << "\n"
                              << norm_inf(b - (A * x)) << "\n\n";
                         return 0;
                    }
                    // Otherwise, use the last command line argument as the file name,
                    // stores the result
                    else
                    {
                         x.output(argv[6]);
                         return 0;
                    }
               }
               // If user want the precision be double
               else if (strcmp(argv[2], "double") == 0)
               {
                    denseM<double> A(row, col, argv[4]);
                    denseM<double> b(row, 1, argv[5]);
                    denseM<double> x = LU_solver(A, b);
                    if (strcmp(argv[6], "N") == 0)
                    {
                         cout << "The solution of the linear system by LU decomposition in double precision is: \n"
                              << x << "\n";
                         cout << "Verify the result by norm(b - Ax): "
                              << "\n"
                              << norm_inf(b - (A * x)) << "\n\n";
                         return 0;
                    }
                    else
                    {
                         x.output(argv[6]);
                         return 0;
                    }
               }
               else
               {
                    cout << "Please use either single or double precision in the 3rd argument.";
               }
          }

          // Using Cholesky-factorization to solve linear system
          // Example: main cholesky_solver double 3 C1.csv C2.csv C_result.csv
          // Example meaning: solve linear system by cholesky in double precision, using a 3x3 matrix
          // from C1.csv and a 3x1 matrix from C2.csv, store the result in C_result.csv
          else if (strcmp(argv[1], "cholesky_solver") == 0)
          {
               // Check if the number of command line argument is correct.
               if (argc != 7)
               {
                    cout << "The number of arguments is incorrect. "
                         << "You need 7 arguments in total. \n"
                         << "Example: main cholesky_solver double 3 C1.csv C2.csv C_result.csv";
                    return 1;
               }

               // Check if the input size a valid positive integer
               // it can also prevent negative integers
               for (uint64_t i = 0; i < strlen(argv[3]); i++)
               {
                    if (isdigit(argv[3][i]) == 0)
                    {
                         cout << "The size of the matrix needs be an positive integer!";
                         return 1;
                    }
               }

               uint64_t size = (uint64_t)(stoi(argv[3]));
               // If user want the precision be float
               if (strcmp(argv[2], "single") == 0)
               {
                    // A is a nxn matrix
                    denseM<float> A(size, size, argv[4]);
                    A.change_is_pos_def();
                    // b is a nx1 matrix
                    denseM<float> b(size, 1, argv[5]);
                    denseM<float> x = cholesky_solver(A, b);
                    // If the last command line argument is "N", don't output a file but print
                    // the result on the terminal
                    if (strcmp(argv[6], "N") == 0)
                    {
                         cout << "The solution of the linear system by Cholesky decomposition in single precision is: \n"
                              << x << "\n";
                         cout << "Verify the result by norm(b - Ax): "
                              << "\n"
                              << norm_inf(b - (A * x)) << "\n\n";
                         return 0;
                    }
                    // Otherwise, use the last command line argument as the file name,
                    // stores the result
                    else
                    {
                         x.output(argv[6]);
                         return 0;
                    }
               }
               // If user want the precision be double
               else if (strcmp(argv[2], "double") == 0)
               {
                    denseM<double> A(size, size, argv[4]);
                    A.change_is_pos_def();
                    denseM<double> b(size, 1, argv[5]);
                    denseM<double> x = cholesky_solver(A, b);
                    if (strcmp(argv[6], "N") == 0)
                    {
                         cout << "The solution of the linear system by Cholesky decomposition in double precision is: \n"
                              << x << "\n";
                         cout << "Verify the result by norm(b - Ax): "
                              << "\n"
                              << norm_inf(b - (A * x)) << "\n\n";
                         return 0;
                    }
                    else
                    {
                         x.output(argv[6]);
                         return 0;
                    }
               }
               else
               {
                    cout << "Please use either single or double precision in the 3rd argument.";
               }
          }

          // Using GMRES with initial guess to solve linear system
          // Example: main GMRES double 3 G1.csv G2.csv guess.csv G_result.csv
          // Example meaning: solve linear system by GMRES in double precision, using a 3x3 matrix
          // from G1.csv, a 3x1 matrix from G2.csv, and a 3x1 initial guess. Store the result in G_result.csv
          else if (strcmp(argv[1], "GMRES") == 0)
          {
               // Check if the number of command line argument is correct.
               if (argc != 8)
               {
                    cout << "The number of arguments is incorrect. "
                         << "You need 8 arguments in total. \n"
                         << "Example: main GMRES double 3 G1.csv G2.csv guess.csv G_result.csv";
                    return 1;
               }

               // Check if the input size a valid positive integer
               // it can also prevent negative integers
               for (uint64_t i = 0; i < strlen(argv[3]); i++)
               {
                    if (isdigit(argv[3][i]) == 0)
                    {
                         cout << "The size of the matrix needs be an positive integer!";
                         return 1;
                    }
               }

               uint64_t row = (uint64_t)(stoi(argv[3]));
               uint64_t col = (uint64_t)(stoi(argv[3]));
               // If user want the precision be float
               if (strcmp(argv[2], "single") == 0)
               {
                    // A is a nxn matrix
                    denseM<float> A(row, col, argv[4]);
                    // b is a nx1 matrix
                    denseM<float> b(row, 1, argv[5]);

                    // If the 6th argument is "N", it means not giving initial guess, assume
                    // it is a nx1 zero matrix
                    denseM<float> init_guess(row, 1);
                    // Else, take the matrix from the file at the 6th argument as initial guess
                    if (strcmp(argv[6], "N") != 0)
                    {
                         denseM<float> new_guess(row, 1, argv[6]);
                         init_guess = new_guess;
                    }
                    denseM<float> x = GMRES(A, b, init_guess);
                    // If the last command line argument is "N", don't output a file but print
                    // the result on the terminal
                    if (strcmp(argv[7], "N") == 0)
                    {
                         cout << "The solution of the linear system by GMRES in single precision is: \n"
                              << x << "\n";
                         cout << "Verify the result by norm(b - Ax): "
                              << "\n"
                              << norm_inf(b - (A * x)) << "\n\n";
                         return 0;
                    }
                    // Otherwise, use the last command line argument as the file name,
                    // stores the result
                    else
                    {
                         x.output(argv[7]);
                         return 0;
                    }
               }
               // If user want the precision be double
               else if (strcmp(argv[2], "double") == 0)
               {
                    denseM<double> A(row, col, argv[4]);
                    denseM<double> b(row, 1, argv[5]);

                    // If the 6th argument is "N", it means not giving initial guess, assume
                    // it is a nx1 zero matrix
                    denseM<double> init_guess(row, 1);
                    // Else, take the matrix from the file at the 6th argument as initial guess
                    if (strcmp(argv[6], "N") != 0)
                    {
                         denseM<double> new_guess(row, 1, argv[6]);
                         init_guess = new_guess;
                    }
                    denseM<double> x = GMRES(A, b, init_guess);

                    if (strcmp(argv[7], "N") == 0)
                    {
                         cout << "The solution of the linear system by GMRES in double precision is: \n"
                              << x << "\n";
                         cout << "Verify the result by norm(b - Ax): "
                              << "\n"
                              << norm_inf(b - (A * x)) << "\n\n";
                         return 0;
                    }
                    else
                    {
                         x.output(argv[7]);
                         return 0;
                    }
               }
               else
               {
                    cout << "Please use either single or double precision in the 3rd argument.";
               }
          }

          // Using iterative refinement to solve linear system
          // Example: main IR single double double 3 IR_LU1.txt IR_LU2.txt IR_result.txt LU
          // Example meaning: solve linear system by Iterative Refinement in mixed precision,
          // (single, double, double). Using a 3x3 matrix from IR_LU1.txt and a 3x1 matrix from
          // IR_LU2.txt, store the result in IR_result.txt. Using LU factorization.
          else if (strcmp(argv[1], "IR") == 0)
          {
               // Check if the number of command line argument is correct.
               if (argc != 10)
               {
                    cout << "The number of arguments is incorrect. "
                         << "You need 10 arguments in total. \n"
                         << "Example: main IR single double double 3 IR_LU1.txt IR_LU2.txt IR_result.txt LU";
                    return 1;
               }

               // Check if the input size a valid positive integer
               // it can also prevent negative integers
               for (uint64_t i = 0; i < strlen(argv[5]); i++)
               {
                    if (isdigit(argv[5][i]) == 0)
                    {
                         cout << "The size of the matrix needs be an positive integer!";
                         return 1;
                    }
               }
               uint64_t row = (uint64_t)(stoi(argv[5]));
               uint64_t col = (uint64_t)(stoi(argv[5]));

               // Using LU factorization in iterative refinement
               // If the second precision is single
               if (strcmp(argv[3], "single") == 0)
               {
                    denseM<float> A(row, col, argv[6]);
                    // If the user knows the matrix A is positive definite, the user can
                    // choose Cholesky factorization in IR
                    if (strcmp(argv[9], "C") == 0)
                    {
                         A.change_is_pos_def();
                    }
                    // If the user wants IR in LU factorization, nothing changes
                    else if (strcmp(argv[9], "LU") == 0)
                    {
                    }
                    else
                    {
                         cout << "You need to choose either 'LU' or 'C' in the last argument";
                         return 1;
                    }

                    denseM<float> b(row, 1, argv[7]);
                    denseM<float> x(row, 1);
                    // List all combinations
                    if (strcmp(argv[2], "single") == 0 && strcmp(argv[4], "single") == 0)
                    {
                         x = IR<float, float, float>(A, b);
                    }
                    else if (strcmp(argv[2], "single") == 0 && strcmp(argv[4], "double") == 0)
                    {
                         x = IR<float, float, double>(A, b);
                    }
                    else if (strcmp(argv[2], "double") == 0 && strcmp(argv[4], "single") == 0)
                    {
                         x = IR<double, float, float>(A, b);
                    }
                    else if (strcmp(argv[2], "double") == 0 && strcmp(argv[4], "double") == 0)
                    {
                         x = IR<double, float, double>(A, b);
                    }
                    else
                    {
                         cout << "Please use either single or double precision in the 3rd/4th/5th argument.";
                         return 1;
                    }
                    // If user does not want to output a file
                    if (strcmp(argv[8], "N") == 0)
                    {
                         cout << "The solution of the linear system by IR is: \n"
                              << x << "\n";
                         cout << "Verify the result by norm(b - Ax): "
                              << "\n"
                              << norm_inf(b - (A * x)) << "\n\n";
                         return 0;
                    }
                    else
                    {
                         x.output(argv[8]);
                         return 0;
                    }
               }
               // If the second precision is double
               else if (strcmp(argv[3], "double") == 0)
               {
                    denseM<double> A(row, col, argv[6]);

                    // If the user knows the matrix A is positive definite, the user can
                    // choose Cholesky factorization in IR
                    if (strcmp(argv[9], "C") == 0)
                    {
                         A.change_is_pos_def();
                    }
                    // If the user wants IR in LU factorization, nothing changes
                    else if (strcmp(argv[9], "LU") == 0)
                    {
                    }
                    else
                    {
                         cout << "You need to choose either 'LU' or 'C' in the last argument";
                         return 1;
                    }

                    denseM<double> b(row, 1, argv[7]);
                    denseM<double> x(row, 1);
                    // List all combinations
                    if (strcmp(argv[2], "single") == 0 && strcmp(argv[4], "single") == 0)
                    {
                         x = IR<float, double, float>(A, b);
                    }
                    else if (strcmp(argv[2], "single") == 0 && strcmp(argv[4], "double") == 0)
                    {
                         x = IR<float, double, double>(A, b);
                    }
                    else if (strcmp(argv[2], "double") == 0 && strcmp(argv[4], "single") == 0)
                    {
                         x = IR<double, double, float>(A, b);
                    }
                    else if (strcmp(argv[2], "double") == 0 && strcmp(argv[4], "double") == 0)
                    {
                         x = IR<double, double, double>(A, b);
                    }
                    else
                    {
                         cout << "Please use either single or double precision in the 3rd/4th/5th argument.";
                         return 1;
                    }
                    // If user does not want to output a file
                    if (strcmp(argv[8], "N") == 0)
                    {
                         cout << "The solution of the linear system by IR is: \n"
                              << x << "\n";
                         cout << "Verify the result by norm(b - Ax): "
                              << "\n"
                              << norm_inf(b - (A * x)) << "\n\n";
                         return 0;
                    }
                    else
                    {
                         x.output(argv[8]);
                         return 0;
                    }
               }
               else
               {
                    cout << "Please use either single or double precision in the 3rd/4th/5th argument.";
                    return 1;
               }
          }

          // Using preconditioned system GMRES iterative refinement to solve linear system
          // Example: main GMRES_IR single double double 3 G_IR1.csv G_IR2.csv G_IR_Result.csv
          // Example meaning: solve linear system by Iterative Refinement in mixed precision,
          // (single, double, double). Using a 3x3 matrix from G_IR1.csv and a 3x1 matrix from
          // G_IR2.csv, store the result in G_IR_Result.csv.
          else if (strcmp(argv[1], "GMRES_IR") == 0)
          {
               // Check if the number of command line argument is correct.
               if (argc != 9)
               {
                    cout << "The number of arguments is incorrect. "
                         << "You need 9 arguments in total. \n"
                         << "Example: main GMRES_IR single double double 3 G_IR1.csv G_IR2.csv G_IR_Result.csv";
                    return 1;
               }

               // Check if the input size a valid positive integer
               // it can also prevent negative integers
               for (uint64_t i = 0; i < strlen(argv[5]); i++)
               {
                    if (isdigit(argv[5][i]) == 0)
                    {
                         cout << "The size of the matrix needs be an positive integer!";
                         return 1;
                    }
               }
               uint64_t row = (uint64_t)(stoi(argv[5]));
               uint64_t col = (uint64_t)(stoi(argv[5]));

               // Using LU factorization in iterative refinement
               // If the second precision is single
               if (strcmp(argv[3], "single") == 0)
               {
                    denseM<float> A(row, col, argv[6]);
                    denseM<float> b(row, 1, argv[7]);
                    denseM<float> x(row, 1);
                    // List out all combinations
                    if (strcmp(argv[2], "single") == 0 && strcmp(argv[4], "single") == 0)
                    {
                         x = GMRES_IR<float, float, float>(A, b);
                    }
                    else if (strcmp(argv[2], "single") == 0 && strcmp(argv[4], "double") == 0)
                    {
                         x = GMRES_IR<float, float, double>(A, b);
                    }
                    else if (strcmp(argv[2], "double") == 0 && strcmp(argv[4], "single") == 0)
                    {
                         x = GMRES_IR<double, float, float>(A, b);
                    }
                    else if (strcmp(argv[2], "double") == 0 && strcmp(argv[4], "double") == 0)
                    {
                         x = GMRES_IR<double, float, double>(A, b);
                    }
                    else
                    {
                         cout << "Please use either single or double precision in the 3rd/4th/5th argument.";
                         return 1;
                    }

                    // If user does not want to output a file
                    if (strcmp(argv[8], "N") == 0)
                    {
                         cout << "The solution of the linear system by GMRES_IR is: \n"
                              << x << "\n";
                         cout << "Verify the result by norm(b - Ax): "
                              << "\n"
                              << norm_inf(b - (A * x)) << "\n\n";
                         return 0;
                    }
                    else
                    {
                         x.output(argv[8]);
                         return 0;
                    }
               }
               // If the second precision is double
               else if (strcmp(argv[3], "double") == 0)
               {
                    denseM<double> A(row, col, argv[6]);
                    denseM<double> b(row, 1, argv[7]);
                    denseM<double> x(row, 1);
                    // List out all combinations
                    if (strcmp(argv[2], "single") == 0 && strcmp(argv[4], "single") == 0)
                    {
                         x = GMRES_IR<float, double, float>(A, b);
                    }
                    else if (strcmp(argv[2], "single") == 0 && strcmp(argv[4], "double") == 0)
                    {
                         x = GMRES_IR<float, double, double>(A, b);
                    }
                    else if (strcmp(argv[2], "double") == 0 && strcmp(argv[4], "single") == 0)
                    {
                         x = GMRES_IR<double, double, float>(A, b);
                    }
                    else if (strcmp(argv[2], "double") == 0 && strcmp(argv[4], "double") == 0)
                    {
                         x = GMRES_IR<double, double, double>(A, b);
                    }
                    else
                    {
                         cout << "Please use either single or double precision in the 3rd/4th/5th argument.";
                         return 1;
                    }

                    // If user does not want to output a file
                    if (strcmp(argv[8], "N") == 0)
                    {
                         cout << "The solution of the linear system by GMRES_IR is: \n"
                              << x << "\n";
                         cout << "Verify the result by norm(b - Ax): "
                              << "\n"
                              << norm_inf(b - (A * x)) << "\n\n";
                         return 0;
                    }
                    else
                    {
                         x.output(argv[8]);
                         return 0;
                    }
               }
               else
               {
                    cout << "Please use either single or double precision in the 3rd/4th/5th argument.";
                    return 1;
               }
          }

          else
          {
               cout << "The second argument should be one of the functions below: \n"
                    << "LU_solver \ncholesky_solver \nGMRES \nIR \nGMRES_IR";
               return 1;
          }
     }
     catch (const denseM<double>::size_mismatch &e)
     {
          cout << "The matrix vector does not have the same amount of elements as the given size.\n";
     }
     catch (const denseM<double>::invalid_size &e)
     {
          cout << "The matrix size can not be zero!\n";
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
