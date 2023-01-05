
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