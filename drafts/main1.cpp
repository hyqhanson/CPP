#include <iostream>
#include <vector>

using namespace std;

// Function to perform forward substitution
vector<double> forwardSubstitution(vector<vector<double>> &U, vector<double> &b)
{
    int n = b.size();
    vector<double> x(n, 0);

    // perform forward substitution
    for (int i = 0; i < n; i++)
    {
        double sum = 0;
        for (int j = 0; j < i; j++)
        {
            sum += U[i][j] * x[j];
        }

        x[i] = (b[i] - sum) / U[i][i];
    }

    return x;
}

int main()
{
    // test the forward substitution function
    vector<vector<double>> U = {{1, 2, 3},
                                {0, 4, 5},
                                {0, 0, 6}};
    vector<double> b = {7, 8, 9};

    vector<double> x = forwardSubstitution(U, b);

    // print the solution
    for (double xi : x)
    {
        cout << xi << " ";
    }

    return 0;
}
