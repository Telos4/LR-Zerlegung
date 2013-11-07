#include "MatrixVector.hpp"

void printMatrix(double ** A, int n)
{
    int l = 7;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            std::cout << std::setw(l) << A[i][j] << "  ";
        }
        std::cout << std::endl;
    }
}


void printVector(double * vec, int n)
{
    int l = 7;
    for (int i = 0; i < n; i++)
    {
        std::cout << std::setw(l) << vec[i] << std::endl;
    }
}
