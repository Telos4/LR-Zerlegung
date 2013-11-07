#include <iostream>

#include "LR_Fact.hpp"
#include "test.hpp"
#include "MatrixVector.hpp"

using namespace std;

int main()
{
    double ** mat = new double*[2];
    mat[0] = new double[2];
    mat[1] = new double[2];

    mat[0][0] = 1.0;
    mat[0][1] = 0.99;
    mat[1][0] = 0.99;
    mat[1][1] = 0.98;

    double * b = new double[2];
    b[0] = 1.989903;
    b[1] = 1.970106;

    CLR_Fact lr(mat,2);
    lr.printLR();

    double * x = new double[2];
    lr.solve(x,b);

    printVector(x,2);

    double ** mat1 = new double*[2];
    mat1[0] = new double[2];
    mat1[1] = new double[2];

    lr.inv(mat1);
    printMatrix(mat1,2);

    test();
    /**

    printMatrix(mat,2);
    test();

    cout << endl << endl;

    // Alternativer Test mit Zufallsgenerierter Matrix
    random_test(3);
    */

    return 0;
}
