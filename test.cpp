#include "test.hpp"


void test()
{
    int n = 4;
    double ** A = new double*[n];
    for (int i = 0; i < n; i++)
        A[i] = new double[n];

    A[0][0] = 2;
    A[0][1] = 1;
    A[0][2] = 3;
    A[0][3] = 0;

    A[1][0] = 4;
    A[1][1] = 2;
    A[1][2] = 1;
    A[1][3] = 0;

    A[2][0] = 1;
    A[2][1] = 1;
    A[2][2] = 2;
    A[2][3] = 0;

    A[3][0] = 4;
    A[3][1] = 1;
    A[3][2] = 0;
    A[3][3] = 1;

    std::cout << "Matrix A: " << std::endl;
    printMatrix(A,n);
    std::cout << std::endl;

    // Aufruf der LR-Zerlegung
    CLR_Fact lr(A,n);

    std::cout << "LR-Zerlegung von A: " << std::endl;
    lr.printLR();
    std::cout << std::endl;


    //###########################
    //#####     Test 1      #####
    //###########################

    double * b = new double[n];
    double * x = new double[n];
//    for (int i = 0; i < n; i++)
//        b[i] = i+1;
    b[0] = 9;
    b[1] = 3;
    b[2] = 7;
    b[3] = 2;

    std::cout << "Vektor b: " << std::endl;
    printVector(b,n);
    std::cout << std::endl;

    // Lösung des Gleichungssystems A * x = b
    lr.solve(x,b);

    std::cout << "Loesung von A * x = b: " << std::endl;
    printVector(x,n);
    std::cout << std::endl;



    //###########################
    //#####     Test 2      #####
    //###########################

    std::cout << "Determinante von A: " << lr.det() << std::endl << std::endl;



    //###########################
    //#####     Test 3      #####
    //###########################

    double ** A_inv = new double*[n];
    for (int i = 0; i < n; i++)
        A_inv[i] = new double[n];

    // Berechnung der Inversen von A
    lr.inv(A_inv);

    std::cout << "Inverse Matrix von A: " << std::endl;
    printMatrix(A_inv,n);

    for (int i = 0; i < n; i++)
    {
        delete [] A[i];
        A[i] = NULL;
        delete [] A_inv[i];
        A_inv[i] = NULL;
    }
    delete [] A;
    A = NULL;
    delete [] A_inv;
    A_inv = NULL;
    delete [] x;
    x = NULL;
    delete [] b;
    b = NULL;

}

void random_test(int n)
{
    srand((unsigned)time(0));

    double ** A;
    A = new double*[n];
    for (int i = 0; i < n; i++)
    {
        A[i] = new double[n];
    }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
        {
            if (rand() % 2 == 0)
                A[i][j] = (double)(rand()%10);
            else
                A[i][j] = -(double)(rand()%10);
        }

    std::cout << "Zufallsmatrix: " << std::endl;
    printMatrix(A,n);
    std::cout << std::endl;

    CLR_Fact lr(A,n);

    std::cout << "LR-Zerlegung der Zufallsmatrix: " << std::endl;
    lr.printLR();
    std::cout << std::endl;

    double ** A_inv;
    A_inv = new double*[n];
    for (int i = 0; i < n; i++)
    {
        A_inv[i] = new double[n];
    }

    lr.inv(A_inv);

    std::cout << "Inverse Matrix: " << std::endl;
    printMatrix(A_inv,n);

    for (int i = 0; i < n; i++)
    {
        delete [] A[i];
        A[i] = NULL;
        delete [] A_inv[i];
        A_inv[i] = NULL;
    }
    delete [] A;
    A = NULL;
    delete [] A_inv;
    A_inv = NULL;
}
