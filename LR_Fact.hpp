#ifndef LR_FACT_H
#define LR_FACT_H

#include <iostream>
#include <iomanip>
#include <cmath>

class CLR_Fact
{
private:

    int m_n;                        // Dimension der Matrix/des Vektors
    double ** m_matrix;             // Matrix
    int * m_indizes;                // Indexvektor, der die Vertauschungen w�hrend der LR-Zerlegung speichert
    bool m_inv;                     // Variable, die speichert, ob die Matrix invertierbar ist
    static const double s_epsilon;  // Statisches Datenelement zum Testen auf 0

    CLR_Fact();         // Default-Konstruktor soll nicht aufgerufen werden k�nnen -> private
    bool factorize();   // Methode, die die LR-Zerlegung durchf�hrt. Diese wird beim Konstruktoraufruf mit aufgerufen. Ein Aufruf von au�erhalb der Klasse macht keinen Sinn -> private


public:

    CLR_Fact(double ** matrix, int dimension);  // Konstruktor, dem eine Matrix (inkl. Dimension) �bergeben wird, f�r die die LR-Zerlegung durchgef�hrt wird
    ~CLR_Fact();    // Destruktor

    double det();                           // Methode zur Berechnung der Determinante
    void inv(double ** inv_mat);            // Methode zur Berechnung der inversen Matrix
    void solve(double * x, double * b);     // Methode zum L�sen des Gleichungssystems A * x = b

    void add_rows(int i1, int i2, double scalar);   // Hilfsmethode, die das scalar-fache der i2-ten Zeile zur i1-ten Zeile addiert
    void printLR();     // Ausgabemethode, die die Matrizen L und R ausgibt
};
#endif
