#include "LR_Fact.hpp"

const double CLR_Fact::s_epsilon = 1e-9;     // Epsilon wird mit 1e-9 initialisiert

CLR_Fact::CLR_Fact(double ** matrix, int dimension)
{
    m_n = dimension;
    m_matrix = matrix;

    m_indizes = new int[m_n];

    for (int i = 0; i < m_n; i++)
        m_indizes[i] = i;

    m_inv = factorize();
}


CLR_Fact::~CLR_Fact()
{
    m_matrix = NULL;

    delete [] m_indizes;
    m_indizes = NULL;
}


void CLR_Fact::add_rows(int i1, int i2, double scalar)
{
    for (int j = i2; j < m_n; j++)  // Damit die Einträge unterhalb der Diagonalen (die ja die Einträge von L speichern) nicht verändert werden, läuft die Schleife erst bei i2 los.
        if (std::abs(m_matrix[i1][j] + scalar * m_matrix[i2][j]) <= s_epsilon)
            m_matrix[i1][j] = 0.0;
        else
            m_matrix[i1][j] += scalar * m_matrix[i2][j];
}


double CLR_Fact::det()
{
    if (m_inv)
    {
        double det = 1.0;
        for (int i = 0; i < m_n; i++)   // Multiplikation der Diagonalelemente der Matrix (= Diagonaleinträge von R) liefert die Determinante
            det *= m_matrix[i][i];
        return det;
    }
    else
        return 0.0;
}


void CLR_Fact::inv(double ** inv_mat)
{
    if (m_inv)    // Check, ob die Matrix überhaupt invertierbar ist
    {
        double * e = new double[m_n];

        for (int i = 0; i < m_n; i++)
            e[i]  = 0.0;

        // Lösen des Gleichungssystems A * x = e_i für alle Einheitsvektoren e_i
        for (int i = 0; i < m_n; i++)
        {
            e[i] = 1.0;
            solve(inv_mat[i], e);
            e[i] = 0.0;
        }

        // Matrix muss noch transponiert werden
        double temp;
        for (int i = 0; i < m_n; i++)
            for (int j = i; j < m_n; j++)
            {
                temp = inv_mat[j][i];
                inv_mat[j][i] = inv_mat[i][j];
                inv_mat[i][j] = temp;
            }

        delete [] e;
        e = NULL;
    }
    else
        std::cout << "Matrix ist nicht invertierbar!" << std::endl;
}


void CLR_Fact::solve(double * x, double * b)
{
    if (m_inv)
    {
        double * temp = new double[m_n];

        // Die Einträge im Vektor b müssen erst umgeordnet werden, da bei der LR-Zerlegung ggfs. Zeilen vertauscht wurden
        for (int i = 0; i < m_n; i++)
            temp[i] = b[i];
        for (int i = 0; i < m_n; i++)
            b[i] = temp[m_indizes[i]];

        delete [] temp;
        temp = NULL;

        // Lösen des Gleichungssystems  L * y = b (mit y = R * x)
        // = Vorwärtseinsetzen
        double * y = new double [m_n];
        for (int i = 0; i < m_n; i++)
        {
            y[i] = b[i];
            for (int  j = 0; j < i; j++)
                y[i] -= m_matrix[i][j] * y[j];
        }

        // Lösen des Gleichungssystems  R * x = y
        // = Rückwärtseinsetzen
        double a;
        for (int j = m_n-1; j >= 0; j--)
        {
            a = 0.0;
            for (int k = j+1; k < m_n; k++)
                a += m_matrix[j][k] * x[k];
            x[j] = 1/m_matrix[j][j] * (y[j] - a);
        }
        delete [] y;
        y = NULL;
    }
    else
        std::cout << "Matrix ist nicht invertierbar! Gleichungssystem nicht eindeutig loesbar!" << std::endl;
}


/*  Die factorize-Methode führt die LR-Zerlegung für eine Matrix aus.
    Die Methode wird im Konstruktor der Klasse mit aufgerufen und ist außerhalb der Klasse NICHT aufrufbar.
    Das garantiert, dass für die weiterführenden Methoden(det, inv, solve) immer schon eine gültige LR-Zerlegung vorhanden ist (falls man die Matrix A außerhalb der Klasse nicht noch irgendwie verändert).
    Die LR-Zerlegung wird über die Matrix A gespeichert (insbes. gehen die Einträge von A dabei verloren!).
    Die Matrix R sind dabei die Einträge oberhalb der Diagonalen, die Matrix L die Einträge unterhalb. Für die Diagonale der Matrix L (die ja nur aus Einsen besteht) werden die Einträge nicht extra gespeichert.
*/
bool CLR_Fact::factorize()
{
    double max_elem, scalar;
    int k, temp2;
    double * temp1;

    for (int i = 0; i < m_n; i++)   // Die Matrix wird Spaltenweise durchlaufen
    {
        max_elem = std::abs(m_matrix[i][i]); // Pivotelement wird zunächst mit dem Diagonalelement der aktuellen Spalte initialisiert
        k = i;                                  // Speichert die Zeile, in der das Pivotelement gefunden wird
        for (int j = i; j < m_n; j++)
            if (std::abs(m_matrix[j][i]) > max_elem)    // Falls in der aktuellen Spalte ein betragsmäßig größeres Element gefunden wird, wird dieses als Pivotelement gesetzt
            {
                max_elem = std::abs(m_matrix[j][i]);
                k = j;
            }
        if (max_elem <= s_epsilon)  // Falls in der aktuellen Spalte kein Pivotelement gefunden werden kann, so ist A nicht invertierbar.
            return false;

        if (k != i)    // Vertauschen der k-ten mit der i-ten Zeile
        {
            // Die Einträge der Zeilen müssen nicht kopiert werden, stattdessen kopiert man einfach die Zeiger auf die i-te und k-te Zeile
            temp1 = m_matrix[i];
            m_matrix[i] = m_matrix[k];
            m_matrix[k] = temp1;

            // Die Vertauschung wird im Indexvektor gespeichert
            temp2 = m_indizes[i];
            m_indizes[i] = m_indizes[k];
            m_indizes[k] = temp2;
        }

        for (int j = i+1; j < m_n; j++)     // Elimination der weiter unten stehenden Einträge der aktuellen Spalte
        {
            scalar = m_matrix[j][i]/m_matrix[i][i];
            add_rows(j,i,-scalar);
            m_matrix[j][i] = scalar;    // Der Skalar, der zur Elimination des Eintrags einer Spalte benutzt wurde, wird an genau diese Stelle gespeichert -> Matrix L
        }
    }

    return true;    // Wenn die Faktorisierung erfolgreich war, wird true zurückgegeben, da die Matrix A dann invertierbar ist
}

void CLR_Fact::printLR()
{
    if (m_inv)
    {
        int l = 7;
        std::cout << "Matrix L: " << std::endl;
        for (int i = 0; i < m_n; i++)
        {
            for (int j = 0; j < m_n; j++)
                if (i > j)
                    std::cout << std::setw(l) << m_matrix[i][j];
                else if (i < j)
                    std::cout << std::setw(l) << 0.0;
                else
                    std::cout << std::setw(l) << 1.0;
            std::cout << std::endl;
        }

        std::cout << std::endl << "Matrix R: " << std::endl;
        for (int i = 0; i < m_n; i++)
        {
            for (int j = 0; j < m_n; j++)
                if (i <= j)
                    std::cout << std::setw(l) << m_matrix[i][j];
                else if (i > j)
                    std::cout << std::setw(l) << 0.0;
            std::cout << std::endl;
        }
    }
    else
        std::cout << "Matrix hat keine LR-Zerlegung!" << std::endl;
}
