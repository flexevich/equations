#include "Equations.h"
#include <iostream>

using namespace std;

Equations::Equations() {      // хуета которая все чистит
    array = nullptr;
    equally = nullptr;
}
Equations::Equations(int n) 
{   
    this->n = n;
    array = new double*[n];
    equally = new double[n];
    for (int i=0;i<n;i++)
    {
        array[i] = new double[n];
    }
}

Equations::~Equations() {
    if (array != nullptr) {
        for (int i = 0; i < n; i++) {
            delete[] array[i];
        }
        delete[] array;
    }

    if (equally != nullptr) {
        delete[] equally;
    }
}

void Equations::wEquations() {           //хуета которая создает массив с СЛАУ
    cout << "Number of equations: " << endl;
    cin >> n;

    array = new double*[n];
    equally = new double[n];

    for (int i = 0; i < n; i++) {
        array[i] = new double[n];
        for (int j = 0; j < n; j++) {
            cout << "Coefficient (" << i << ")" << "(" << j <<") :" << endl;
            cin >> array[i][j];
        }
    }

    for (int i = 0; i < n; i++) {  // массив чему равен СЛАУ
        cout << "Equally (" << i << ") :" << endl;
        cin >> equally[i];
    }
}

void Equations::print() {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << array[i][j] << "*x" << j;
            if (j < n - 1)
                cout << " + ";
        }
        cout << " = " << equally[i] << endl;
    }
}

void Equations::slae_gauss()
{   
    
    //создаём временные массивы, выделяем им память и заполняем соответственными элементами исходных массивов
    double **var_a;
    double *var_e;
    var_a = new double*[n];
    var_e = new double[n];
    if (det()==0)
    {
        cout << "det is equal to zero, the system either has no solutions, or there are an infinite number of them\n";
        return;
    }
    for (int i=0;i<n;i++)
    {   
        var_a[i] = new double[i];
        for (int j=0;j<n;j++)
        {
            var_a[i][j] = array[i][j];
        }
        var_e[i] = equally[i];
    }
    //зануление нижнего левого угла матрицы
    for (int k = 0; k < n; k++)
    {
        if (fabs(array[k][k]) < 0.00000001)
        {   
            for(int i=0;i<n;i++)
            {   
                auto ij = 0;
                for (int j=0;j<n;j++)
                {
                    ij += array[i][j];
                    if (ij == 0&& equally[i]!=0)
                    {
                        cout << "there is no solutions "<<ij<<"!="<<equally[i]<<"\n";
                        return;
                    }
                }
            }
            
            return;
        }
        for (int i = 0; i < n; i++) //нормируем элементы строки k - делим на первый ненулевой элемент строки
            var_a[k][i] = var_a[k][i]/array[k][k];
        var_e[k] = var_e[k]/array[k][k];

        for (int i = k + 1; i < n; i++) // i - следующая строка после строки k
        {
            double var_coeff = var_a[i][k]/var_a[k][k];
            for (int j = 0; j < n; j++)
            {
                var_a[i][j] = var_a[i][j]- var_a[k][j] * var_coeff;
            }
            var_e[i] = var_e[i]- var_coeff * var_e[k];
        }
        
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                array[i][j] = var_a[i][j];
            }
            equally[i] = var_e[i];
        }

    }

    //зануление верхнего нижнего угла матрицы (обратная операция)
    for (int k=n-1;k>=0;k--)
    {
        if (fabs(array[k][k]) < 0.00000001)
        {
            for (int i = 0; i < n; i++)
            {
                auto ij = 0;
                for (int j = 0; j < n; j++)
                {
                    ij += array[i][j];
                    if (ij == 0 && equally[i] != 0)
                    {
                        cout << "there is no solutions " << ij << "!=" << equally[i]<<"\n";
                        return;
                    }
                }
            }
            return;
        }
        for (int i = n - 1; i >= 0; i--)
            var_a[k][i] = var_a[k][i] / array[k][k];
        var_e[k] = var_e[k] / array[k][k];
        
        for (int i=k-1;i>=0;i--)
        {
            double var_coeff = var_a[i][k] / var_a[k][k];
            for (int j = n - 1; j >= 0; j--)
            {
                var_a[i][j] = var_a[i][j] - var_a[k][j] * var_coeff;
            }
            var_e[i] = var_e[i] - var_e[k] * var_coeff;
        }
        
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                array[i][j] = var_a[i][j];
            }
            equally[i] = var_e[i];
        }
        
    }
}
double Equations::minor(int i,int j)
{       
    Equations minor(n-1);
    int a = 0;
    int b = 0;
    for (int x = 0; x < n; x++)
    {   
        if (x != i)
        {
            for (int y = 0; y < n; y++)
            {
                if (y!=j)
                {
                    minor.array[a][b] = array[x][y];
                    b += 1;
                }
            }
            b = 0;
            a += 1;
        }
    }
    return minor.det();
}

double Equations::det()
{
    if (n==2)
    {
        return array[0][0] * array[1][1] - array[0][1] * array[1][0];
    }
    if (n>2)
    {   
        double k = 0;
        int m1 = 1;
        for (int j=0;j<n;j++)
        {  
            k += m1*array[0][j]*minor(0,j);
            m1 *= -1;
            
        }
        return k;
    }
}