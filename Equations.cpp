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

Equations::Equations(const Equations& e)
{   
    n = e.n;
    array = new double* [n];
    equally = new double[n];
    for (int i = 0; i < n; i++)
    {
        array[i] = new double[n];
        equally[i] = e.equally[i];
        for (int j = 0; j < n; j++)
            array[i][j] = e.array[i][j];
    }
};
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
    Equations vareq = *this;
    if (det()==0)
    {
        cout << "Determinant is equal to zero, the system either has no solutions, or there are an infinite number of them.\nTry another method.\n"; //тут я думаю всё понятно
        return;
    }

    for (int oof=0;oof<n;oof++)//проверка чтобы бим бим бам бам по главной диагонали не был равен нулю, а если равен прибавим к нему другую строку где он не будет равен нулю
    {
        if (vareq.array[oof][oof] == 0)
        {
            for (int oofi=0;oofi<n;oofi++)
            {
                if(vareq.array[oofi][oof]!=0)
                {
                    for(int oofj=0;oofj<n;oofj++)
                    {
                        vareq.array[oof][oofj] += vareq.array[oofi][oofj];
                    }
                    vareq.equally[oof] += vareq.equally[oofi];
                    break;
                }
            }
        }
    }
    for (int u = n - 1; u > 0; u--)//перемещаем строки чтобы исключить деление на ноль и всякую беду ящеров
    {
        if (vareq.array[u - 1][0] < vareq.array[u][0])
        {
            for (int v = 0; v < n; v++)
            {
                double temp = vareq.array[u][v];
                vareq.array[u][v] = vareq.array[u - 1][v];
                vareq.array[u - 1][v] = temp;
            }
            double temp = vareq.equally[u];
            vareq.equally[u] = vareq.equally[u - 1];
            vareq.equally[u - 1] = temp;
        }
    }
    Equations var_one = vareq;
    //зануление нижнего левого угла матрицы
    for (int k = 0; k < n; k++)
    {
        
        for (int i = 0; i < n; i++) //нормируем элементы строки k - делим на первый ненулевой элемент строки
            var_one.array[k][i] = var_one.array[k][i]/ vareq.array[k][k];
        var_one.equally[k] = var_one.equally[k]/ vareq.array[k][k];

        for (int i = k + 1; i < n; i++) // i - следующая строка после строки k
        {
            double var_coeff = var_one.array[i][k]/ var_one.array[k][k];
            for (int j = 0; j < n; j++)
            {
                var_one.array[i][j] = var_one.array[i][j]- var_one.array[k][j] * var_coeff;
            }
            var_one.equally[i] = var_one.equally[i]- var_coeff * var_one.equally[k];
        }
        
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                vareq.array[i][j] = var_one.array[i][j];
            }
            vareq.equally[i] = var_one.equally[i];
        }

    }

    //зануление верхнего нижнего угла матрицы (обратная операция)
    for (int k=n-1;k>=0;k--)
    {
        for (int i = n - 1; i >= 0; i--)
            var_one.array[k][i] = var_one.array[k][i] / vareq.array[k][k];
        var_one.equally[k] = var_one.equally[k] / vareq.array[k][k];
        
        for (int i=k-1;i>=0;i--)
        {
            double var_coeff = var_one.array[i][k] / var_one.array[k][k];
            for (int j = n - 1; j >= 0; j--)
            {
                var_one.array[i][j] = var_one.array[i][j] - var_one.array[k][j] * var_coeff;
            }
            var_one.equally[i] = var_one.equally[i] - var_one.equally[k] * var_coeff;
        }
        
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                vareq.array[i][j] = var_one.array[i][j];
            }
            vareq.equally[i] = var_one.equally[i];
        }
        
    }
    cout << "Solution:\n";
    for (int i=0;i<n;i++)
    {
        cout << "x" << i+1 << " = " << vareq.equally[i] << "\n";
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