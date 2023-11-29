#include "Equations.h"
#include <iostream>

using namespace std;

Equations::Equations() {      // хуета которая все чистит
    array = nullptr;
    equally = nullptr;
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