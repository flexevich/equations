#ifndef EQUATIONS_H
#define EQUATIONS_H

class Equations {
private:
    double** array; // массив коэффициентов
    double* equally; // массив чисел которые идут после =
    int n; // кол-во уравнений

public:
    Equations();
    ~Equations();
    void wEquations();
    void print();
};
#endif