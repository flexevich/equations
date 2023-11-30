#ifndef EQUATIONS_H
#define EQUATIONS_H

class Equations {
private:
    double** array; // массив коэффициентов
    double* equally; // массив чисел которые идут после =
    int n; // кол-во уравнений

public:
    Equations();
    Equations(int n);
    Equations(const Equations& e);
    ~Equations();
    void wEquations();
    void print();
    void slae_gauss();
    double minor(int i, int j);
    double det();
    
};
#endif