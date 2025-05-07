// f(x, y) = x^3 - x*y + y^2 - 2x + 3y - 4
// f(x, y) -> min
// Начальная точка X0 = (0, 0)

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include "../math.cpp"

using namespace std;

// Константы алгоритма
const double STEP = pow(2.0, -10);  // Шаг градиентного спуска
const double EPSILON = pow(2.0, -60);  // Погрешность для условия остановки

// Строка с коэффициентами полинома
// f(x, y) = x^3 - x*y + y^2 - 2x + 3y - 4
const string FUNCTION_COEFFICIENTS = "6 3 0 1 1 1 -1 0 2 1 1 0 -2 0 1 3 0 0 -4";


int main(){
    // Инициализация целевой функции
    Poly targetFunction(2, 3);  // 2 переменные, 3-я степень
    targetFunction.setCoefs(FUNCTION_COEFFICIENTS);

    // Инициализация точек
    const size_t dimension = targetFunction.getVarc();
    vector<double> lastPoint(dimension, 0.0);   // Предыдущая точка
    vector<double> nextPoint(dimension, 0.0);   // Текущая точка

    // Вычисление градиента функции
    vector<Poly> gradient = Nabla(targetFunction);

    // Настройка вывода в научной нотации
    cout.setf(ios::scientific);

    // Вывод параметров алгоритма
    cout << "Параметры алгоритма:\n";
    cout << "Начальная точка: " << lastPoint << endl;
    cout << "Погрешность: " << EPSILON << endl;

    cout << "Нажмите Enter для начала оптимизации..." << endl;
    cin.get();

    size_t j = 0;
    vector<double> direction = (-1) * evaluate2list(gradient, lastPoint);

    while (true)
    {
        
        /* code */
        if ((calculateNorm(direction) < EPSILON) || calculateNorm(lastPoint - nextPoint) < EPSILON)
        {
            cout << nextPoint << endl;
            cout << j << endl;
            cout << targetFunction.evaluate(nextPoint) << endl;
            break;
        }
        
    }
    


    return 0;
}
