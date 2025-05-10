// f(x, y) = x^3 - x*y + y^2 - 2x + 3y - 4
// f(x, y) -> min
// Начальная точка X0 = (0, 0)

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include "../libs/math.cpp"

using namespace std;

const double MAX_STEP = pow(2, 0);  // Максимальный размер шага
const float EPSILON = pow(2, -60);  // Точность вычислений
Poly poly(2, 3);

// Проверка критериев остановки алгоритма
bool checkStoppingCriteria(const vector<double>& lastPoint,
                          const vector<double>& nextPoint,
                          const Poly& func,
                          const vector<Poly>& grad) {
    if (calculateNorm(lastPoint - nextPoint) < EPSILON) {
        return true;
    }
    if (calculateNorm(evaluate2list(grad, nextPoint)) < EPSILON) {
        return true;
    }
    if (abs(func.evaluate(lastPoint) - func.evaluate(nextPoint)) < EPSILON) {
        return true;
    }
    return false;
}

// Метод золотого сечения для поиска минимума функции одной переменной
double goldenSectionSearch(double a, double b, double epsilon, const Poly& func, const vector<double>& point) {
    const double PHI = (sqrt(5) - 1) / 2;  // Золотое сечение
    size_t dim = point.size();

    vector<double> arg1(dim + 1, 0);
    vector<double> arg2(dim + 1, 0);
    copy(point.begin(), point.end(), arg1.begin());
    copy(point.begin(), point.end(), arg2.begin());

    arg1[dim] = b - PHI * (b - a);
    arg2[dim] = a + PHI * (b - a);
    double y1 = func.evaluate(arg1);
    double y2 = func.evaluate(arg2);
    
    while ((b - a) > epsilon) {
        if (y1 < y2) {
            b = arg2[dim];
            arg2[dim] = arg1[dim];
            y2 = y1;

            arg1[dim] = b - PHI * (b - a);
            y1 = func.evaluate(arg1);
        } else {
            a = arg1[dim];
            arg1[dim] = arg2[dim];
            y1 = y2;

            arg2[dim] = a + PHI * (b - a);
            y2 = func.evaluate(arg2);
        }
    }

    return (a + b) / 2;
}

// Получение оптимального размера шага методом Ньютона
double getOptimalStep(const vector<double>& point, const Poly& func, const Poly& der1, const Poly& der2) {
    double startStep = goldenSectionSearch(0, MAX_STEP, pow(2, 0), func, point);

    size_t dim = point.size();
    const double currentEpsilon = pow(2, -5);
    double delta = currentEpsilon + 1;
    vector<double> arg(dim + 1, 0);
    copy(point.begin(), point.end(), arg.begin());
    arg[dim] = startStep;

    while (abs(delta) > currentEpsilon) {
        if (arg[dim] > MAX_STEP) {
            arg[dim] = MAX_STEP;
            break;
        }
        delta = der1.evaluate(arg) / der2.evaluate(arg);
        arg[dim] -= delta;
    }
     
    return arg[dim];
}

int main(int argc, char* argv[]) {
    // Инициализация
    poly.setCoefs("6 3 0 1 1 1 -1 0 2 1 1 0 -2 0 1 3 0 0 -4"); // f(x, y) = x^3 - x*y + y^2 - 2x + 3y - 4
    
    const size_t variableCount = poly.getVarc();
    vector<double> lastPoint(variableCount, 0);
    vector<double> nextPoint(variableCount, 0);
    
    // Вычисление градиента
    vector<Poly> grad = Nabla(poly);
    cout << "Функция: " << poly << endl << endl;
    
    // Вспомогательная функция
    Poly stepPoly(variableCount + 1, 1, 0);
    vector<int> indices(variableCount + 1, 0);
    indices[variableCount] = 1;
    stepPoly(indices) = 1;
    
    vector<Poly> emptyPoint(variableCount, Poly(variableCount, 1, 0));
    for (size_t i = 0; i < variableCount; ++i) {
        indices = vector<int>(variableCount, 0);
        indices[i] = 1;
        emptyPoint[i](indices) = 1;
    }
    
    vector<Poly> func = emptyPoint - stepPoly * grad; // Pk - step*grad f(Pk)
    Poly extraPoly = poly.compose(func); // f(Pk - step*grad f(Pk))
    
    // Производные вспомогательной функции
    Poly der1 = derivate(extraPoly, variableCount);
    Poly der2 = derivate(der1, variableCount);

    cout.setf(ios::scientific);
    ofstream outputFile("output.txt");
    cout << "Допустимая погрешность: " << EPSILON << endl;
    cout << "Нажмите Enter для начала вычислений" << endl;
    
    cin.get();

    size_t iteration = 1;
    double step = 0;
    outputFile << "Итерация | Точка (x, y) | Дельта | Значение функции | Размер шага" << endl;
    
    while (true) {
        step = getOptimalStep(lastPoint, extraPoly, der1, der2);
        nextPoint = lastPoint - (step * evaluate2list(grad, lastPoint));

        outputFile << iteration << " | (" << lastPoint[0] << ", " << lastPoint[1] << ") | ";
        outputFile << calculateNorm(nextPoint - lastPoint) << " | " << poly.evaluate(nextPoint) << " | " << step << endl;
        
        cout << "Итерация " << iteration << ": (" << lastPoint[0] << ", " << lastPoint[1] << ") | ";
        cout << "Дельта: " << calculateNorm(nextPoint - lastPoint) << " | Значение: " << poly.evaluate(nextPoint) << " | Шаг: " << step << endl;

        
        double currentValue = poly.evaluate(nextPoint);
        if (std::isinf(currentValue) || std::isnan(currentValue)) {
            cout << "\nПредупреждение: Достигнуто некорректное значение функции!" << endl;
            cout << "Последнее корректное значение: " << poly.evaluate(lastPoint) << endl;
            cout << "Последняя корректная точка: (" << lastPoint[0] << ", " << lastPoint[1] << ")" << endl;
            return 1;
        }
        

        if (checkStoppingCriteria(lastPoint, nextPoint, poly, grad)) {
            cout << "Вычисление остановлено на итерации " << iteration << endl;
            break;
        }
        
        ++iteration;
        lastPoint = nextPoint;
    }
    
    cout << "Результат: X = (" << nextPoint[0] << ", " << nextPoint[1] << ")" << endl;
    cout << "Значение функции: " << poly.evaluate(nextPoint) << endl;
    
    return 0;
}