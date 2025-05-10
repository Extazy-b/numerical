// f(x, y) = x^3 - x*y + y^2 - 2x + 3y - 4
// f(x, y) -> min
// Начальная точка X0 = (0, 0)

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>
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
    double normDiff = calculateNorm(lastPoint - nextPoint);
    double normGrad = calculateNorm(evaluate2list(grad, nextPoint));
    double funcDiff = abs(func.evaluate(lastPoint) - func.evaluate(nextPoint));
    
    ofstream logFile("output.log", ios::app);
    logFile << "Проверка критериев остановки:" << endl;
    logFile << "Норма разности точек: " << scientific << normDiff << endl;
    logFile << "Норма градиента: " << scientific << normGrad << endl;
    logFile << "Разность значений функции: " << scientific << funcDiff << endl;
    logFile.close();

    if (normDiff < EPSILON) {
        return true;
    }
    if (normGrad < EPSILON) {
        return true;
    }
    if (funcDiff < EPSILON) {
        return true;
    }
    return false;
}

// Метод золотого сечения для поиска минимума функции одной переменной
double goldenSectionSearch(double a, double b, double epsilon, const Poly& func, const vector<double>& point) {
    const double PHI = (sqrt(5) - 1) / 2;  // Золотое сечение
    size_t dim = point.size();

    ofstream logFile("output.log", ios::app);
    logFile << "\tМетод золотого сечения:" << endl;
    logFile << "\tНачальный интервал: [" << a << ", " << b << "]" << endl;

    vector<double> arg1(dim + 1, 0);
    vector<double> arg2(dim + 1, 0);
    copy(point.begin(), point.end(), arg1.begin());
    copy(point.begin(), point.end(), arg2.begin());

    arg1[dim] = b - PHI * (b - a);
    arg2[dim] = a + PHI * (b - a);
    double y1 = func.evaluate(arg1);
    double y2 = func.evaluate(arg2);
    
    while ((b - a) > epsilon) {
        logFile << "\tИнтервал: [" << a << ", " << b << "], длина: " << b-a << endl;
        
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

    double result = (a + b) / 2;
    logFile << "\tНайденный оптимальный шаг: " << result << endl;
    logFile.close();
    return result;
}

// Получение оптимального размера шага методом Ньютона
double getOptimalStep(const vector<double>& point, const Poly& func, const Poly& der1, const Poly& der2) {
    ofstream logFile("output.log", ios::app);
    logFile << "\nПоиск оптимального шага методом Ньютона:" << endl;
    
    double startStep = goldenSectionSearch(0, MAX_STEP, pow(2, 0), func, point);
    logFile << "Начальное приближение шага: " << startStep << endl;

    size_t dim = point.size();
    const double currentEpsilon = pow(2, -5);
    double delta = currentEpsilon + 1;
    vector<double> arg(dim + 1, 0);
    copy(point.begin(), point.end(), arg.begin());
    arg[dim] = startStep;

    while (abs(delta) > currentEpsilon) {
        if (arg[dim] > MAX_STEP) {
            arg[dim] = MAX_STEP;
            logFile << "Достигнут максимальный размер шага" << endl;
            break;
        }
        delta = der1.evaluate(arg) / der2.evaluate(arg);
        arg[dim] -= delta;
        logFile << "Текущий шаг: " << arg[dim] << ", delta: " << delta << endl;
    }
     
    logFile << "Итоговый оптимальный шаг: " << arg[dim] << endl;
    logFile.close();
    return arg[dim];
}

int main(int argc, char* argv[]) {
    ofstream logFile("output.log");

    // Инициализация
    poly.setCoefs("6 3 0 1 1 1 -1 0 2 1 1 0 -2 0 1 3 0 0 -4"); // f(x, y) = x^3 - x*y + y^2 - 2x + 3y - 4
    
    const size_t variableCount = poly.getVarc();
    vector<double> lastPoint(variableCount, 0);
    vector<double> nextPoint(variableCount, 0);
    
    // Вычисление градиента
    vector<Poly> grad = Nabla(poly);
    
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


    logFile << "Начало оптимизации" << endl;
    logFile << "Целевая функция: " << poly << endl;
    logFile << "Градиент: " << grad << endl;
    logFile << "Целевая функция: " << poly << endl;
    logFile << "Максимальный шаг: " << MAX_STEP << endl;
    logFile << "Погрешность: " << EPSILON << endl;
    logFile.close();

    cout << "Начало оптимизации" << endl;
    cout << "Градиент: " << grad << endl;
    cout << "Целевая функция: " << poly << endl;
    cout << "Максимальный шаг: " << MAX_STEP << endl;
    cout << "Погрешность: " << EPSILON << endl;
    cout.setf(ios::scientific);

    cin.get();

    size_t iteration = 1;
    double step = 0;
    
    while (true) {
        cout << "\nИтерация " << iteration << ":" << endl;
        cout << "Текущая точка: (" << lastPoint[0] << ", " << lastPoint[1] << ")" << endl;
        
        ofstream logFile("output.log", ios::app);
        logFile << "\n=== Итерация " << iteration << " ===" << endl;
        logFile << "Текущая точка: (" << lastPoint[0] << ", " << lastPoint[1] << ")" << endl;
        logFile << "Значение функции: " << poly.evaluate(lastPoint) << endl;
        logFile << "Градиент в точке: (" << evaluate2list(grad, lastPoint)[0] << ", " 
                << evaluate2list(grad, lastPoint)[1] << ")" << endl;
        logFile.close();

        step = getOptimalStep(lastPoint, extraPoly, der1, der2);
        nextPoint = lastPoint - (step * evaluate2list(grad, lastPoint));

        double currentValue = poly.evaluate(nextPoint);
        if (std::isinf(currentValue) || std::isnan(currentValue)) {
            cout << "Ошибка: значение функции стало бесконечным или NaN" << endl;
            ofstream logFile("output.log", ios::app);
            logFile << "Ошибка: значение функции стало бесконечным или NaN" << endl;
            logFile.close();
            return 1;
        }
        
        cout << "Найденный шаг: " << step << endl;
        cout << "Новая точка: (" << nextPoint[0] << ", " << nextPoint[1] << ")" << endl;
        cout << "Значение функции: " << currentValue << endl;

        if (checkStoppingCriteria(lastPoint, nextPoint, poly, grad)) {
            cout << "\nОптимизация завершена!" << endl;
            cout << "Найден минимум в точке: (" << nextPoint[0] << ", " << nextPoint[1] << ")" << endl;
            cout << "Значение функции в минимуме: " << poly.evaluate(nextPoint) << endl;
            
            ofstream logFile("output.log", ios::app);
            logFile << "\nОптимизация завершена!" << endl;
            logFile << "Найден минимум в точке: (" << nextPoint[0] << ", " << nextPoint[1] << ")" << endl;
            logFile << "Значение функции в минимуме: " << poly.evaluate(nextPoint) << endl;
            logFile << "Общее количество итераций: " << iteration << endl;
            logFile.close();
            break;
        }
        
        ++iteration;
        lastPoint = nextPoint;
    }
    return 0;
}