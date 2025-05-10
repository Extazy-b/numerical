// f(x, y) = x^3 - x*y + y^2 - 2x + 3y - 4
// grad f = [3*x^2 - y - 2; -x + 2y + 3]
// f(x, y) -> min
// Начальная точка X0 = (0, 0)

#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>
#include "../libs/math.cpp"

using namespace std;

// Константы алгоритма
const double STEP = pow(2.0, -10);  // Шаг градиентного спуска
const double EPSILON = pow(2.0, -60);  // Погрешность для условия остановки

// Строка с коэффициентами полинома
// f(x, y) = x^3 - x*y + y^2 - 2x + 3y - 4
const string FUNCTION_COEFFICIENTS = "6 3 0 1 1 1 -1 0 2 1 1 0 -2 0 1 3 0 0 -4";

/**
 * Проверка условий остановки алгоритма
 * @param lastPoint Предыдущая точка
 * @param nextPoint Текущая точка
 * @param function Целевая функция
 * @param gradient Градиент функции
 * @return true если выполнено хотя бы одно условие остановки
 */
bool checkStopConditions(const vector<double>& lastPoint,
                        const vector<double>& nextPoint,
                        const Poly& function,
                        const vector<Poly>& gradient,
                        ofstream& logFile) {
    // Проверяем три условия остановки
    if (calculateNorm(lastPoint - nextPoint) < EPSILON) {
        logFile << "Достигнута требуемая точность по аргументу\n";
        logFile << "Норма разности точек: " << calculateNorm(lastPoint - nextPoint) << "\n";
        std::cout << "Достигнута требуемая точность по аргументу\n";
        return true;
    }
    if (calculateNorm(evaluate2list(gradient, nextPoint)) < EPSILON) {
        logFile << "Достигнута требуемая точность по градиенту\n";
        logFile << "Норма градиента: " << calculateNorm(evaluate2list(gradient, nextPoint)) << "\n";
        std::cout << "Достигнута требуемая точность по градиенту\n";
        return true;
    }
    if (std::abs(function.evaluate(lastPoint) - function.evaluate(nextPoint)) < EPSILON) {
        logFile << "Достигнута требуемая точность по значению функции\n";
        logFile << "Разность значений функции: " << std::abs(function.evaluate(lastPoint) - function.evaluate(nextPoint)) << "\n";
        std::cout << "Достигнута требуемая точность по значению функции\n";
        return true;
    }
    return false;
}

int main(int argc, char* argv[]) {
    // Открываем файл для логирования
    ofstream logFile("output.log");
    logFile << fixed << setprecision(10);

    // Инициализация целевой функции
    Poly targetFunction(2, 3);  // 2 переменные, 3-я степень
    targetFunction.setCoefs(FUNCTION_COEFFICIENTS);
    
    logFile << "Начало работы алгоритма\n";
    logFile << "Шаг градиентного спуска: " << STEP << "\n";
    logFile << "Погрешность для условия остановки: " << EPSILON << "\n\n";

    // Инициализация точек
    const size_t dimension = targetFunction.getVarc();
    vector<double> lastPoint(dimension, 0.0);   // Предыдущая точка
    vector<double> nextPoint(dimension, 0.0);   // Текущая точка

    // Вычисление градиента функции
    vector<Poly> gradient = Nabla(targetFunction);

    cout << "Параметры работы алгоритма" << endl;
    cout << "Функция: " << targetFunction << endl;
    cout << "Градиент: " << gradient << endl;
    cout << "Шаг градиентного спуска: " << STEP << endl;
    cout << "Погрешность для условия остановки: " << EPSILON << endl;
    cin.get();
    cout << endl;

    // Основной цикл градиентного спуска
    size_t iterationCount = 1;
    
    while (true) {
        // Вычисление следующей точки x_{k+1} = x_k - h*grad(f)(x_k)
        nextPoint = lastPoint - STEP * evaluate2list(gradient, lastPoint);
        
        // Логирование текущего состояния
        logFile << "Итерация " << iterationCount << ":\n";
        logFile << "Текущая точка: (";
        for (size_t i = 0; i < dimension; ++i) {
            logFile << nextPoint[i] << (i < dimension-1 ? ", " : ")\n");
        }
        logFile << "Значение функции: " << targetFunction.evaluate(nextPoint) << "\n";
        logFile << "Норма градиента: " << calculateNorm(evaluate2list(gradient, nextPoint)) << "\n\n";

        cout << "Итерация " << iterationCount << ", f(x) = " << targetFunction.evaluate(nextPoint) << endl;
        
        // Проверка условий остановки
        if (checkStopConditions(lastPoint, nextPoint, targetFunction, gradient, logFile)) {
            break;
        }
        
        iterationCount++;
        lastPoint = nextPoint;
    }

    // Итоговый вывод
    logFile << "\nРезультаты работы алгоритма:\n";
    logFile << "Количество итераций: " << iterationCount << "\n";
    logFile << "Найденная точка минимума: (";
    for (size_t i = 0; i < dimension; ++i) {
        logFile << nextPoint[i] << (i < dimension-1 ? ", " : ")\n");
    }
    logFile << "Значение функции в точке минимума: " << targetFunction.evaluate(nextPoint) << "\n";
    
    cout << "\nОптимизация завершена. Подробности в файле gradient_descent.log\n";
    
    logFile.close();
    return 0;
}