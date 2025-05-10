// f(x, y) = x^3 - x*y + y^2 - 2x + 3y - 4
// grad f = [3*x^2 - y - 2; -x + 2y + 3]
// f(x, y) -> min
// Начальная точка X0 = (0, 0)

#include <vector>
#include <iostream>
#include <cmath>
#include <string>
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
                        const vector<Poly>& gradient) {
    // Проверяем три условия остановки
    if (calculateNorm(lastPoint - nextPoint) < EPSILON) {
        std::cout << "Достигнута требуемая точность по аргументу\n";
        return true;
    }
    if (calculateNorm(evaluate2list(gradient, nextPoint)) < EPSILON) {
        std::cout << "Достигнута требуемая точность по градиенту\n";
        return true;
    }
    if (std::abs(function.evaluate(lastPoint) - function.evaluate(nextPoint)) < EPSILON) {
        std::cout << "Достигнута требуемая точность по значению функции\n";
        return true;
    }
    return false;
}

int main(int argc, char* argv[]) {
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
    cout << "Шаг: " << STEP << endl;
    cout << "Погрешность: " << EPSILON << endl;

    cout << "Нажмите Enter для начала оптимизации..." << endl;
    cin.get();

    // Основной цикл градиентного спуска
    size_t iterationCount = 1;
    cout << "\nНачало итерационного процесса:\n";
    
    while (true) {
        // Вычисление следующей точки x_{k+1} = x_k - h*grad(f)(x_k)
        nextPoint = lastPoint - STEP * evaluate2list(gradient, lastPoint);

        // Вывод информации о текущей итерации
        cout << "Итерация " << iterationCount << ":\n";
        cout << "Текущая точка: " << nextPoint;
        cout << "Изменение: " << calculateNorm(nextPoint - lastPoint) 
             << " | Значение функции: " << targetFunction.evaluate(nextPoint) << endl;
        
        // Проверка условий остановки
        if (checkStopConditions(lastPoint, nextPoint, targetFunction, gradient)) {
            cout << "\nОптимизация завершена на итерации: " << iterationCount << endl;
            break;
        }
        
        iterationCount++;
        lastPoint = nextPoint;
    }

    // Вывод результатов оптимизации
    cout << "\nРезультаты оптимизации:\n";
    cout << "Найденная точка минимума: [" << nextPoint[0] << ", " << nextPoint[1] << "]\n";
    cout << "Значение функции в точке минимума: " << targetFunction.evaluate(nextPoint) << endl;
    
    return 0;
}