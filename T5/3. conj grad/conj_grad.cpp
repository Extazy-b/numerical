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
const double MAX_STEP = pow(2, 0);  // Максимальный размер шага
const double EPSILON = pow(2.0, -60);  // Погрешность для условия остановки
Poly targetFunction(2, 3);  // 2 переменные, 3-я степень

// Строка с коэффициентами полинома
// f(x, y) = x^3 - x*y + y^2 - 2x + 3y - 4
const string FUNCTION_COEFFICIENTS = "6 3 0 1 1 1 -1 0 2 1 1 0 -2 0 1 3 0 0 -4";



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



int main(){
    // Инициализация целевой функции
    targetFunction.setCoefs(FUNCTION_COEFFICIENTS);

    // Инициализация точек
    const size_t dimension = targetFunction.getVarc();
    vector<double> lastPoint(dimension, 0.0);   // Предыдущая точка
    vector<double> tmpPoint(dimension, 0.0);   // Промежуточкая точка
    vector<double> nextPoint(dimension, 0.0);   // Текущая точка

    // Вычисление градиента функции
    vector<Poly> gradient = Nabla(targetFunction);

    // Функция вдоль направления
    Poly stepPoly(dimension + 1, 1, 0);
    vector<int> indices(dimension + 1, 0);
    indices[dimension] = 1;
    stepPoly(indices) = 1;
    
    vector<Poly> emptyPoint(dimension, Poly(dimension, 1, 0));
    for (size_t i = 0; i < dimension; ++i) {
        indices = vector<int>(dimension, 0);
        indices[i] = 1;
        emptyPoint[i](indices) = 1;
    }
    
    vector<Poly> DirectionLine = emptyPoint - stepPoly * gradient; // Pk - step*grad f(Pk)
    Poly functionByDirection = targetFunction.compose(DirectionLine); // f(Pk - step*grad f(Pk))
    
    // Производные вспомогательной функции
    Poly der1 = derivate(functionByDirection, dimension);
    Poly der2 = derivate(der1, dimension);

    // Настройка вывода в научной нотации
    cout.setf(ios::scientific);

    // Вывод параметров алгоритма
    cout << "Параметры алгоритма:\n";
    cout << "Начальная точка: " << lastPoint << endl;
    cout << "Погрешность: " << EPSILON << endl;

    cout << "Нажмите Enter для начала оптимизации..." << endl;
    cin.get();

    size_t k = 0;
    vector<double> direction = (-1) * evaluate2list(gradient, lastPoint);
    double step;
    double omega;
    bool isDone = false;

    while (true)
    {
        lastPoint = nextPoint;
        vector<double> direction = (-1) * evaluate2list(gradient, lastPoint);
        for (size_t j = 0; j < dimension; ++j){
            step = getOptimalStep(lastPoint, functionByDirection, der1, der2);
            
            tmpPoint = lastPoint + (step * direction);

            omega = calculateNorm(evaluate2list(gradient, tmpPoint)) / calculateNorm(evaluate2list(gradient, lastPoint));
            omega = omega * omega;

            direction = (-1)*evaluate2list(gradient, tmpPoint) + omega * direction;

            cout << "j: " << j << endl;
            cout << "x_k^j: " << lastPoint << endl;
            cout << "x_l^j+1: " << tmpPoint << endl;
            cout << "S_k^j+1: " << direction << endl;
            cout << "step: " << step << endl;
            cout << "value: " << targetFunction.evaluate(tmpPoint) << endl; 
            cout << endl;

            if ((calculateNorm(direction) < EPSILON) || calculateNorm(lastPoint - tmpPoint) < EPSILON)
            {
                cout << "LastPoint: " << lastPoint << endl;
                cout << "NextPoint: " << nextPoint << endl;
                cout << "j: " << j << endl;
                cout << "value: " << targetFunction.evaluate(nextPoint) << endl;
                isDone = true;
                break;
            }
        }
        if (isDone){
            cout << k << endl;
            break;
        }
        nextPoint = tmpPoint;
        cout << "k: " << k << endl;
        ++k;
    }
    return 0;
}
