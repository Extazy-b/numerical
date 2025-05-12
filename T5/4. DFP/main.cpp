#include "../libs/math.cpp"
#include <fstream>
#include <iostream>

using namespace std;

ofstream logFile("output.log");

int dimension = 2;
Poly targetFunction(dimension, 3, 0);
vector<double> start_point(dimension, 0);
vector<vector<double>> hess(dimension, vector <double> (dimension, 0));

const double EPSILON = pow(2, -60);
const double MAX_STEP = pow(2, -5);
const double MIN_STEP = pow(2, -25);

Poly getDirectionFunction(vector<Poly>& gradient){
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
    
    vector<Poly> DirectionLine = emptyPoint - stepPoly * gradient;
    Poly functionByDirection = targetFunction.compose(DirectionLine);

    return functionByDirection;
}

double goldenSectionSearch(double a, double b, double epsilon, const Poly& func, const vector<double>& point) {
    const double PHI = (sqrt(5) - 1) / 2;
    size_t dim = point.size();

    logFile << "\n\t\t\t\t == Запуск метода золотого сечения ==" << endl;
    logFile << "Начальный интервал: [" << a << ", " << b << "]" << endl;
    logFile << "Требуемая точность: " << epsilon << endl;

    vector<double> arg1(dim + 1, 0);
    vector<double> arg2(dim + 1, 0);
    copy(point.begin(), point.end(), arg1.begin());
    copy(point.begin(), point.end(), arg2.begin());

    arg1[dim] = b - PHI * (b - a);
    arg2[dim] = a + PHI * (b - a);
    double y1 = func.evaluate(arg1);
    double y2 = func.evaluate(arg2);
    
    int iteration = 1;
    while ((b - a) > epsilon) {
        logFile << "\nИтерация золотого сечения " << iteration << ":" << endl;
        logFile << "Текущий интервал: [" << a << ", " << b << "]" << endl;
        logFile << "Точки деления: x1 = " << arg1[dim] << ", x2 = " << arg2[dim] << endl;
        logFile << "Значения функции: f(x1) = " << y1 << ", f(x2) = " << y2 << endl;

        if (y1 < y2) {
            b = arg2[dim];
            arg2[dim] = arg1[dim];
            y2 = y1;

            arg1[dim] = b - PHI * (b - a);
            y1 = func.evaluate(arg1);
            logFile << "Сужение интервала справа" << endl;
        } else {
            a = arg1[dim];
            arg1[dim] = arg2[dim];
            y1 = y2;

            arg2[dim] = a + PHI * (b - a);
            y2 = func.evaluate(arg2);
            logFile << "Сужение интервала слева" << endl;
        }
        iteration++;
    }

    double result = (a + b) / 2;
    logFile << "\nРезультат метода золотого сечения:" << endl;
    logFile << "Найденная точка: " << result << endl;
    logFile << "Длина финального интервала: " << (b - a) << endl;
    return result;
}

double getOptimalStep(const vector<double>& point, const Poly& func, const Poly& der1, const Poly& der2) {
    
    logFile << "\n\t\t\t == Поиск оптимального шага ==" << endl;
    
    double startStep = goldenSectionSearch(0, MAX_STEP, pow(2, -7), func, point);
    logFile << "Начальный шаг (метод золотого сечения): " << startStep << endl;

    size_t dim = point.size();
    const double currentEpsilon = pow(2, -5);
    double delta = currentEpsilon + 1;
    vector<double> arg(dim + 1, 0);
    copy(point.begin(), point.end(), arg.begin());
    arg[dim] = startStep;

    int iteration = 1;
    while (abs(delta) > currentEpsilon) {
        logFile << "\nИтерация поиска шага " << iteration << ":" << endl;
        logFile << "Текущий шаг: " << arg[dim] << endl;
        
        if (arg[dim] > MAX_STEP) {
            logFile << "Достигнут максимальный размер шага" << endl;
    
            return MAX_STEP;
        }
        if (arg[dim] < MIN_STEP){
            logFile << "Достигнут минимальный размер шага" << endl;
    
            return MIN_STEP;
        }

        double der1_val = der1.evaluate(arg);
        double der2_val = der2.evaluate(arg);
        logFile << "Значение первой производной: " << der1_val << endl;
        logFile << "Значение второй производной: " << der2_val << endl;

        delta = der1_val / abs(der2_val);

        if (isinf(delta)){
            logFile << "Обнаружено деление на ноль, возвращаем минимальный шаг" << endl;
    
            return MIN_STEP;
        }
        if (isnan(delta)){
            logFile << "Некорректное значение, возвращаем максимальный шаг" << endl;
    
            return MAX_STEP;
        }

        logFile << "Текущая дельта: " << delta << endl;
        arg[dim] -= delta;
        logFile << "Новое значение шага: " << arg[dim] << endl;
        iteration++;
    }
     
    logFile << "\nРезультат поиска оптимального шага:" << endl;
    logFile << "Найденный шаг: " << arg[dim] << endl;
    logFile << "Финальная дельта: " << delta << endl;
    return arg[dim];
}

int main(){
    targetFunction.setCoefs("6 3 0 1 1 1 -1 0 2 1 1 0 -2 0 1 3 0 0 -4");
    vector<Poly> gradient = Nabla(targetFunction);
    Poly functionByDirection = getDirectionFunction(gradient);
    Poly der1 = derivate(functionByDirection, dimension);
    Poly der2 = derivate(der1, dimension);

    for (size_t i = 0; i < dimension; i++) hess[i][i] = 1;
    
    size_t k = 1;
    vector <double> gradienValue(dimension, 0);
    vector <double> last_point(dimension, 0);
    vector <double> next_point(dimension, 0);
    vector <double> middle_point(dimension, 0);

    vector <double> direction(dimension, 0);
    vector <double> P(dimension, 0);
    vector <double> q(dimension, 0);
    vector <vector <double>> hessDelta(dimension, vector<double>(dimension, 0));
    gradienValue = evaluate2list(gradient, middle_point);

    copy(start_point.begin(), start_point.end(), last_point.begin());
    copy(start_point.begin(), start_point.end(), next_point.begin());
    copy(next_point.begin(), next_point.end(), middle_point.begin());


    cout << "Функция: " << targetFunction << endl;
    cout << "Градиент: " << gradient << endl;
    cout << "Погрешность: " << EPSILON << endl;
    cout << "Границы шага: [" << MIN_STEP << " - " << MAX_STEP << ']' << endl;
    cout << "Начальная точка: " << start_point << endl;
    cout << "Градиент в начальной точке: " << gradienValue << endl;

    logFile << "Функция: " << targetFunction << endl;
    logFile << "Градиент: " << gradient << endl;
    logFile << "Погрешность: " << EPSILON << endl;
    logFile << "Границы шага: [" << MIN_STEP << " - " << MAX_STEP << ']' << endl;
    logFile << "Начальная точка: " << start_point << endl;
    logFile << "Градиент в начальной точке: " << gradienValue << endl;

    cin.get();

    cout << "\t\t == Начало оптимизации ==" << endl;
    logFile << "\t\t == Начало оптимизации ==" << endl;

    while (calculateNorm(gradienValue) >= EPSILON){  
        cout << "\nИтерация " << k << endl;
        logFile << "\nИтерация " << k << endl;  
        
        copy(next_point.begin(), next_point.end(), last_point.begin());
                
        for (size_t j = 0; j<dimension; ++j){
            cout << "Подитерация " << j + 1 << ":" << endl;
            logFile << "Подитерация " << j + 1 << ":" << endl;

            gradienValue = evaluate2list(gradient, middle_point);
            direction = (-1) * (hess * gradienValue);           

    
            
            double step = getOptimalStep(middle_point, functionByDirection, der1, der2);
            
            P = step * direction;

            logFile << "Оптимальный шаг: " << step << endl;

            logFile << "Значение градиента: " << gradienValue << endl;

            logFile << "Направление: " << direction << endl;

            logFile << "Вектор P: " << P << endl;

            logFile << "Матрица H: " << endl << hess << endl;

            middle_point = middle_point + P;
            
            cout << "Текущая точка: " << middle_point << endl;
            logFile << "Текущая точка: " << middle_point << endl;
            
            q = evaluate2list(gradient, middle_point) - gradienValue;
            hessDelta = ((1/(P*q)) * vector2matrix(P)) - ((1/(q * (hess * q))) * (hess * vector2matrix(q) * hess));
            hess = hess + hessDelta;
            
            logFile << "Вектор q: " << q << endl;

            logFile << "Изменение матрицы H: " << endl << hessDelta << endl;

    
        }

          
        copy(middle_point.begin(), middle_point.end(), next_point.begin());
        
        double norm = calculateNorm(next_point - last_point);
        double value = targetFunction.evaluate(next_point);

        cout << "Норма разности точек: " << norm << endl;
        logFile << "Норма разности точек: " << norm << endl;
        cout << "Значение функции: " << value << endl;
        logFile << "Значение функции: " << value << endl;
        cout << "Изменение значения функции: " << abs(value - targetFunction.evaluate(last_point)) << endl;
        logFile << "Изменение значения функции: " << abs(value - targetFunction.evaluate(last_point)) << endl;
        
        if ((norm < EPSILON) || (abs(value - targetFunction.evaluate(last_point)) < EPSILON)){
            cout << "Достигнута требуемая точность. Завершение." << endl;
            logFile << "Достигнута требуемая точность. Завершение." << endl;
            break;

        }

        k++;
    }
    
    cout << "\nРезультат оптимизации:" << endl;
    logFile << "\nРезультат оптимизации:" << endl;
    cout << "Найденная точка: ";
    logFile << "Найденная точка: ";
    for(const auto& val : next_point) {
        cout << val << " ";
        logFile << val << " ";
    }
    cout << endl;
    logFile << endl;
    logFile.close();

    return 0;
}