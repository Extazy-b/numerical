/*
Алгоритм оптимизации методом сопряженных градиентов:

1. Вычисляется градиент целевой функции в текущей точке
2. Определяется направление спуска на основе градиента:
   - На первой итерации направление спуска равно антиградиенту
   - На последующих итерациях направление спуска вычисляется как линейная комбинация текущего антиградиента и предыдущего направления спуска
   - Коэффициент для линейной комбинации рассчитывается по формуле Флетчера-Ривса или Полака-Рибьера
3. Находится оптимальный шаг вдоль направления спуска
4. Выполняется перемещение в новую точку
5. Процесс повторяется до достижения условия сходимости
*/

#include "../libs/math.cpp"
#include <fstream>
#include <iostream>

using namespace std;

ofstream logFile("output.log");

int dimension = 2;
Poly targetFunction(dimension, 3, 0);
vector<double> start_point(dimension, 0);

const double EPSILON = pow(2, -60);
const double MIN_STEP = pow(2, -60);
const double MAX_STEP = pow(2, 0);

vector<double> old_point(dimension, 0);
vector<double> new_point(dimension, 0);
vector<Poly> gradient(dimension, Poly());
vector<double> gradientValue(dimension, 0);
double step;

size_t k;

Poly getDirectionFunction(vector<double> point){
    Poly function(1, 0, 0);
    vector<Poly> extraFunction(dimension, Poly(1, 0, 0));

    for (size_t i = 0; i < dimension; i++)
    {
        extraFunction[i] = Poly(0, 0, point[i]) + gradientValue[i] * (Poly(1, 1, 1) - Poly(1, 0, 1));
    }

    for (size_t ind = 0; ind < targetFunction.getLen(); ind++)
    {
        Poly summand(0, 0, targetFunction[ind]);
        vector<int> indices = targetFunction.ind2multyind(ind);
        for (size_t subind = 0; subind < dimension; subind++){
            summand = summand * pow(extraFunction[subind], indices[subind]);
        }
        function  = function + summand;
    }

    return function;
}
double getOptimalStep(vector <double> point){
    double step;
    Poly DirectionFunction = getDirectionFunction(point);
    double delta;

    double startStep = goldenSectionSearch(MIN_STEP, MAX_STEP, pow(2, -10), DirectionFunction);
    step = Newton(DirectionFunction, startStep, MIN_STEP, MAX_STEP, EPSILON);

    return step;
}
void logging(size_t counter){
    cout << endl << "\tsubIter: " << counter;
    cout << "\tdirection: " << gradientValue << endl;
    cout << "\tstep: " << step << endl;
    cout << "\tnew_point: " << new_point << endl;
    cout << "\tpoints delta: " << calculateNorm(old_point - new_point) << endl;
    cout << "\tvalue: " << targetFunction.evaluate(new_point) << endl;
    cout << "\tvalue delta: " << targetFunction.evaluate(new_point) - targetFunction.evaluate(old_point) << endl;

    logFile << endl << "\tsubIter: " << counter;
    logFile << "\tdirection: " << gradientValue << endl;
    logFile << "\tstep: " << step << endl;
    logFile << "\tnew_point: " << new_point << endl;
    logFile << "\tpoints delta: " << calculateNorm(old_point - new_point) << endl;
    logFile << "\tvalue: " << targetFunction.evaluate(new_point) << endl;
    logFile << "\tvalue delta: " << targetFunction.evaluate(new_point) - targetFunction.evaluate(old_point) << endl;
}

int main(){
    targetFunction.setCoefs("6 3 0 1 1 1 -1 0 2 1 1 0 -2 0 1 3 0 0 -4");
    gradient = Nabla(targetFunction);

    copy(start_point.begin(), start_point.end(), old_point.begin());
    copy(start_point.begin(), start_point.end(), new_point.begin());

    double omega;
    bool isDone = false;

    cout << "Function:" << endl << targetFunction << endl;
    cout << "EPSILON: " << EPSILON << endl;
    cout << "MIN STEP: " << MIN_STEP << endl;
    cout << "MAX STEP: " << MAX_STEP << endl;
    cout << "press enter to start calculating" << endl;

    logFile << "Function:" << endl << targetFunction << endl;
    logFile << "EPSILON: " << EPSILON << endl;
    logFile << "MIN STEP: " << MIN_STEP << endl;
    logFile << "MAX STEP: " << MAX_STEP << endl;

    while (true)
    {
        gradientValue = (-1) * evaluate2list(gradient, old_point);

        cout << endl << "== Iter: " << k << " ==" << endl;
        logFile << endl << "== Iter: " << k << " ==" << endl;

        for (size_t j = 0; j < dimension; ++j){            
            step = getOptimalStep(old_point);
            new_point = old_point + (step * gradientValue);
            
            omega = calculateNorm(evaluate2list(gradient, new_point)) / calculateNorm(evaluate2list(gradient, old_point));
            omega = omega * omega;

            gradientValue = (-1)*evaluate2list(gradient, new_point) + omega * gradientValue;

            logging(j);

            if (calculateNorm(gradientValue) < EPSILON) {isDone = true; break;}
            if (calculateNorm(old_point - new_point) < EPSILON) {isDone = true; break;}
            if (abs(targetFunction.evaluate(new_point) - targetFunction.evaluate(old_point)) < EPSILON) {isDone = true; break;}


            copy(new_point.begin(), new_point.end(), old_point.begin());
        }
        cout << isDone;
        if (isDone) break;
        
        ++k;
    }

    cout << endl << "END" << endl;
    logFile.close(); 
    return 0;
}