/*
Алгоритм градиентного спуска для поиска минимума полиномиальной функции:
1. Задается начальная точка и целевая функция
2. На каждой итерации:
   - Вычисляется градиент функции в текущей точке
   - Делается шаг в направлении, противоположном градиенту
   - Проверяются условия остановки (малое изменение точки или градиента)
3. Процесс продолжается до достижения локального минимума
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
const double STEP = pow(2, -10);

vector<double> old_point(dimension, 0);
vector<double> new_point(dimension, 0);
vector<Poly> gradient(dimension, Poly());
vector<double> gradientValue(dimension, 0);

size_t i;

void logging(){
    cout << endl << "== Iter: " << i << " ==" << endl;
    cout << "new point: " << new_point << endl;
    cout << "points delta: " << calculateNorm(new_point - old_point) << endl;
    cout << "value: " << targetFunction.evaluate(new_point) << endl;
    cout << "gradientValue: " << gradientValue << endl;
    cout << "gradientNotm: " << calculateNorm(gradientValue) << endl;

    logFile << endl << "== Iter: " << i << " ==" << endl;
    logFile << "new point: " << new_point << endl;
    logFile << "points delta: " << calculateNorm(new_point - old_point) << endl;
    logFile << "value: " << targetFunction.evaluate(new_point) << endl;
    logFile << "gradientValue: " << gradientValue << endl;
    logFile << "gradientNotm: " << calculateNorm(gradientValue) << endl;
}

int main(){
    targetFunction.setCoefs("6 3 0 1 1 1 -1 0 2 1 1 0 -2 0 1 3 0 0 -4");
    vector<Poly> gradient = Nabla(targetFunction);
    copy(start_point.begin(), start_point.end(), old_point.begin());
    copy(start_point.begin(), start_point.end(), new_point.begin());

    cout << "Function:" << endl << targetFunction << endl;
    cout << "EPSILON: " << EPSILON << endl;
    cout << "STEP: " << STEP << endl;
    cout << "press enter to start calculating" << endl;

    logFile << "Function:" << endl << targetFunction << endl;
    logFile << "EPSILON: " << EPSILON << endl;
    logFile << "STEP: " << STEP << endl;
    
    cin.get();

    while (true){
        i++;
    
        gradientValue = evaluate2list(gradient, old_point);
        new_point = old_point - STEP * gradientValue;

        logging();

        if (calculateNorm(new_point - old_point) < EPSILON) break;
        if (calculateNorm(gradientValue) < EPSILON) break;
        if (abs(targetFunction.evaluate(old_point) - targetFunction.evaluate(new_point)) < EPSILON) break;

        copy(new_point.begin(), new_point.end(), old_point.begin());
        
    }

    cout << endl << "END" << endl;
    logFile.close(); 
    return 0;
}