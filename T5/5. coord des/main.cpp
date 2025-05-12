#include "../libs/math.cpp"
#include <fstream>
#include <iostream>

using namespace std;

ofstream logFile("output.log");

int dimension = 2;
Poly targetFunction(dimension, 3, 0);
vector<double> new_point(dimension, 0);
vector<double> old_point(dimension, 0);

const double EPSILON = pow(2, -60);
const double STEP = pow(2, -5);

int main(){
    bool flag = false;

    targetFunction.setCoefs("6 3 0 1 1 1 -1 0 2 1 1 0 -2 0 1 3 0 0 -4");
    vector<Poly> gradient = Nabla(targetFunction);

    cout << "Начало поиска минимума функции" << endl;
    logFile << "Начало поиска минимума функции" << endl;

    double k=1;

    while (true){
        for (size_t i = 0; i < dimension; i++)
        {
            cout << "Итерация: " << k << endl;
            logFile << "Итерация: " << k << endl;

            copy(new_point.begin(), new_point.end(), old_point.begin());
            new_point[i] = old_point[i] - STEP * gradient[i].evaluate(old_point);
            
            cout << "Итерация по координате " << i + 1 << ":" << endl;
            cout << "Старая точка: " << old_point << endl;
            cout << "Новая точка: " << new_point << endl;
            cout << "Значение: " << targetFunction.evaluate(new_point) << endl;
            cout << "Изменение значения: " << abs(targetFunction.evaluate(new_point) - targetFunction.evaluate(old_point)) << endl;
            
            logFile << "Итерация по координате " << i + 1 << ":" << endl;
            logFile << "Старая точка: " << old_point << endl;
            logFile << "Новая точка: " << new_point << endl;
            logFile << "Значение: " << targetFunction.evaluate(new_point) << endl;
            logFile << "Изменение значения: " << abs(targetFunction.evaluate(new_point) - targetFunction.evaluate(old_point)) << endl;
            
            k++;
            cout << endl;
            logFile << endl;

            if (calculateNorm(old_point - new_point) < EPSILON){
                cout << "Достигнута требуемая точность по норме разности точек" << endl;
                logFile << "Достигнута требуемая точность по норме разности точек" << endl;
                flag = true;
                break;
            }
            if (abs(targetFunction.evaluate(old_point) - targetFunction.evaluate(new_point)) < EPSILON){
                cout << "Достигнута требуемая точность по значению функции" << endl;
                logFile << "Достигнута требуемая точность по значению функции" << endl;
                flag = true;
                break;
            }
        }
        if (flag) break;
    }

    cout << "Найден минимум в точке: " << new_point << endl;
    cout << "Значение функции в минимуме: " << targetFunction.evaluate(new_point) << endl;
    
    logFile << "Найден минимум в точке: " << new_point << endl;
    logFile << "Значение функции в минимуме: " << targetFunction.evaluate(new_point) << endl;

    logFile.close();

    return 0;
}