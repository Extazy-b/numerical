/*
Метод координатного спуска:
1. Начинаем с произвольной точки
2. Поочередно двигаемся вдоль каждой координатной оси
3. В каждом направлении ищем минимум функции с заданным шагом
4. Процесс повторяется, пока не достигнем точки, где изменения меньше заданной точности
*/

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

size_t k;

void logging(){
    cout << endl << "== Iter: " << k << " ==" << endl;
    cout << "new point: " << new_point << endl;
    cout << "points delta: " << calculateNorm(new_point - old_point) << endl;
    cout << "value: " << targetFunction.evaluate(new_point) << endl;

    logFile << endl << "== Iter: " << k << " ==" << endl;
    logFile << "new point: " << new_point << endl;
    logFile << "points delta: " << calculateNorm(new_point - old_point) << endl;
    logFile << "value: " << targetFunction.evaluate(new_point) << endl;
}

int main(){
    targetFunction.setCoefs("6 3 0 1 1 1 -1 0 2 1 1 0 -2 0 1 3 0 0 -4");
    vector<Poly> gradient = Nabla(targetFunction);

    srand(time(0));

    cout << "Function:" << endl << targetFunction << endl;
    cout << "EPSILON: " << EPSILON << endl;
    cout << "STEP: " << STEP << endl;
    cout << "press enter to start calculating" << endl;

    logFile << "Function:" << endl << targetFunction << endl;
    logFile << "EPSILON: " << EPSILON << endl;
    logFile << "STEP: " << STEP << endl;

    cin.get();

    size_t i = 0;
    while (true){  
        // i = (i < dimension - 1) ? i+1 : 0;
        i = rand() % dimension;
        
        new_point[i] = old_point[i] - STEP * gradient[i].evaluate(old_point);
        
        k++;
        
        logging();

        if (calculateNorm(old_point - new_point) < EPSILON) break;
        if (abs(targetFunction.evaluate(old_point) - targetFunction.evaluate(new_point)) < EPSILON) break;

        copy(new_point.begin(), new_point.end(), old_point.begin());
    }

    cout << endl << "END" << endl;
    logFile.close(); 
    return 0;
}