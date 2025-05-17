/*
Алгоритм метода Ньютона для поиска минимума функции:
1. Вычисляется матрица Гессе (вторых производных) целевой функции
2. Вычисляется антиградиент функции в текущей точке
3. Выполняется разложение Холецкого матрицы Гессе для решения системы линейных уравнений
4. На основе полученных данных находится следующая точка приближения к минимуму функции
*/


#include "../libs/math.cpp"
#include <fstream>
#include <iostream>

using namespace std;

ofstream logFile("output.log");

int dimension = 2;
Poly targetFunction(dimension, 2, 0);
vector <double> point(dimension, 0);

vector <double> getNextPointByNewton(vector <double> point, Poly function){
    vector <vector <double>> hess(dimension, vector <double>(dimension, 0));
    for (size_t i = 0; i < dimension; i++)
    {
        for (size_t j = i; j < dimension; j++)
        {
            hess[i][j] = derivate(derivate(targetFunction, i), j).evaluate(point);
            hess[j][i] = hess[i][j];
        }
    }

    vector <double> antigradient = (-1) * evaluate2list(Nabla(targetFunction), point);

    cout << "\nantigradient:\n" << antigradient << endl;
    cout << "\npoint:\n" << point << endl;
    cout << "\nhess matrix:\n" << hess << endl;

    logFile << "\nantigradient:\n" << antigradient << endl;
    logFile << "\npoint:\n" << point << endl;
    logFile << "\nhess matrix:\n" << hess << endl;

    // Cholesky decomposition
    vector <vector <double>> L(dimension, vector <double>(dimension, 0));
    for (size_t i = 0; i < dimension; i++)
    {
        for (size_t j = i; j < dimension; j++)
        {
            double tmp = hess[j][i];
            for (size_t k = 0; k < i; k++)
            {
                tmp -= L[i][k] * L[j][k];
            }
            if (i == j){

                L[i][i] = sqrt(tmp);
            }
            else{
                L[j][i] = tmp / L[i][i];
            }
        }
    }

    vector <double> middleValue(dimension, 0);
    vector <double> delta(dimension, 0);

    for (size_t i = 0; i < dimension; i++)
    {
        double val = antigradient[i];
        middleValue[i] = (val - middleValue * L[i]) / L[i][i];
    }

    for (int i = dimension - 1; i >= 0; --i)
    {   
        double val = (middleValue)[i];
        delta[i] = (val - delta * transpose(L)[i]) / transpose(L)[i][i];
    }

    vector <double> result = point + delta;

/*     cout << "+++++++++++++++++++++++++" << endl;
    
    cout << "L: \n" << L << endl;
    cout << "middleValue: " << middleValue << endl;
    cout << "L*middleValue: " << L*middleValue << endl;
    cout << "antigradient: " << antigradient << endl;

    cout << "+++++++++++++++++++++++++" << endl;

    cout << "transpose(L): \n" << transpose(L) << endl;
    cout << "nextPoint: " << nextPoint << endl;
    cout << "transpose(L)*nextPoint: " << transpose(L)*nextPoint << endl;
    cout << "middleValue: " << middleValue << endl;

    cout << "+++++++++++++++++++++++++" << endl;

    cout << "hess: \n" << hess << endl;
    cout << "nextPoint: " << nextPoint << endl;
    cout << "hess*nextPoint: " << hess*nextPoint << endl;
    cout << "antigradient: " << antigradient << endl; */

    return result; 
}

int main(){
    targetFunction.setCoefs("2 2 0 100 0 2 1");
    point[1] = 10;
    size_t k = 0;
    cout << "\nfunction:\n" << targetFunction << endl;

    // int dimension = 4;
    // Poly targetFunction(dimension, 2, 0);
    // targetFunction.setCoefs("10 2 0 0 0 2 0 2 0 0 2.5 0 0 2 0 3 0 0 0 2 3.5 1 1 0 0 1 1 0 1 0 1 1 0 0 1 1 0 1 1 0 1 0 1 0 1 1 0 0 1 1 1");
    // vector <double> point(dimension, 5);

    while(targetFunction.evaluate(point) > 0)
    {
        vector <double> nextPoint = getNextPointByNewton(point, targetFunction);
        
        cout << "iter: " << k << endl;
        cout << "old point:\n" << point << endl;
        cout << "old value: " << targetFunction.evaluate(point) << endl << endl;
        cout << "new point:\n" << nextPoint << endl;
        cout << "new value: " << targetFunction.evaluate(nextPoint) << endl;
        cout << "++++++++++++++++++++++++++++++++++++" << endl;

        logFile << "iter: " << k << endl;
        logFile << "old point:\n" << point << endl;
        logFile << "old value: " << targetFunction.evaluate(point) << endl << endl;
        logFile << "new point:\n" << nextPoint << endl;
        logFile << "new value: " << targetFunction.evaluate(nextPoint) << endl;
        logFile << "++++++++++++++++++++++++++++++++++++" << endl;
        

        copy(nextPoint.begin(), nextPoint.end(), point.begin());
        k++;
    }



    return 0;
}