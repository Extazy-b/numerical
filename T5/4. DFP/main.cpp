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
vector <double> direction(dimension, 0);

vector<vector<double>> hess(dimension, vector <double> (dimension, 0));
vector <double> P(dimension, 0);
vector <double> q(dimension, 0);
vector <vector <double>> hessDelta(dimension, vector<double>(dimension, 0));

double step;
size_t k;

Poly getDirectionFunction(vector<double> point){
    Poly function(1, 0, 0);
    vector<Poly> extraFunction(dimension, Poly(1, 0, 0));

    for (size_t i = 0; i < dimension; i++)
    {
        extraFunction[i] = Poly(0, 0, point[i]) + direction[i] * (Poly(1, 1, 1) - Poly(1, 0, 1));
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
    logFile << "\tP: " << P << endl;
    logFile << "\tQ: " << q << endl;
    logFile << "\thess delta: " << hessDelta << endl;
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
    for (size_t i = 0; i < dimension; i++) hess[i][i] = 1;
    
    k = 1;

    cout << "Function:" << endl << targetFunction << endl;
    cout << "EPSILON: " << EPSILON << endl;
    cout << "MIN STEP: " << MIN_STEP << endl;
    cout << "MAX STEP: " << MAX_STEP << endl;
    cout << "press enter to start calculating" << endl;

    logFile << "Function:" << endl << targetFunction << endl;
    logFile << "EPSILON: " << EPSILON << endl;
    logFile << "MIN STEP: " << MIN_STEP << endl;
    logFile << "MAX STEP: " << MAX_STEP << endl;

    cin.get();

    while (true){  
        cout << endl << "== Iter: " << k << " ==" << endl;
        logFile << endl << "== Iter: " << k << " ==" << endl;
                
        for (size_t j = 0; j<dimension; ++j){
            gradientValue = evaluate2list(gradient, new_point);
            direction = (-1) * (hess * gradientValue);

            step = getOptimalStep(new_point);

            P = step * direction;
            new_point = new_point + P;
            q = evaluate2list(gradient, new_point) - gradientValue;

            hessDelta = ((1/(P*q)) * vector2matrix(P)) - ((1/(q * (hess * q))) * (hess * vector2matrix(q) * hess));
            hess = hess + hessDelta;

            logging(j);
        }

        if (calculateNorm(gradientValue) < EPSILON) break;
        if (calculateNorm(new_point - old_point) < EPSILON) break;
        if (abs(targetFunction.evaluate(new_point) - targetFunction.evaluate(old_point)) < EPSILON) break;

        copy(new_point.begin(), new_point.end(), old_point.begin());
        k++;
    }

    cout << endl << "END" << endl;
    logFile.close(); 
    return 0;
}