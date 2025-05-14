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


size_t i;

Poly getDirectionFunction(vector<double> point){
    Poly function(1, 0, 0);
    vector<Poly> extraFunction(dimension, Poly(1, 0, 0));

    for (size_t i = 0; i < dimension; i++)
    {
        extraFunction[i] = Poly(0, 0, point[i]) - gradientValue[i] * (Poly(1, 1, 1) - Poly(1, 0, 1));
    }

    cout << gradient << endl;
    cout << gradientValue << endl;
    cout << extraFunction << endl;

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

void logging(){
    i++;

    cout << endl << "== Iter: " << i << " ==" << endl;
    cout << "step: " << step << endl;
    cout << "new point: " << new_point << endl;
    cout << "points delta: " << calculateNorm(new_point - old_point) << endl;
    cout << "value: " << targetFunction.evaluate(new_point) << endl;
    cout << "gradientValue: " << gradientValue << endl;
    cout << "gradientNotm: " << calculateNorm(gradientValue) << endl;

    logFile << endl << "== Iter: " << i << " ==" << endl;
    logFile << "step: " << step << endl;
    logFile << "new point: " << new_point << endl;
    logFile << "points delta: " << calculateNorm(new_point - old_point) << endl;
    logFile << "value: " << targetFunction.evaluate(new_point) << endl;
    logFile << "gradientValue: " << gradientValue << endl;
    logFile << "gradientNotm: " << calculateNorm(gradientValue) << endl;
}



int main(){
    targetFunction.setCoefs("6 3 0 1 1 1 -1 0 2 1 1 0 -2 0 1 3 0 0 -4");
    gradient = Nabla(targetFunction);
    copy(start_point.begin(), start_point.end(), old_point.begin());
    copy(start_point.begin(), start_point.end(), new_point.begin());

    cout << "Function:" << endl << targetFunction << endl;
    cout << "EPSILON: " << EPSILON << endl;
    cout << "MIN STEP: " << MIN_STEP << endl;
    cout << "MAX STEP: " << MAX_STEP << endl;
    cout << "press enter to start calculating" << endl;

    logFile << "Function:" << endl << targetFunction << endl;
    logFile << "EPSILON: " << EPSILON << endl;
    logFile << "MIN STEP: " << MIN_STEP << endl;
    logFile << "MAX STEP: " << MAX_STEP << endl;
    

    while (true){
        gradientValue = evaluate2list(gradient, old_point);
        step = getOptimalStep(old_point);
        new_point = old_point - step * gradientValue;

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