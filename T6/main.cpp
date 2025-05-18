#include "../T5/libs/math.cpp"
#include <iostream>
#include <fstream>
#include <functional>

using namespace std;

// x1 - 2 * x2^2 + 4 * x2 --> max
// - 3x1 - 2x2  = 6

size_t dimension = 2;
size_t equalLimitsCount = 1;
size_t inequalLimitsCount = 0;

const double EPSILON = pow(2, -60);
const double MINSTEP = pow(2, -60);
const double MAXSTEP = pow(2, 0);
const double C = pow(2, 0.5);
double r = pow(2, 0);
double k = 0;

Poly targetFunction(dimension, 2, 0);
Poly extraFunction(dimension, 2, 0);
vector <Poly> equalLimits(equalLimitsCount, Poly(dimension, 2, 0));
vector <Poly> inequalLimits(inequalLimitsCount, Poly(dimension, 2, 0));


vector<double> annFunctGradien(function<double(vector <double>)> targetFunction, vector <double> point){
    vector <double> res(dimension, 0);
    vector <double> prevPoint(dimension, 0);
    vector <double> nextPoint(dimension, 0);
    const double STEP = pow(10, -6);
    for (size_t i=0; i<dimension; i++){
        copy(point.begin(), point.end(), prevPoint.begin());
        copy(point.begin(), point.end(), nextPoint.begin());
        
        prevPoint[i] = prevPoint[i] - STEP;
        nextPoint[i] = nextPoint[i] + STEP;
        // cout << targetFunction(prevPoint) << endl;
        // cout << targetFunction(nextPoint) << endl;
        // cout << STEP << endl;
        // cout << (targetFunction(nextPoint) - targetFunction(prevPoint))/(2*STEP) << endl;
        // cin.get();
        res[i] = (targetFunction(nextPoint) - targetFunction(prevPoint))/(2*STEP);
    }
    return res;
}

double goldenSectionSearch(double a, double b, double epsilon, function <double(double)> func) {
    const double PHI = (sqrt(5) - 1) / 2;

    double var1 = b - PHI * (b - a);
    double var2 = a + PHI * (b - a);
    double y1 = func(var1);
    double y2 = func(var2);
    

    while ((b - a) > epsilon) {
        if (y1 < y2) {
            b = var2;
            var2 = var1;
            y2 = y1;

            var1 = b - PHI * (b - a);
            y1 = func(var1);
        } else {
            a = var1;
            var1 = var2;
            y1 = y2;

            var2 = a + PHI * (b - a);
            y2 = func(var2);
        }
    }

    double result = (a + b) / 2;
    return result;
}

double Newton(function <double(double)> func, double startArg, double minValue, double maxValue){
    double der1 = 0;
    double der2 = 0;

    double arg = startArg;
    double delta;
    const double STEP = pow(10, -6);
    const double EPSILON = pow(10, -10);
    
 
    cout << "---------------- Newton -----------------\n";
    cout << "EPSILON: " << EPSILON << endl;
    cout << "STEP: " << STEP << endl;
    cout << endl;
    while (true) {
        
        der1 = (func(arg + STEP) - func(arg - STEP))/(2 * STEP);
        der2 = (func(arg + STEP) - 2 * func(arg) + func(arg - STEP)) / pow(STEP, 2);

        cout << "arg: " << arg << endl;
        cout << "value: " << func(arg) << endl;
        cout << "der1: " << der1 << endl;
        cout << "der2: " << der2 << endl;

        if (der2 == 0) {cout << "NEWTON.DER2.NULL" << endl; break;}

        delta = der1 / abs(der2);

        cout << "val del: " << (func(arg - delta) - func(arg)) << endl;

        arg -= delta;
        
        cout << "delta: " << delta << endl;
        cout << endl;
        // cin.get();

        if (arg > maxValue) {cout << "NEWTON.maxValue" << endl; arg = maxValue; break;}
        
        if (arg < minValue) {cout << "NEWTON.minValue" << endl; arg = minValue; break;}
        
        if (abs(delta) < EPSILON) break;
        if (abs(func(arg - delta) - func(arg)) < EPSILON) break;
    }
    cout << "---------------------------------\n";
    return arg;
}


double getOptimalStep(function<double(vector <double>)> targetFunction, vector <double> point){
    auto oneDimFunction {[point, targetFunction](double step) {return targetFunction(point - step * annFunctGradien(targetFunction, point));} };

    double startPoint = goldenSectionSearch(MINSTEP, MAXSTEP, pow(2, -5), oneDimFunction);
    double step = Newton(oneDimFunction, startPoint, MINSTEP, MAXSTEP);

    return step;
}

vector <double> fastGradDest(function<double(vector <double>)> targetFunction, vector <double> startPoint){
    vector <double> nextPoint(dimension, 0);
    vector <double> prevPoint(dimension, 0);

    double step = 0;
    const double EPSILON = pow(10, -10);

    copy(startPoint.begin(), startPoint.end(), prevPoint.begin());
    copy(startPoint.begin(), startPoint.end(), nextPoint.begin());

    cout << "++++++ Gradient destination start ++++" << endl;

    while (true) {
        step = getOptimalStep(targetFunction, prevPoint);
        
        nextPoint = prevPoint - step * annFunctGradien(targetFunction, prevPoint);

        cout << "step: " << step << endl;
        cout << "nextPoint: " << nextPoint << endl;
        cout << "grad: " << annFunctGradien(targetFunction, prevPoint) << endl;
        cout << "value: " << targetFunction(nextPoint) << endl;
        cout << "value delta: " << targetFunction(prevPoint) - targetFunction(nextPoint) << endl;
        cout << endl;

        if (calculateNorm(prevPoint - nextPoint) < EPSILON) break;
        if (abs(targetFunction(prevPoint) - targetFunction(nextPoint)) < EPSILON) break;
        if (calculateNorm(annFunctGradien(targetFunction, nextPoint)) < EPSILON) break;

        copy(nextPoint.begin(), nextPoint.end(), prevPoint.begin());
        // cin.get();
    }
    cout << "++++++++++++++++++++++++++++++++++++++++++++" << endl;

    return nextPoint;
}

double newFuntion(vector <double> point){
    double res = 0;

    res -= targetFunction.evaluate(point);
    res += (r/2) * extraFunction.evaluate(point);

    for (size_t i = 0; i < inequalLimitsCount; i++)
    {
        res += (r/2) * max(.0, 
                           pow(inequalLimits[i], 2).evaluate(point));
    }
    return res;
}


int main(){
    targetFunction.setCoefs("3 1 0 1 0 2 -2 0 1 4");
    equalLimits[0].setCoefs("3 1 0 -3 0 1 -2 0 0 -6");
    
    vector<double> startPoint(dimension, 0.1);
    vector<double> oldPoint(dimension, 0);
    vector<double> newPoint(dimension, 0);

    
    for (size_t i=0; i<equalLimitsCount; i++){
        extraFunction = extraFunction + pow(equalLimits[i], 2);
    }

    copy(startPoint.begin(), startPoint.end(), oldPoint.begin());

    cout << "function: " << targetFunction << endl;
    cout << "limit: " << equalLimits[0] << endl;
    cout << "MINSTEP: " << MINSTEP << endl;
    
    while (true){
        newPoint = fastGradDest(newFuntion, oldPoint);
        cout << endl;
        cout << " == global iter: " << k << " ==\n";
        cout << "new point:" << newPoint << endl;
        cout << "r: " << r << endl;
        cout << "value: " << targetFunction.evaluate(newPoint) << endl;
        cout << "extra value: " << newFuntion(newPoint) << endl;
        cout << "P: " << newFuntion(newPoint) + targetFunction.evaluate(newPoint) << endl;
        cout << endl; 
        // cin.get();
        if ((newFuntion(newPoint) + targetFunction.evaluate(newPoint)) <= EPSILON) break;

        r = C * r;
        copy(newPoint.begin(), newPoint.end(), oldPoint.begin());
        k++;
    }
    
    return 0;
}