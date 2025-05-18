#include <iostream>
#include <cmath>
#include <functional>

using namespace std;

const double EPSILON = pow(2, -60);
const double MINSTEP = pow(2, -60);
const double MAXSTEP = pow(2, 0);

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
    while (true) {
        
        der1 = (func(arg + STEP) - func(arg - STEP))/(2 * STEP);
        der2 = (func(arg + STEP) - 2 * func(arg) + func(arg - STEP)) / pow(STEP, 2);

        cout << "arg: " << arg << endl;
        cout << "value: " << func(arg) << endl;
        cout << "der1: " << der1 << endl;
        cout << "der2: " << der2 << endl;

        if (der2 == 0) {cout << "NEWTON.DER2.NULL" << endl; break;}

        cout << "val del: " << (func(arg + STEP) - 2 * func(arg) + func(arg - STEP)) << endl;

        delta = der1 / abs(der2);
        arg -= delta;
        
        cout << "delta: " << delta << endl;
        cout << endl;
        cin.get();

        if (arg > maxValue) {cout << "NEWTON.maxValue" << endl; arg = maxValue; break;}
        
        if (arg < minValue) {cout << "NEWTON.minValue" << endl; arg = minValue; break;}
        
        if (abs(delta) < EPSILON) break;
    }
    return arg;
}

int main(){
    auto func {[](double x){return sin(x) * (pow(x, 2) + sqrt(x));} };
    double startPoint = 2.5;
    double optimum = Newton(func, startPoint, 0, 10);
    cout << optimum << "; " << func(optimum) << endl;

    return 0;
}