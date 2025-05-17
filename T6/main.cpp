#include "../T5/libs/math.cpp"
#include <iostream>
#include <fstream>

using namespace std;

// x1 - 2 * x2^2 + 4 * x2 --> max
// - 3x1 - 2x2  = 6

size_t dimension = 2;
size_t equalLimitsCount = 1;
size_t inequalLimitsCount = 0;

const double EPSILON = pow(2, -60);
const double MINSTEP = pow(2, -60);
const double MAXSTEP = pow(2, 0);
const double C = pow(2, 1);
double r = pow(2, 0);
double k = 0;

Poly targetFunction(dimension, 2, 0);
Poly extraFunction(dimension, 2, 0);
vector <Poly> equalLimits(equalLimitsCount, Poly(dimension, 2, 0));
vector <Poly> inequalLimits(inequalLimitsCount, Poly(dimension, 2, 0));


double newFuntion(vector <double> point){
    double res = 0;

    res += targetFunction.evaluate(point);
    res += (r/2) * extraFunction.evaluate(point);

    for (size_t i = 0; i < inequalLimitsCount; i++)
    {
        res += (r/2) * max(.0, 
                           pow(inequalLimits[i], 2).evaluate(point));
    }
    return res;
}

vector <double>newFuntionGradient(vector <double> point){
    vector<double> res(dimension, 0);

    res = res + evaluate2list(Nabla(targetFunction), point);
    res = res + (r/2) * evaluate2list(Nabla(extraFunction), point);

    for (size_t j = 0; j < inequalLimitsCount; j++)
    {
        vector <double> val = evaluate2list(
                                        Nabla(
                                            pow(inequalLimits[j], 2)), point);

        for (size_t ind = 0; ind < dimension; ind++)
        {
            res[ind] = res[ind] + (r/2) * max(.0, val[ind]);
        }   
    }
    return res;
}


double getOptimalStep(vector<double> point){
    const double t = pow(2, -10);
    const double dStep = pow(2, -50);
    
    double oldStep = MINSTEP;
    double newStep = 0;
    double derivate = 0;
    vector <double> arg1(dimension, 0);
    vector <double> arg2(dimension, 0);
    
    while (true) {
        arg1 = point - oldStep * newFuntionGradient(point);
        arg2 = point - (oldStep + dStep) * newFuntionGradient(point);

        derivate = (newFuntion(arg2) - newFuntion(arg1))/dStep;
        
        newStep = oldStep - t * derivate;

        // cout << arg1 << "|||" << arg2 << endl;
        // cout << newFuntion(arg2) - newFuntion(arg1) << endl; 
        // cout << abs(oldStep - newStep) << endl;
        // cout << abs(newFuntion(point - oldStep * newFuntionGradient(point)) - newFuntion(point - newStep * newFuntionGradient(point))) << endl;
        // cout << abs(derivate) << endl;

        if (abs(oldStep - newStep) < EPSILON) break;
        if (abs(newFuntion(point - oldStep * newFuntionGradient(point)) - newFuntion(point - newStep * newFuntionGradient(point))) < EPSILON) break;
        if (abs(derivate) < EPSILON) break;
        
        oldStep = newStep;

        // cin.get();
    }
    return newStep;
}

vector <double>fastGradientDest(vector <double> startPoint){
    vector <double> nextPoint(dimension, 0);
    double step = 0;
    size_t iter = 1;

    while (true)
    {
        step = getOptimalStep(startPoint);
        nextPoint = startPoint - step * newFuntionGradient(startPoint);
        
        cout << endl << "Grad destination iteration: " << iter << endl;
        cout << "Step: " << step << endl;
        cout << "New point: " << nextPoint << endl;
        cout << "extra value: " << newFuntion(nextPoint) << endl;
        cout << "function value: " << (-1) * targetFunction.evaluate(nextPoint) << endl;
        cout << "points delta: " << calculateNorm(nextPoint - startPoint) << endl;
        cout << "extra value delta: " << newFuntion(nextPoint) - newFuntion(startPoint) << endl; 
        
        if (newFuntion(nextPoint) < EPSILON) break;
        if (calculateNorm(newFuntionGradient(nextPoint)) < EPSILON) break;
        if (calculateNorm(nextPoint - startPoint) < EPSILON) break;

        copy(nextPoint.begin(), nextPoint.end(), startPoint.begin());
        iter++;
        // cin.get();
    }

    return nextPoint;
    
}
int main(){
    double gradStep = pow(2, -5);
    double k = 1;
    // x1 - 2 * x2^2 + 4 * x2
    targetFunction.setCoefs("3 1 0 1 0 2 -2 0 1 4");
    targetFunction = (-1) * targetFunction; //maximisation
    equalLimits[0].setCoefs("3 1 0 -3 0 1 -2 0 0 -6");
    
    vector<double> startPoint(dimension, 0);
    vector<double> oldPoint(dimension, 0);
    vector<double> newPoint(dimension, 0);

    
    for (size_t i=0; i<equalLimitsCount; i++){
        extraFunction = extraFunction + pow(equalLimits[i], 2);
    }

    copy(startPoint.begin(), startPoint.end(), oldPoint.begin());

    cout << "function: " << targetFunction << endl;
    cout << "limit: " << equalLimits[0] << endl;
    cout << "MINSTEP: " << MINSTEP << endl;
    
    cin.get();

    while (true){
        newPoint = fastGradientDest(oldPoint);
        cout << endl;
        cout << " == global iter: " << k << " ==\n";
        cout << "new point:" << newPoint << endl;
        cout << "r: " << r << endl;
        cout << "value: " << (-1) * targetFunction.evaluate(newPoint) << endl;
        cout << "P: " << newFuntion(newPoint) - targetFunction.evaluate(newPoint) << endl;
        cout << endl; 
        if (newFuntion(newPoint) - targetFunction.evaluate(newPoint) <= EPSILON) break;

        r = C * r;
        copy(newPoint.begin(), newPoint.end(), oldPoint.begin());
        k++;
    }
    
    return 0;
}

