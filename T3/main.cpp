#include <iostream>
#include <fstream>
#include "../T5/libs/math.cpp"
#include <functional>

using namespace std;

pair<vector<double>, vector<double>> Euler(

    function<double(double, double)> equation,
    double ragngeStart,
    double rangeEnd,
    double startValue,
    double sizeOfCrushing)

{
    pair<vector<double>, vector<double>> res;
    double step = (rangeEnd - ragngeStart) / sizeOfCrushing;
    vector<double> functionArgumets(sizeOfCrushing + 1, 0);
    vector<double> functionValues(sizeOfCrushing + 1, 0);

    functionValues[0] = startValue;
    functionArgumets[0] = ragngeStart;

    for (size_t i = 1; i <= sizeOfCrushing; i++)
    {
        functionArgumets[i] = functionArgumets[i - 1] + step;
        functionValues[i] = functionValues[i - 1] + step * (equation)(functionValues[i - 1], functionArgumets[i - 1]);
    }

    res.first = functionArgumets;
    res.second = functionValues;

    return res;
}

pair<vector<double>, vector<double>> Euler_Cauchy(

    function<double(double, double)> equation,
    double ragngeStart,
    double rangeEnd,
    double startValue,
    double sizeOfCrushing)

{
    pair<vector<double>, vector<double>> res;
    double step = (rangeEnd - ragngeStart) / sizeOfCrushing;
    vector<double> functionArgumets(sizeOfCrushing + 1, 0);
    vector<double> functionValues(sizeOfCrushing + 1, 0);

    functionValues[0] = startValue;
    functionArgumets[0] = ragngeStart;

    for (size_t i = 1; i <= sizeOfCrushing; i++)
    {
        functionArgumets[i] = functionArgumets[i - 1] + step;

        functionValues[i] = functionValues[i - 1] + step * equation(functionValues[i - 1], functionArgumets[i - 1]);

        functionValues[i] = functionValues[i - 1] + 0.5 * step *
                                                        (equation(functionArgumets[i - 1], functionValues[i - 1]) +
                                                         equation(functionArgumets[i], functionValues[i]));
    }

    res.first = functionArgumets;
    res.second = functionValues;

    return res;
}

pair<vector<double>, vector<double>> Euler_improved(

    function<double(double, double)> equation,
    double ragngeStart,
    double rangeEnd,
    double startValue,
    double sizeOfCrushing)

{
    double tmp = 0;

    pair<vector<double>, vector<double>> res;
    double step = (rangeEnd - ragngeStart) / sizeOfCrushing;
    vector<double> functionArgumets(sizeOfCrushing + 1, 0);
    vector<double> functionValues(sizeOfCrushing + 1, 0);

    functionValues[0] = startValue;
    functionArgumets[0] = ragngeStart;

    for (size_t i = 1; i <= sizeOfCrushing; i++)
    {
        functionArgumets[i] = functionArgumets[i - 1] + step;

        tmp = functionValues[i - 1] + 0.5 * step * equation(functionValues[i - 1], functionArgumets[i - 1]);

        functionValues[i] = functionValues[i - 1] + step * equation(functionArgumets[i] - 0.5 * step, tmp);
    }

    res.first = functionArgumets;
    res.second = functionValues;

    return res;
}

int main()
{
    auto myEquation = [](double x, double y)
    { return pow(x + y, 2); };
    pair<vector<double>, vector<double>> result;
    result = Euler(myEquation, 0, 0.5, 0, 5);
    cout << result.first << endl;
    cout << result.second << endl;

    cout << "++++++++++++++++++++++++" << endl;

    result = Euler_Cauchy(myEquation, 0, 0.5, 0, 5);
    cout << result.first << endl;
    cout << result.second << endl;

    cout << "++++++++++++++++++++++++" << endl;

    result = Euler_improved(myEquation, 0, 0.5, 0, 5);
    cout << result.first << endl;
    cout << result.second << endl;

    return 0;
}
