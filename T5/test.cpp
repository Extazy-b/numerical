#include <iostream>
#include "math.cpp"
#include <string>
#include <vector>

using namespace std;

//6 3 0 1 1 1 -1 0 2 1 1 0 -2 0 1 3 0 0 -4

int test_grad(){
    Poly poly(2, 3);
    poly.setCoefs("6 3 0 1 1 1 -1 0 2 1 1 0 -2 0 1 3 0 0 -4");
    cout << poly << endl << endl;
    int varc = poly.getVarc();

    vector <Poly> grad = Nabla(poly);
    cout << grad << endl << endl;

    Poly step(varc+1, 1, 0);
    vector<int> ind(varc+1, 0);
    ind[varc] = 1;
    step(ind) = 1;
    cout << step << endl << endl;

    vector<Poly>empty_point(varc, Poly(varc, 1, 0));
    for (int i=0; i<varc; i++){
        ind = vector<int>(varc, 0);
        ind[i] = 1;
        empty_point[i](ind) = 1;
    }
    cout << empty_point << endl << endl;

    vector<Poly> func = empty_point - step*grad;
    cout << func << endl << endl;

    Poly extra_poly = poly.compose(func);
    cout << extra_poly << endl << endl;

    Poly der = derivate(extra_poly, varc);
    cout << der << endl << endl; 

    return 0;
}

int test_sum(){
    // Create first polynomial with 2 variables and degree 3
        Poly p1(2, 3);
        p1.setCoefs("5 0 0 2 0 2 -1 1 1 3 2 1 2 3 0 1");  // x^3 + 2x^2y + 3xy - y^2 + 2
    
        // Create second polynomial with 3 variables and degree 2
        Poly p2(3, 2);
        p2.setCoefs("5 0 0 0 1 0 0 1 1 0 1 1 -1 1 1 0 1 2 0 0 2");  // 2x^2 + xy - yz + z + 1
    
        // Sum the polynomials
        Poly sum = p1.resize(3, 3) + p2.resize(3, 3);
        cout << "First polynomial (2 vars, degree 3): " << p1 << endl;
        cout << "Second polynomial (3 vars, degree 2): " << p2 << endl;
    
        // Output results
        
        cout << p1.getVarc() << endl;
        cout << p2.getDeg() << endl;
        // cout << p1.resize(3, 3)("3 2 1") << " + " << p2.resize(3, 3)("3 2 1") << " = " << sum("3 2 1") << endl;

        // cout << sum("3 2 1") << endl;
        cout << "Sum result: " << sum << endl;
    
        return 0;
    
}

int test_prod(){
    // Test polynomial multiplication
        Poly p1(2, 2);  // First polynomial
        p1.setCoefs("2 1 0 1 0 1 1");  // x + y
        
        Poly p2(1, 1);  // Second polynomial
        p2.setCoefs("2 1 2 0 1");    // 2x + 1
        
        Poly result = pow(p1, 3) * p2;  // Multiply polynomials
        
        cout << "First polynomial: " << p1 << endl;
        cout << "Second polynomial: " << p2 << endl;
        cout << "Result of multiplication: " << result << endl;
        
    
    return 0;
}

int test_compose(){
    Poly p1(2, 2, 1);
    Poly p2(1, 1, 1);
    
    cout << p1 << endl;
    cout << p2 << endl;

    vector <Poly> arg(2, p2);

    Poly p = p1.compose(arg);

    cout << p << endl;

    return 0;
}

int main(){
    test_grad();
    return 0;
}