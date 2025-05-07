// Подключение необходимых заголовочных файлов
#include <cmath>
#include <algorithm>
#include <numeric>
#include "Poly.cpp"

// Оператор сложения двух векторов чисел
std::vector<double> operator+(std::vector<double> self, std::vector<double> other) {
    if (size(self) == size(other)) {
        std::vector<double> res(size(self), 0);
        for (int i = 0; i < size(self); i++) {
            res[i] = self[i] + other[i];
        }
        return res;
    } else {
        std::cerr << "Размерности различны" << std::endl;
        return std::vector<double>();
    }
}

// Оператор сложения двух векторов полиномов
std::vector<Poly> operator+(std::vector<Poly> self, std::vector<Poly> other) {
    if (size(self) == size(other)) {
        std::vector<Poly> res(size(self), 0);
        for (int i = 0; i < size(self); i++) {
            res[i] = self[i] + other[i];
        }
        return res;
    } else {
        std::cerr << "Размерности различны" << std::endl;
        return std::vector<Poly>();
    }
}

// Оператор умножения вектора чисел на число
std::vector<double> operator*(double num, std::vector<double> self) {
    std::vector<double> res(size(self), 0);
    for (int i = 0; i < size(self); i++) {
        res[i] = self[i] * num;
    }
    return res;
}

// Оператор умножения вектора полиномов на число
std::vector<Poly> operator*(double num, std::vector<Poly> self) {
    std::vector<Poly> res(size(self), 0);
    for (int i = 0; i < size(self); i++) {
        res[i] = num * self[i];
    }
    return res;
}

// Оператор умножения вектора полиномов на полином
std::vector<Poly> operator*(Poly pol, std::vector<Poly> self) {
    std::vector<Poly> res(size(self), 0);
    for (int i = 0; i < size(self); i++) {
        res[i] = pol * self[i];
    }
    return res;
}

// Оператор вычитания двух векторов чисел
std::vector<double> operator-(std::vector<double> self, std::vector<double> other) {
    return self + (-1 * other);
}

// Оператор вычитания двух векторов полиномов
std::vector<Poly> operator-(std::vector<Poly> self, std::vector<Poly> other) {
    return self + (-1 * other);
}


// Оператор скалярного произведения векторов чисел
double operator*(std::vector<double> self, std::vector<double>other){
    if (self.size() == other.size()){
        double res;
        for (size_t i=0; i<self.size(); ++i){
            res = res + self[i] * other[i];
        }
        return res;
    }
    else{
        std::cerr << "Размерности должны совпадать";
        return (double) 0.0;
    }
    
}

// Оператор скалярного произведения векторов полиномов
Poly operator*(std::vector<Poly> self, std::vector<Poly>other){
    if (self.size() == other.size()){
        Poly res;
        for (size_t i=0; i<self.size(); ++i){
            res = res + self[i] * other[i];
        }
        return res;
    }
    else{
        std::cerr << "Размерности должны совпадать";
        return Poly();
    }
    
}

// Оператор вывода вектора полиномов
std::ostream& operator<<(std::ostream& os, std::vector<Poly> self) {
    for (int i = 0; i < self.size(); i++) {
        os << self[i] << ", ";
    }
    return os;
}

// Оператор вывода вектора целых чисел
std::ostream& operator<<(std::ostream& os, std::vector<int> self) {
    for (int i = 0; i < self.size(); i++) {
        os << self[i] << ", ";
    }
    return os;
}

// Оператор вывода вектора вещественных чисел
std::ostream& operator<<(std::ostream& os, std::vector<double> self) {
    for (int i = 0; i < self.size(); i++) {
        os << self[i] << ", ";
    }
    return os;
}

// Вычисление евклидовой нормы вектора
double calculateNorm(const std::vector<double>& vec) {
    double sumSquares = 0.0;
    for (const auto& component : vec) {
        sumSquares += component * component;
    }
    return sqrt(sumSquares);
}

std::vector<std::vector<size_t>> permutations(unsigned int N) {
    std::vector<size_t> perm(N);
    std::iota(perm.begin(), perm.end(), 0);

    std::vector<std::vector<size_t>> result;
    result.reserve(1);
    
    do {
        result.push_back(perm);
    } while (std::next_permutation(perm.begin(), perm.end()));

    return result;
}

bool isEvenPermutation(const std::vector<size_t>& perm) {
    size_t invCount = 0;
    for (size_t i = 0; i < perm.size(); ++i) {
        for (size_t j = i + 1; j < perm.size(); ++j) {
            if (perm[i] > perm[j]) {
                ++invCount;
            }
        }
    }
    return (invCount % 2) == 0;
}

Poly Hessian(const Poly& function){
    size_t dim = function.getVarc();
    std::vector <std::vector <size_t>> indices = permutations(dim);
    Poly res(dim, 0, 0);
    for (std::vector <size_t> perm: indices){
        Poly tmp = pow(Poly(dim, 0, -1), (float) isEvenPermutation(perm));
        for (size_t i = 0; i < dim; i++)
        {
            tmp = tmp * derivate(derivate(function, i), perm[i]);
        }
        res = res + tmp;
    }
    return res;
}