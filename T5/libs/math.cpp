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

std::vector <double> operator*(const std::vector<std::vector<double>>& matrix, const std::vector<double>& vector){
    if (matrix[0].size() != vector.size()){
        std::cerr << "Длина вектора должна совпадать с шириной матрицы";
        return std::vector<double>(0, 0);
    }
    std::vector<double> res(matrix.size(), 0);
    for (size_t i=0; i<res.size(); ++i){
        for (size_t k=0; k<vector.size(); ++k){
            res[i] += vector[k] * matrix[i][k];
        }
    }
    return res;
}

std::vector <Poly> operator*(const std::vector<std::vector<double>>& matrix, const std::vector<Poly>& vector){
    if (matrix[0].size() != vector.size()){
        std::cerr << "Длина вектора должна совпадать с шириной матрицы";
        return std::vector<Poly>(0, Poly());
    }
    std::vector<Poly> res(matrix.size(), 0);
    for (size_t i=0; i<res.size(); ++i){
        for (size_t k=0; k<vector.size(); ++k){
            res[i] = res[i] + matrix[i][k] * vector[i];
        }
    }
    return res;
}

std::vector<std::vector<double>> vector2matrix(const std::vector<double>& vector){
    size_t dim = vector.size();
    std::vector<std::vector<double>> res(dim, std::vector<double>(dim, 0));
    for (size_t i = 0; i < dim; i++)
    {
        for (size_t j = i; j < dim; j++)
        {
            res[i][j] = vector[i] * vector[j];
            res[j][i] = res[i][j];
        }
    }
    return res;
}

std::vector<std::vector<double>> operator*(const double& scalar, const std::vector<std::vector<double>>& matrix){
    size_t dim1 = matrix.size();
    size_t dim2 = matrix[0].size();
    std::vector<std::vector<double>> res(dim1, std::vector<double> (dim2, 0));
    for (size_t i = 0; i < dim1; i++)
    {
        for (size_t j = i; j < dim2; j++)
        {
            res[i][j] = scalar * matrix[i][j];
        }
    }
    return res;
}

std::vector<std::vector<double>> operator+(const std::vector<std::vector<double>>& lhs, const std::vector<std::vector<double>>& rhs){
    if ((lhs.size() != rhs.size()) || (rhs[0].size() != lhs[0].size())){
        std::cerr << "Размеры матриц должн совпадать";
        return std::vector<std::vector<double>>(lhs.size(), std::vector<double>(lhs[0].size(), 0));
    }
    std::vector<std::vector<double>> res(lhs.size(), std::vector<double>(lhs[0].size(), 0));
    for (size_t i = 0; i < lhs.size(); i++)
    {
        for (size_t j = i; j < rhs.size(); j++)
        {
            res[i][j] = lhs[i][j] + rhs[i][j];
        }
    }
    return res;
}

std::vector<std::vector<double>> operator-(const std::vector<std::vector<double>>& lhs, const std::vector<std::vector<double>>& rhs){
    return lhs + (-1) * rhs;
}

std::vector<std::vector<double>> operator*(const std::vector<std::vector<double>>& lhs, const std::vector<std::vector<double>>& rhs){
    if (lhs[0].size() != rhs.size()){
        std::cerr << "Высота правой матрицы должна совпадать с шириной левой матрицы";
        return std::vector<std::vector<double>>(lhs.size(), std::vector<double>(rhs[0].size(), 0));
    }
    std::vector<std::vector<double>> res(lhs.size(), std::vector<double>(rhs[0].size(), 0));
    for (size_t i = 0; i < lhs.size(); i++)
    {
        for (size_t j = i; j < rhs[0].size(); j++)
        {
            for (size_t k = 0; k < lhs[0].size(); k++)
            {
                res[i][j] += lhs[i][k] * rhs[k][j];
            }
        }
    }
    return res;
}

std::ostream& operator<<(std::ostream& os, const std::vector<std::vector<double>>& self) {
    for (int i = 0; i < self.size(); i++) {
        os << self[i] << std::endl;
    }
    return os;
}