#include <cmath>
#include "Poly.cpp"


std::vector <Poly> operator+(std::vector <Poly> self, std::vector <Poly> other){
    if (size(self) == size(other)){
        std::vector <Poly> res(size(self), 0);
        for (int i=0; i<size(self); i++){
            res[i] = self[i] + other[i];
        }
        return res;
    }
    else {
        std::cerr << "Размерности различны" << std::endl;
        return std::vector<Poly>();
    }
}

std::vector <Poly> operator*(double num, std::vector <Poly> self){
    std::vector <Poly> res(size(self), 0);
    for (int i=0; i<size(self); i++){
        res[i] = num * self[i];
    }
    return res;
}
std::vector <Poly> operator*(Poly pol, std::vector <Poly> self){
    std::vector <Poly> res(size(self), 0);
    for (int i=0; i<size(self); i++){
        res[i] = pol * self[i];
    }
    return res;
}
std::vector <Poly> operator-(std::vector <Poly> self, std::vector <Poly> other){
    return self + (-1 * other);
}

std::ostream& operator<<(std::ostream& os, std::vector <Poly> self){
    for (int i=0; i<self.size(); i++){
        os << self[i] << ", ";
    }
    return os;
}

std::vector <double> operator+(std::vector <double> self, std::vector <double> other){
    if (size(self) == size(other)){
        std::vector <double> res(size(self), 0);
        for (int i=0; i<size(self); i++){
            res[i] = self[i] + other[i];
        }
        return res;
    }
    else {
        std::cerr << "Размерности различны" << std::endl;
        return std::vector<double>();
    }
}

std::vector <double> operator*(double num, std::vector <double> self){
    std::vector <double> res(size(self), 0);
    for (int i=0; i<size(self); i++){
        res[i] = self[i] * num;
    }
    return res;
}

std::vector <double> operator-(std::vector <double> self, std::vector <double> other){
    return self + (-1 * other);
}


std::ostream& operator<<(std::ostream& os, std::vector <int> self){
    for (int i=0; i<self.size(); i++){
        os << self[i] << ", ";
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, std::vector <double> self){
    for (int i=0; i<self.size(); i++){
        os << self[i] << ", ";
    }
    return os;
}