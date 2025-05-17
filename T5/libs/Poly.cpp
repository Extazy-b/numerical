// Подключаем необходимые библиотеки
#include <vector>
#include <sstream>
#include <cmath>
#include <ostream>
#include <iostream>
#include <string>

/**
 * @brief Класс для работы с многомерными полиномами
 */
class Poly {
private:
    int deg_;    // Степень полинома + 1
    int varc_;   // Количество переменных
    std::size_t len_;    // Общая длина массива коэффициентов
    std::vector<double> coefs_;  // Вектор коэффициентов

public:
    /**
     * @brief Конструктор класса
     * @param variablesCount - количество переменных
     * @param degree - степень полинома
     * @param defaultValue - значение по умолчанию для коэффициентов
     */
    Poly(int variablesCount = 0, int degree = 0, double defaultValue = 0.0) 
        : varc_(variablesCount), deg_(degree + 1) {
        len_ = static_cast<std::size_t>(std::pow(deg_, varc_));
        coefs_.resize(len_, defaultValue);
    }

    /**
     * @brief Получить степень полинома
     * @return Степень полинома
     */
    int getDeg() const {
        return deg_ - 1;
    }

    /**
     * @brief Получить количество переменных
     * @return Количество переменных
     */
    int getVarc() const {
        return varc_;
    }

    /**
     * @brief Получить длину массива коэффициентов
     * @return Длина массива коэффициентов
     */
    std::size_t getLen() const {
        return len_;
    }

    /**
     * @brief Преобразование линейного индекса в мультииндекс
     * @param index - линейный индекс
     * @return Вектор индексов по каждой переменной
     */
    std::vector<int> ind2multyind(int index) const {
        std::vector<int> res(varc_, 0);
        for (int k = 0; k < varc_; k++) {
            res[k] = static_cast<int>(index / std::pow(deg_, k)) % deg_;
        }
        return res;
    }

    /**
     * @brief Преобразование мультииндекса в линейный индекс
     * @param indices вектор индексов по каждой переменной
     * @return Линейный индекс
     */
    int multyind2ind(const std::vector<int>& indices) const {
        int idx = 0;
        for (int i = 0; i < varc_; ++i) {
            idx += indices[i] * static_cast<int>(std::pow(deg_, i));
        }
        return idx;
    }

    /**
     * @brief Ввод коэффициентов полинома из строки
     * @param str Строка с коэффициентами в формате:
     *        <количество вводов>
     *        [<координата 1> ... <координата n> <значение>]
     *        где все числа разделены пробелами
     */
    void setCoefs(const std::string& str) {
        std::stringstream ss(str);
        // пример ввода //6 3 0 1 1 1 -1 0 2 1 1 0 -2 0 1 3 0 0 -4
        double coef;
        int count = 0;
        
        ss >> count;
        for (size_t i = 0; i < count; i++)
        {
            std::vector <int> indices(varc_, 0);
            for (size_t ind = 0; ind < varc_; ind++)
            {
                ss >> indices[ind];
            }
            ss >> coef;
            coefs_[multyind2ind(indices)] = coef;
        }
    } 
  
    /**
     * @brief Увеличение количества переменных в полиноме
     * @param new_vars Новое количество переменных
     */
    Poly resize(int new_vars, int new_deg) const {
        if (new_vars <= varc_ && new_deg <= getDeg()) {
            return *this;
        }

        int actual_new_vars = std::max(varc_, new_vars);
        int actual_new_deg = std::max(getDeg(), new_deg);
        
        Poly res(actual_new_vars, actual_new_deg);
        
        // Копируем существующие коэффициенты
        for (std::size_t i = 0; i < len_; ++i) {
            if (coefs_[i] != 0.0) {
                std::vector<int> old_indices = ind2multyind(i);
                std::vector<int> new_indices(actual_new_vars, 0);
                
                // Копируем старые индексы в начало нового вектора
                std::copy(old_indices.begin(), old_indices.end(), new_indices.begin());
                
                // Вычисляем новый линейный индекс и копируем коэффициент
                int new_idx = res.multyind2ind(new_indices);
                res[new_idx] = coefs_[i];
            }
        }
        
        return res;
    }
    
    /**
     * @brief Доступ к коэффициентам по линейному индексу
     */
    double& operator[](int idx) {
        return coefs_[idx];
    }

    /**
     * @brief Константный доступ к коэффициентам по линейному индексу
     */
    const double& operator[](int idx) const {
        return coefs_[idx];
    }

    /**
     * @brief Доступ к коэффициентам по вектору индексов
     */
    double& operator()(const std::vector<int>& indices) {
        return coefs_[multyind2ind(indices)];
    }

    /**
     * @brief Константный доступ к коэффициентам по вектору индексов
     */
    const double& operator()(const std::vector<int>& indices) const {
        return coefs_[multyind2ind(indices)];
    }

    /**
     * @brief Доступ к коэффициентам по строке индексов
     */
    double& operator()(const std::string& indices) {
        std::vector<int> idx_vec;
        std::stringstream ss(indices);
        int num;
        while (ss >> num) {
            idx_vec.push_back(num);
        }
        
        return (*this)(idx_vec);
    }

    /**
     * @brief Константный доступ к коэффициентам по строке индексов
     */
    const double& operator()(const std::string& indices) const {
        std::vector<int> idx_vec;
        std::stringstream ss(indices);
        int num;
        while (ss >> num) {
            idx_vec.push_back(num);
        }
        
        return (*this)(idx_vec);
    }

    /**
     * @brief Сложение полиномов
     */
    friend Poly operator+(const Poly& lhs, const Poly& rhs) {
        // Используем логический оператор && вместо побитового &
        if ((lhs.getVarc() == rhs.getVarc()) && (lhs.getDeg() == rhs.getDeg())){
            Poly res(lhs.getVarc(), lhs.getDeg());
            for (std::size_t i = 0; i < res.getLen(); i++) {
                res[i] = lhs[i] + rhs[i];
            }
            return res;
        }
        else{
            int new_varc = std::max(lhs.getVarc(), rhs.getVarc());
            int new_deg = std::max(lhs.getDeg(), rhs.getDeg());
            
            Poly lhs_resized = lhs.resize(new_varc, new_deg);
            Poly rhs_resized = rhs.resize(new_varc, new_deg);

            Poly res(new_varc, new_deg);
            for (std::size_t i = 0; i < res.getLen(); i++) {
                res[i] = lhs_resized[i] + rhs_resized[i];
            }
            return res;
        }
    }
    
    /**
     * @brief Умножение полинома на скаляр
     */
    friend Poly operator*(double scalar, const Poly& poly) {
        Poly res(poly.getVarc(), poly.getDeg());
        for (std::size_t i = 0; i < poly.getLen(); i++) {
            res[i] = scalar * poly[i];
        }
        return res;
    }

    /**
     * @brief Умножение полиномов
     */
    friend Poly operator*(const Poly& lhs, const Poly& rhs) {
        if (lhs.getVarc() != rhs.getVarc()) {
            int new_varc = std::max(lhs.getVarc(), rhs.getVarc());
            return lhs.resize(new_varc, lhs.getDeg()) * rhs.resize(new_varc, rhs.getDeg());
        }
        
        Poly res(lhs.getVarc(), lhs.getDeg() + rhs.getDeg());
        
        for (std::size_t ind1 = 0; ind1 < lhs.getLen(); ind1++) {
            if (lhs[ind1] == 0.0) continue;
            
            std::vector<int> indices1 = lhs.ind2multyind(ind1);
            
            for (std::size_t ind2 = 0; ind2 < rhs.getLen(); ind2++) {
                if (rhs[ind2] == 0.0) continue;
                
                std::vector<int> indices2 = rhs.ind2multyind(ind2);
                std::vector<int> result_indices(lhs.getVarc());
                
                for (int k = 0; k < lhs.getVarc(); k++) {
                    result_indices[k] = indices1[k] + indices2[k];
                }
                
                res(result_indices) += lhs[ind1] * rhs[ind2];
            }
        }
        
        return res;
    }

    /**
     * @brief Вычитание полиномов
     */
    friend Poly operator-(const Poly& lhs, const Poly& rhs) {
        return lhs + (-1.0) * rhs;
    }

    /**
     * @brief Возведение полинома в степень
     */
    friend Poly pow(const Poly& poly, int power) {
        if (power < 0) {
            throw std::invalid_argument("Power must be non-negative");
        }
        if (power == 0) {
            Poly result(poly.getVarc(), 0);
            result(std::vector<int>(poly.getVarc(), 0)) = 1.0;
            return result;
        }
        if (power == 1) {
            return poly;
        }
        
        Poly result = poly;
        for (int i = 1; i < power; i++) {
            result = result * poly;
        }
        return result;
    }
    
    /**
     * @brief Вывод полинома в поток
     */    
    friend std::ostream& operator<<(std::ostream& os, const Poly& poly) {
        bool first = true;
        
        for (std::size_t ind = 0; ind < poly.getLen(); ind++) {
            if (poly[ind] != 0.0) {
                if (!first) {
                    if (poly[ind] > 0.0) {
                        os << " + ";
                    } else {
                        os << " - ";
                    }
                } else {
                    first = false;
                    if (poly[ind] < 0.0) {
                        os << "-";
                    }
                }
                
                double abs_coef = std::abs(poly[ind]);
                if (abs_coef != 1.0) {
                    os << abs_coef;
                }
                
                std::vector<int> indices = poly.ind2multyind(ind);
                bool has_vars = false;
                
                for (int subInd = 0; subInd < poly.getVarc(); subInd++) {
                    if (indices[subInd] > 0) {
                        if (has_vars) {
                            os << "*";
                        }
                        
                        os << "X" << (subInd+1);
                        if (indices[subInd] > 1) {
                            os << "^" << indices[subInd];
                        }
                        
                        has_vars = true;
                    }
                }
                
                if (!has_vars && abs_coef == 1.0) {
                    os << "1";
                }
            }
        }
        
        if (first) {
            os << "0";
        }
        
        return os;
    }

    /**
     * @brief Вычисление значения полинома при заданных значениях переменных
     * @param values - вектор значений переменных
     * @return Значение полинома
     */
    double evaluate(const std::vector<double>& values) const {
        double result = 0.0;
        
        for (std::size_t ind = 0; ind < len_; ind++) {
            if (coefs_[ind] == 0.0) continue;
            
            std::vector<int> indices = ind2multyind(ind);
            double term = coefs_[ind];
            
            for (int var = 0; var < varc_; var++) {
                term *= std::pow(values[var], indices[var]);
            }
            
            result += term;
        }
        
        return result;
    }


    Poly compose(std::vector <Poly> polyes){
        Poly res;

        for (size_t i = 0; i < len_; i++)
        {
            Poly tmp(0, 1, 1);
            std::vector <int> power = ind2multyind(i);
            for (size_t k = 0; k < varc_; k++)
            {
                tmp = tmp * pow(polyes[k], power[k]);
            }
            res = res + coefs_[i] * tmp;
        }
        return res;
    }
};


/**
 * @brief Вычисление значений полиномов в заданной точке
 * @param polys - вектор полиномов
 * @param values - значения переменных
 * @return Вектор значений полиномов
 */
std::vector<double> evaluate2list(std::vector<Poly> polys, const std::vector<double>& values){
    std::vector<double> results;
    results.reserve(polys.size());

    for (const auto& poly : polys) {
        results.push_back(poly.evaluate(values));
    }
    
    return results;
};

/**
 * @brief Вычисление производной полинома по заданной переменной
 * @param arg - полином
 * @param var - индекс переменной, по которой берется производная
 * @return Полином-производная
 */
Poly derivate(const Poly& arg, int var = 0) {
    Poly res(arg.getVarc(), arg.getDeg());
    
    for (std::size_t ind = 0; ind < arg.getLen(); ind++) {
        if (arg[ind] == 0.0) continue;
        
        std::vector<int> indices = arg.ind2multyind(ind);
        if (indices[var] == 0) continue;
        
        double coef = arg[ind] * indices[var];
        indices[var] -= 1;
        
        res(indices) = coef;
    }
    
    return res;
}

/**
 * @brief Вычисление градиента полинома
 * @param arg - полином
 * @param var - индекс переменной (если -1, то по всем переменным)
 * @return Вектор полиномов-производных
 */
std::vector<Poly> Nabla(const Poly& arg) {
    std::vector<Poly> res;
    res.reserve(arg.getVarc());
    
    for (int i = 0; i < arg.getVarc(); i++) {
        res.push_back(derivate(arg, i));
    }
    
    return res;
}

/**
 * @brief Преобразование градиента в строковое представление
 * @param arg - вектор полиномов (градиент)
 * @return Строковое представление градиента
 */
std::string Nabla(const std::vector<Poly>& arg) {
    std::stringstream ss;
    ss << "Gradient: [";
    
    for (size_t i = 0; i < arg.size(); i++) {
        if (i > 0) {
            ss << ", ";
        }
        ss << arg[i];
    }
    
    ss << "]";
    return ss.str();
};


