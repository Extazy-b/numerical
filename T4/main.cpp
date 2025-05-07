#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>
#include <cstdlib> 

using namespace std;

void list2bin(const vector<vector<vector<double>>>& U, double dt) {
    // Открываем бинарный файл для данных
    ofstream dataFile("output.bin", ios::binary);
    if (!dataFile) {
        cerr << "Error opening data file for writing." << endl;
        return;
    }
    
    // Открываем индексный файл
    ofstream indexFile("output.idx");
    if (!indexFile) {
        cerr << "Error opening index file for writing." << endl;
        dataFile.close();
        return;
    }
    
    // Записываем размеры массива в индексный файл
    indexFile << U.size() << " " << U[0].size() << " " << U[0][0].size() << endl;
    
    // Размер одного временного шага в байтах
    size_t timeStepSize = U[0].size() * U[0][0].size() * sizeof(double);
    
    // Для каждого временного шага
    for (size_t n = 0; n < U.size(); ++n) {
        // Записываем позицию в файле для этого временного шага
        streampos position = dataFile.tellp();
        indexFile << n << " " << n*dt << " " << position << endl;
        
        // Записываем данные временного шага в бинарном формате
        for (size_t i = 0; i < U[n].size(); ++i) {
            dataFile.write(reinterpret_cast<const char*>(U[n][i].data()), 
                          U[n][i].size() * sizeof(double));
        }
        
        // Выводим прогресс в консоль
        if (n % (U.size() / 10) == 0 || n == U.size() - 1) {
            cout << "Progress: " << n + 1 << "/" << U.size() << " time steps processed" << endl;
        }
    }
    
    dataFile.close();
    indexFile.close();
    cout << "Data saved to output.bin with index file for fast access" << endl;
}


void list2text(const vector<vector<vector<double>>>& U, double dt) {
    ofstream outFile("output.txt");
    if (!outFile) {
        cerr << "Error opening file for writing." << endl;
        return;
    }

    // Определяем ширину поля для каждого элемента (например, 15 символов)
    const int fieldWidth = 15;
    // Устанавливаем точность для чисел с плавающей точкой
    const int precision = 6;
    
    // Настраиваем форматирование вывода
    outFile << fixed << setprecision(precision);
    
    // Записываем размеры массива в первую строку с фиксированной длиной
    outFile << setw(fieldWidth) << U.size() << setw(fieldWidth) << U[0].size() 
            << setw(fieldWidth) << U[0][0].size() << endl;
    
    // Для каждого временного шага
    for (size_t n = 0; n < U.size(); ++n) {
        // Записываем информацию о временном шаге
        outFile << setw(fieldWidth) << n << setw(fieldWidth) << n*dt;
        
        // Заполняем оставшуюся часть строки пробелами до фиксированной длины
        for (size_t k = 0; k < U[0][0].size(); ++k) {
            outFile << setw(fieldWidth) << "";
        }
        outFile << endl;
        
        // Записываем данные для каждой строки матрицы
        for (size_t i = 0; i < U[n].size(); ++i) {
            for (size_t j = 0; j < U[n][i].size(); ++j) {
                outFile << setw(fieldWidth) << U[n][i][j];
            }
            outFile << endl;
        }
        
        // Выводим прогресс в консоль
        if (n % (U.size() / 10) == 0 || n == U.size() - 1) {
            cout << "Progress: " << n + 1 << "/" << U.size() << " time steps processed" << endl;
        }
    }

    outFile.close();
    cout << "Data saved to output.txt with fixed width formatting" << endl;
}


void printUsage() {
    cout << "Использование: ./program [параметры]" << endl;
    cout << "Параметры:" << endl;
    cout << "  -I <число>     : Двоичный прядок кооличества точек по оси X (по умолчанию 4)" << endl;
    cout << "  -J <число>     : Двоичный прядок кооличества точек по оси Y (по умолчанию 4)" << endl;
    cout << "  -N <число>     : Двоичный прядок кооличества точек временных шагов (по умолчанию 10)" << endl;
    cout << "  -L <число>     : Длина области по X (по умолчанию 1.0)" << endl;
    cout << "  -M <число>     : Длина области по Y (по умолчанию 1.0)" << endl;
    cout << "  -T <число>     : Конечное время (по умолчанию 1.0)" << endl;
    cout << "  -o <имя_файла> : Имя выходного файла (по умолчанию output.txt)" << endl;
    cout << "  -h             : Показать эту справку" << endl;
}

int main(int argc, char* argv[]) {
    // Значения по умолчанию
    double L = 1.0;
    double M = 1.0;
    double T = 1.0;
    int I = (1 << 4); // 16
    int J = (1 << 4); // 16
    int N = (1 << 12); // 1024
    double delta = (1 << 31);

    string outputFile = "output.txt";

    // Обработка аргументов командной строки
    for (int i = 1; i < argc; i++) {
        string arg = argv[i];
        
        if (arg == "-h" || arg == "--help") {
            printUsage();
            return 0;
        } else if (arg == "-L" && i + 1 < argc) {
            L = atof(argv[++i]);
        } else if (arg == "-M" && i + 1 < argc) {
            M = atof(argv[++i]);
        } else if (arg == "-T" && i + 1 < argc) {
            T = atof(argv[++i]);
        } else if (arg == "-I" && i + 1 < argc) {
            I = L * (1 << atoi(argv[++i]));
        } else if (arg == "-J" && i + 1 < argc) {
            J = M * (1 << atoi(argv[++i]));
        } else if (arg == "-N" && i + 1 < argc) {
            N = T * (1 << atoi(argv[++i]));
        } else if (arg == "-o" && i + 1 < argc) {
            outputFile = argv[++i];
        } else {
            cout << "Неизвестный аргумент: " << arg << endl;
            printUsage();
            return 1;
        }
    }

    // Проверка корректности параметров
    if (I <= 0 || J <= 0 || N <= 0 || L <= 0 || M <= 0 || T <= 0) {
        cout << "Ошибка: все числовые параметры должны быть положительными" << endl;
        return 1;
    }

    double dx = L / I;
    double dy = M / J;
    double dt = T / N;
    bool const_flag = 1;
    int const_n = N;


    cout << "Параметры расчета:" << endl;
    cout << "Размерность сетки: " << I << " x " << J << " x " << N << endl;
    cout << "Область: [0, " << L << "] x [0, " << M << "] x [0, " << T << "]" << endl;
    cout << "Шаги: dx=" << dx << ", dy=" << dy << ", dt=" << dt << endl;
    cout << "Порядок схемы: " << dt/(dx*dx + dy*dy) << endl;
    cout << "Прогрешность стабильности: " << 1/delta << endl;
    cout << "Выходной файл: " << outputFile << endl;
    cout << "Требуемая память (ГБ): " << sizeof(double) * I * J * N / (1024.0 * 1024.0 * 1024.0) << endl;

    // Проверка условия устойчивости CFL
    double cfl_limit = 1.0 / 1.0/(dx*dx) + 3.0/(dy*dy);
    if (dt > cfl_limit) {
        cout << "Предупреждение: нарушено условие устойчивости CFL!" << endl;
        cout << "Текущий шаг по времени: " << dt << endl;
        cout << "Максимально допустимый шаг: " << cfl_limit << endl;
        cout << "Рекомендуется увеличить N до " << ceil(T/cfl_limit) << " или больше" << endl;
    }

    cout << "Нажмите Enter чтобы начать расчёт" << endl;
    
    cin.get();

    vector <vector <vector <double>>> U(N + 1, vector <vector <double>> (I + 1, vector <double> (J + 1, 0)));


    for (int i = 0; i <= I; i++) {
        for (int j = 0; j <= J; j++) {
            U[0][i][j] = 2*i*dx + j*dy;
        }
    }

    for (int n = 0; n < N; n++) {
        if (const_flag) {
            const_flag = 0;

            for (int i = 0; i <= I; i++) {
                U[n][i][0] = 2;
                U[n][i][J] = 2;
            }
            for (int j = 0; j <= J; j++) {
                U[n][0][j] = 2;
                U[n][I][j] = 2;
            }

            for (int i = 1; i < I; i++) {
                for (int j = 1; j < J; j++) {
                    double Uxx = (U[n][i + 1][j] - 2 * U[n][i][j] + U[n][i - 1][j]) / (dx * dx);
                    double Uyy = (U[n][i][j + 1] - 2 * U[n][i][j] + U[n][i][j - 1]) / (dy * dy);
                    U[n + 1][i][j] = U[n][i][j] + dt * (Uxx + 3 * Uyy);
                    if (U[n+1][i][j] - U[n][i][j] > 1/delta) const_flag = 1;
                }
            }

            if (n % (N / 10) == 0) {
                cout << n << endl;
            }
        }
        else {
            cout << "Ситстема уравновесилась (шаг: " << n << ", время:" << n*dt << ")" << endl;
            const_n = n;
            break;
        }
        
    }

    if (const_flag) {
        for (int i = 0; i <= I; i++) {
            U[N][i][0] = 2;
            U[N][i][J] = 2;
        }
        for (int j = 0; j <= J; j++) {
            U[N][0][j] = 2;
            U[N][I][j] = 2;
        }  
    }
    else{
        U.resize(const_n + 1);
        U[const_n] = U[const_n - 1];
    }

    cout << "done" << endl;

//TODO Динамическая запись неиспользуемой части матрицы

//    for (size_t i = 0; i < U[const_n].size(); ++i) {
//        for (size_t j = 0; j < U[const_n][i].size(); ++j) {
//            cout << U[const_n][i][j] - U[const_n - 1][i][j] << " ";
//        }
//        cout << "\n";
//    }
    list2text(U, dt);
    cout << "done" << endl;
    list2bin(U, dt);
    cout << "done" << endl;

    return 0;
}
