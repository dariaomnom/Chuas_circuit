#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
//#include <math.h>
#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>

//const double C1 = 1.0 / 9.0;
//const double C2 = 1.0;
//const double L = 1.0 / 7.0;
//const double G = 0.7;
//const double Ga = -0.8;
//const double Gb = -0.5;
//
//const double m0 = Ga / G;
//const double m1 = Gb / G;
//const double a = C2 / C1;
//const double b = C2 / (G * G * L);
//
//const double T = 10000.0;
//const double h = 0.01;

double C1 = 1.0 / 9.0;
double C2 = 1.0;
double L = 1.0 / 7.0;
double G = 0.7;
double Ga = -0.8;
double Gb = -0.5;

double m0 = Ga / G;
double m1 = Gb / G;
double a = C2 / C1;
double b = C2 / (G * G * L);

double T = 10000.0;
double h = 0.01;

class RenderWindow;

double rand(double start, double stop) {
    double res = rand() / (double)RAND_MAX * (stop - start) + start;
    return res;
}

void make_start_points(double** points, int* colors, double num, double start, double stop) {
    for (int i = 0; i < num; i++) {
        for (int j = 0; j < 3; j++) {
            points[i][j] = rand(start, stop);
        }
        if (points[i][0] < 0 && points[i][1] <= 0) {
            colors[i] = 1;
        } else if (points[i][0] < 0 && points[i][1] > 0) {
            colors[i] = 2;
        } else if (points[i][0] >= 0 && points[i][1] <= 0) {
            colors[i] = 3;
        } else if (points[i][0] >= 0 && points[i][1] > 0) {
            colors[i] = 4;
        }
    }
}



double H(double x) {
    return m1 * x + ((m0 - m1) / 2) * (fabs(x + 1) - fabs(x - 1));
}

void f_Chua(double* X, double* dX) {
    dX[0] = a * (X[1] - X[0] - H(X[0]));
    dX[1] = X[0] - X[1] + X[2];
    dX[2] = -b * X[1];
}

void RK4(double* X, double* X_tmp, double* k1, double* k2, double* k3, double* k4) {
    f_Chua(X, k1);
    for (int i = 0; i < 3; i++) {
        X_tmp[i] = X[i] + h * 0.5 * k1[i];
    }

    f_Chua(X_tmp, k2);
    for (int i = 0; i < 3; i++) {
        X_tmp[i] = X[i] + h * 0.5 * k2[i];
    }

    f_Chua(X_tmp, k3);
    for (int i = 0; i < 3; i++) {
        X_tmp[i] = X[i] + h * 1 * k3[i];
    }

    f_Chua(X_tmp, k4);
    for (int i = 0; i < 3; i++) {
        X[i] += h * (1.0 / 6.0 * k1[i] + 1.0 / 3.0 * k2[i] + 1.0 / 3.0 * k3[i] + 1.0 / 6.0 * k4[i]);
    }
}

// Структура для представления вектора
struct Vector3 {
    float x, y, z;

    Vector3(float x_, float y_, float z_) : x(x_), y(y_), z(z_) {}

    // Перегрузка оператора умножения для умножения вектора на матрицу
    Vector3 operator*(const float matrix[3][3]) const {
        return Vector3(
                x * matrix[0][0] + y * matrix[0][1] + z * matrix[0][2],
                x * matrix[1][0] + y * matrix[1][1] + z * matrix[1][2],
                x * matrix[2][0] + y * matrix[2][1] + z * matrix[2][2]
        );
    }
};



int main() {
    int num_points = 10000;

    auto k = new double**[num_points];
    for (int i = 0; i < num_points; i++) {
        k[i] = new double*[4];
        for (int j = 0; j < 4; j++) {
            k[i][j] = new double[3];
        }
    }

    auto start_points = new double*[num_points];
    for (int i = 0; i < num_points; i++) {
        start_points[i] = new double[3];
    }

    int* original_points = new int[num_points];

    make_start_points(start_points, original_points, num_points, -5.0, 5.0);
//    make_start_points(start_points, original_points, 800, -1.0, 1.0);



//          Для создания графика в gnuplot
// Для стандарта C++11 и новее:
//    double X[3] {-1.0, 1.0, 1.0};
//    double X_tmp[3] {};
//    double k1[3] {}, k2[3] {}, k3[3], k4[3] {};

// Если предыдущий вариант вызывает ошибки:
//    double X[3];
    double X_tmp[3];
//    X[0] = -1.0;
//    X[1] = 1.0;
//    X[2] = 1.0;
//    double k1[3];
//    double k2[3];
//    double k3[3];
//    double k4[3];


//          Для создания графика в gnuplot
//    std::ofstream outputFile("../data3d.txt");
//    for (double t = 0; t < T; t += h) {
//        outputFile << std::fixed << std::setprecision(16) << X[0] << " " << X[1]  << " " << X[2] << std::endl;
//        RK4(X, X_tmp, k1, k2, k3, k4);
//    }


//          вывод текста
        sf::Clock clock;
        sf::Font font;
        if (!font.loadFromFile("../Montserrat-Regular.ttf")) {
            return EXIT_FAILURE;
        }
        sf::Text text;
        text.setFont(font);
        text.setCharacterSize(7);
        text.setFillColor(sf::Color::Red);
        text.setString("A");
        text.setPosition(0, 0);
        text.setOutlineThickness(0.5);
//        text.setOutlineColor(sf::Color::Red);
        sf::Text textB;
        textB.setFont(font);
        textB.setCharacterSize(7);
        textB.setFillColor(sf::Color::Red);
        textB.setString("B");
        textB.setPosition(0, 0);
//        textB.setOutlineThickness(0);
        bool showText = false;
        bool showTextB = false;
//          вывод текста


        const unsigned int windowWidth = 800;
        const unsigned int windowHeight = 800;
        sf::RenderWindow window(sf::VideoMode(windowWidth, windowHeight), "Attractor");
        sf::View view(sf::FloatRect(-5.f, -5.f, 10.f, 10.f));
//        sf::View view(sf::FloatRect(-20.f, -20.f, 40.f, 40.f));
        view.setViewport(sf::FloatRect(0.f, 0.f, 1.f, 1.f));
        window.setView(view);

        std::vector<sf::Vector2f> coordinates;
        for (int i = 0; i < num_points; ++i) {
            coordinates.push_back(sf::Vector2f(start_points[i][0], start_points[i][1]));
        }

//        sf::Clock clock;
        const float desiredFrameTime = 10.0f; // 100 мс = 10 кадров в секунду


        float angleX = 0.0f;
        float angleY = 0.0f;
        float angleZ = 0.0f;
        float Rotation[3][3] = {
            {std::cos(angleY)*std::cos(angleZ),
                    -std::sin(angleZ)*std::cos(angleY),
                    std::sin(angleY)},
            {std::sin(angleX)*std::sin(angleY)*std::cos(angleZ)+std::sin(angleZ)*std::cos(angleX),
                    -std::sin(angleX)*std::sin(angleY)*std::sin(angleZ)+std::cos(angleX)*std::cos(angleZ),
                    -std::sin(angleX)*std::cos(angleY)},
            {std::sin(angleX)*std::sin(angleZ)-std::sin(angleY)*std::cos(angleX)*std::cos(angleZ),
                    std::sin(angleX)*std::cos(angleZ)+std::sin(angleY)*std::sin(angleZ)*std::cos(angleX),
                    std::cos(angleX)*std::cos(angleY)}
        };
        bool flag = false;


        double t = 0;
        while (window.isOpen())
        {
            sf::Event event;
            while (window.pollEvent(event))
            {
                if (event.type == sf::Event::Closed)
                    window.close();
                else if (event.type == sf::Event::KeyPressed)
                {
                    if (event.key.code == sf::Keyboard::Escape)
                        window.close();
                    else if (event.key.code == sf::Keyboard::X) {
                        angleX += 0.1;
//                        std::cout << "X-rotation" << std::endl;
                        flag = true;
                    }
                    else if (event.key.code == sf::Keyboard::Y) {
                        angleZ += 0.1;
//                        std::cout << "Y-rotation" << std::endl;
                        flag = true;
                    }
                    else if (event.key.code == sf::Keyboard::Z) {
                        angleY += 0.1;
//                        std::cout << "Z-rotation" << std::endl;
                        flag = true;
                    }
                    else if (event.key.code == sf::Keyboard::A) {
//                        showText = true;
//                        clock.restart();
                        a *= 1.2;
                    }
                    else if (event.key.code == sf::Keyboard::B) {
//                        showTextB = true;
//                        clock.restart();
                        b *= 1.2;
                    }
                }

            }
            if (flag) {
                Rotation[0][0] = std::cos(angleY) * std::cos(angleZ);
                Rotation[0][1] = -(std::sin(angleZ) * std::cos(angleY));
                Rotation[0][2] = std::sin(angleY);
                Rotation[1][0] = std::sin(angleX) * std::sin(angleY) * std::cos(angleZ) +
                                 std::sin(angleZ) * std::cos(angleX);
                Rotation[1][1] = -(std::sin(angleX) * std::sin(angleY) * std::sin(angleZ)) +
                                 std::cos(angleX) * std::cos(angleZ);
                Rotation[1][2] = -(std::sin(angleX) * std::cos(angleY));
                Rotation[2][0] = std::sin(angleX) * std::sin(angleZ) -
                                 (std::sin(angleY) * std::cos(angleX) * std::cos(angleZ));
                Rotation[2][1] = std::sin(angleX) * std::cos(angleZ) +
                                 std::sin(angleY) * std::sin(angleZ) * std::cos(angleX);
                Rotation[2][2] = std::cos(angleX) * std::cos(angleY);
            }

            sf::Time elapsed = clock.restart();
            if (elapsed.asMilliseconds() < desiredFrameTime)
                sf::sleep(sf::milliseconds(desiredFrameTime - elapsed.asMilliseconds()));

            window.clear();

            if (t < T) {
                for (int i = 0; i < num_points; i++) {
                    Vector3 v(start_points[i][0], start_points[i][1], start_points[i][2]);
                    Vector3 v_rotated = v * Rotation;
                    coordinates[i].x = v_rotated.x;
                    coordinates[i].y = v_rotated.y;

//                    coordinates[i].x = start_points[i][0];
//                    coordinates[i].y = start_points[i][1];

                    sf::CircleShape point_new(0.020f);

                    point_new.setPosition(coordinates[i]);
                    if (original_points[i] == 1) {
                        point_new.setFillColor(sf::Color::Red);
                    } else if (original_points[i] == 2) {
                        point_new.setFillColor(sf::Color::Green);
                    } else if (original_points[i] == 3) {
                        point_new.setFillColor(sf::Color::Yellow);
                    } else if (original_points[i] == 4) {
                        point_new.setFillColor(sf::Color::Blue);
                    }
                    else {
                        point_new.setFillColor(sf::Color::Red);
                    }

                    window.draw(point_new);
                    RK4(start_points[i], X_tmp, k[i][0], k[i][1], k[i][2], k[i][3]);
                }
                t += h;
            }
//            if (showText) {
//                window.draw(text);
//            }

            window.display();

//            if (showText && clock.getElapsedTime().asSeconds() > 3.0) {
//                showText = false;
//            }

        }

}
