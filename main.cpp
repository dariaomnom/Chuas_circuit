#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>

#define NUM_POINTS 10000

const double C1 = 1.0 / 9.0;
const double C2 = 1.0;
const double L = 1.0 / 7.0;
const double G = 0.7;
const double Ga = -0.8;
const double Gb = -0.5;

const double m0 = Ga / G;
const double m1 = Gb / G;
double a = C2 / C1;
double b = C2 / (G * G * L);

const double h = 0.01;

class RenderWindow;

double rand(double start, double stop) {
    double res = rand() / (double)RAND_MAX * (stop - start) + start;
    return res;
}

void make_start_points(double* points, const sf::Color* colors_static[], double num, double start, double stop) {
    for (int i = 0; i < num; i++) {
        for (int j = 0; j < 3; j++) {
            points[i * 3 + j] = rand(start, stop);
        }
        if (points[i * 3 + 0] < 0 && points[i * 3 + 1] <= 0) {
            colors_static[i] = &sf::Color::Red;
        } else if (points[i * 3 + 0] < 0 && points[i * 3 + 1] > 0) {
            colors_static[i] = &sf::Color::Green;
        } else if (points[i * 3 + 0] >= 0 && points[i * 3 + 1] <= 0) {
            colors_static[i] = &sf::Color::Yellow;
        } else if (points[i * 3 + 0] >= 0 && points[i * 3 + 1] > 0) {
            colors_static[i] = &sf::Color::Blue;
        }
    }
}

double H(double x) {
    return m1 * x + ((m0 - m1) / 2) * (fabs(x + 1) - fabs(x - 1));
}

void f_Chua(double X, double Y, double Z, double* dX) {
    dX[0] = a * (Y - X - H(X));
    dX[1] = X - Y + Z;
    dX[2] = -b * Y;
}

void RK4(double* X, double* X_tmp, double* k1, double* k2, double* k3, double* k4) {
    for (int j = 0; j < NUM_POINTS; j++) {
        f_Chua(X[j*3 + 0], X[j*3 + 1], X[j*3 + 2], k1);
        for (int i = 0; i < 3; i++) {
            X_tmp[i] = X[j*3 + i] + h * 0.5 * k1[i];
        }

        f_Chua(X_tmp[0], X_tmp[1], X_tmp[2], k2);
        for (int i = 0; i < 3; i++) {
            X_tmp[i] = X[j*3 + i] + h * 0.5 * k2[i];
        }

        f_Chua(X_tmp[0], X_tmp[1], X_tmp[2], k3);
        for (int i = 0; i < 3; i++) {
            X_tmp[i] = X[j*3 + i] + h * 1 * k3[i];
        }

        f_Chua(X_tmp[0], X_tmp[1], X_tmp[2], k4);
        for (int i = 0; i < 3; i++) {
            X[j*3 + i] += h * (1.0 / 6.0 * k1[i] + 1.0 / 3.0 * k2[i] + 1.0 / 3.0 * k3[i] + 1.0 / 6.0 * k4[i]);
        }
    }
}


struct Vector3 {
    float x, y, z;

    Vector3(float x_, float y_, float z_) : x(x_), y(y_), z(z_) {}

    // перегрузка оператора умножения для умножения вектора на матрицу
    Vector3 operator*(const float matrix[3][3]) const {
        return Vector3(
                x * matrix[0][0] + y * matrix[0][1] + z * matrix[0][2],
                x * matrix[1][0] + y * matrix[1][1] + z * matrix[1][2],
                x * matrix[2][0] + y * matrix[2][1] + z * matrix[2][2]
        );
    }
};



int main() {
    auto k1 = new double[3];
    auto k2 = new double[3];
    auto k3 = new double[3];
    auto k4 = new double[3];

    auto points = new double[NUM_POINTS * 3];
    auto X_tmp = new double[3];
    const sf::Color* colors_static[NUM_POINTS]; // массив для выбора цвета в sfml

    make_start_points(points, colors_static, NUM_POINTS, -2.0, 2.0);

    // параметры окна анимации
    const unsigned int windowWidth = 800;
    const unsigned int windowHeight = 800;
    sf::RenderWindow window(sf::VideoMode(windowWidth, windowHeight), "Chua's attractor");
    sf::View view(sf::FloatRect(-5.f, -5.f, 10.f, 10.f));
    view.setViewport(sf::FloatRect(0.f, 0.f, 1.f, 1.f));
    window.setView(view);

    // перенос начальных точек в массив векторов координат
    std::vector<sf::Vector2f> coordinates;
    coordinates.reserve(NUM_POINTS);
    for (int i = 0; i < NUM_POINTS; ++i) {
        coordinates.emplace_back(points[i*3 + 0], points[i*3 + 1]);
    }

    // 10 миллисекунд = 100 кадров в секунду
    sf::Clock clock;
    const float desiredFrameTime = 10.0f;

    // углы вращения и матрица поворота координат
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
    bool is_rotated = false;

    while (window.isOpen()) {
        // считывание нажатий клавиш X, Y, Z, A, B, Escape (не забудьте переключить раскладку)
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
            else if (event.type == sf::Event::KeyPressed) {
                if (event.key.code == sf::Keyboard::Escape)
                    window.close();
                else if (event.key.code == sf::Keyboard::X) {
                    angleX += 0.1;
                    is_rotated = true;
                    std::cout << "X-rotation" << std::endl;
                }
                else if (event.key.code == sf::Keyboard::Y) {
                    angleZ += 0.1;
                    is_rotated = true;
                    std::cout << "Y-rotation" << std::endl;
                }
                else if (event.key.code == sf::Keyboard::Z) {
                    angleY += 0.1;
                    is_rotated = true;
                    std::cout << "Z-rotation" << std::endl;
                }
                else if (event.key.code == sf::Keyboard::A) {
                    a *= 1.2;
                    std::cout << "Alpha changed (a * 1.2)" << std::endl;
                    std::cout << "Alpha = " << a << std::endl;
                }
                else if (event.key.code == sf::Keyboard::B) {
                    b *= 1.2;
                    std::cout << "Beta changed (b * 1.2)" << std::endl;
                    std::cout << "Beta = " << b << std::endl;
                }
            }

        }

        if (is_rotated) {
            // обновление матрицы поворота, если было вращение по осям координат
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

        // проверка времени, прошедшего с последнего кадра
        sf::Time elapsed = clock.restart();
        if (elapsed.asMilliseconds() < desiredFrameTime)
            sf::sleep(sf::milliseconds(desiredFrameTime - elapsed.asMilliseconds()));

        window.clear();

        // обработка всех точек
        RK4(points, X_tmp, k1, k2, k3, k4);

        // вывод точек в цикле
        for (int i = 0; i < NUM_POINTS; i++) {
            Vector3 v(points[i*3 + 0], points[i*3 + 1], points[i*3 + 2]);
            Vector3 v_rotated = v * Rotation;
            coordinates[i].x = v_rotated.x;
            coordinates[i].y = v_rotated.y;

            sf::CircleShape point_new(0.015f);
            point_new.setPosition(coordinates[i]);
            point_new.setFillColor(*colors_static[i]);

            window.draw(point_new);
        }
        window.display();
    }
}
