#pragma once
#include <functional>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>
double eps_f = 1e-3;
int count_f;


struct point
{
    double x;//координата х
    double y;//координата y
    double f;//значение функции в точке

};
typedef std::function<double(point)> OneDimensionFunction;

struct test
{
    OneDimensionFunction func;
    std::string fname;
    double eps;

};

double f1(point p) {
    double x = p.x;
    double y = p.y;
    return (100 * (y - x) * (y - x) + (1 - x) * (1 - x));
}
double f2(point p) {
    double x = p.x;
    double y = p.y;
    return 100 * (y - x * x) * (y - x * x) + (1 - x) * (1 - x);
}
double f3(point p) {
    double x = p.x;
    double y = p.y;
    return -1 * (3 / (1 + pow(x - 2, 2) + pow((y - 2) / 2, 2)) + 2 / (1 + pow((x - 2) / 3, 2) + pow(y - 3, 2)));
}

point operator-(point x, point y) {
    x.x = x.x - y.x;
    x.y = x.y - y.y;
    return x;
}
point operator+(point x, point y) {
    x.x = x.x + y.x;
    x.y = x.y + y.y;
    return x;
}

double operator*(point x, point y) {
    return x.x * y.x + x.y * y.y;
}

point operator*(double a, point y) {
    y.x *= a;
    y.y *= a;
    return y;
}

point operator*(double (&H)[2][2], point y) {
    point m;
    m.x = H[0][0] * y.x + H[0][1] * y.y;
    m.y = H[1][0] * y.x + H[1][1] * y.y;
    return m;
}

double norm (point x) {
    return sqrt(x.x * x.x + x.y * x.y);
}

