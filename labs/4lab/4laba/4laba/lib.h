#pragma once
#include <functional>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>

struct point
{
    point();
    point(double x, double y);
    point(std::vector<double> x);
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
double f1(point p);
double f2(point p);
double f3(point p);
point operator-(point x, point y);
point operator+(point x, point y);
double operator*(point x, point y);
point operator*(double a, point y);
point operator*(double(&H)[2][2], point y);

double norm(point x);
