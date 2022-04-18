#pragma once
#include <iostream>
#include<cmath>
#include"lib.cpp"
using namespace std;

class method_variable_metric
{
	double delta = 1E-8;
	double diff_x(point point1)//маасив аргументов/по какой переменной/в какую функцию подставляем/шаг
	{
		double h = 10E-7;
		double d0 = point1.x;
		point1.x += h;
		double f_right = function(point1.x, point1.y);
		point1.x = d0 - h;
		double f_left = function(point1.x, point1.y);
		point1.x = d0;
		return (f_right - f_left) / (2 * h);
	}
	double diff_y(point point1)//маасив аргументов/по какой переменной/в какую функцию подставляем/шаг
	{
		double h = 10E-7;
		double d0 = point1.y;
		point1.y += h;
		double f_right = function(point1.x, point1.y);
		point1.y = d0 - h;
		double f_left = function(point1.x, point1.y);
		point1.y = d0;
		return (f_right - f_left) / (2 * h);
	}
	point find_grad(point point1)
	{
		point grad;
		grad.x = diff_x(point1);
		grad.y = diff_y(point1);
		return grad;
	}

	void IntervalSearch(double lymbda, double* a, double* b, point point1, point p0)
	{
		double h = 0;
		double lymbda_previous = lymbda;
		//шаг 1. определяем направление поиска.
		if (function(point1.x + p0.x * lymbda, point1.y + p0.y * lymbda) > function(point1.x + p0.x * (lymbda + delta), point1.y + p0.y * (lymbda + delta)))
		{
			lymbda += delta;
			h = delta;
		}
		else
		{
			lymbda -= delta;
			h = -delta;
		}
		h *= 2;
		lymbda_previous = lymbda;
		lymbda = lymbda + h;
		while (function(point1.x + p0.x * lymbda_previous, point1.y + p0.y * lymbda_previous) > function(point1.x + p0.x * lymbda, point1.y + p0.y * lymbda))
		{
			h *= 2;

			lymbda_previous = lymbda;
			lymbda = lymbda + h;
		}

		*a = lymbda_previous;
		*b = lymbda;
	}
	double GoldenSectionMethod(double a, double b, point point1, point pk)
	{
		int i;
		int fcount = 0;
		const double A_COEFF = ((3 - sqrt(5.0)) / 2);
		const double B_COEFF = ((sqrt(5.0) - 1) / 2);
		double x1 = a + A_COEFF * (b - a);
		double x2 = a + B_COEFF * (b - a);
		double  f1 = function(point1.x + pk.x * x1, point1.y + pk.y * x1);
		double  f2 = function(point1.x + pk.x * x2, point1.y + pk.y * x2);
		double a00, b00;//значение на предыдущей итерации
		fcount += 2;

		for (i = 0; abs(b - a) > 1E-8; i++)
		{
			a00 = a;
			b00 = b;
			if (f1 > f2)
			{
				a = x1;
				x1 = x2;
				f1 = f2;
				x2 = a + B_COEFF * (b - a);
				f2 = function(point1.x + pk.x * x2, point1.y + pk.y * x2);
			}

			else
			{
				b = x2;
				x2 = x1;
				f2 = f1;
				x1 = a + A_COEFF * (b - a);
				f1 = function(point1.x + pk.x * x1, point1.y + pk.y * x1);
			}
		}
		return (a + b) / 2;
	}

	double search_lambd(double lambd, point point1, point pk)
	{
		double a = 0, b = 0;
		IntervalSearch(lambd, &a, &b, point1, pk);
		lambd = GoldenSectionMethod(a, b, point1, pk);
		return lambd;
	}
	void Gessian(double(*H0)[2][2], point _sk, point _yk)
	{
		double sk[2] = { _sk.x, _sk.y };
		double yk[2] = { _yk.x, _yk.y };
		double multiplier[2] = { 0 }, divider;
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
				multiplier[i] += (*H0)[i][j] * yk[j];
			multiplier[i] = sk[i] - multiplier[i];
		}
		divider = multiplier[0] * yk[0] + multiplier[1] * yk[1];
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				(*H0)[i][j] += multiplier[i] * multiplier[j] / divider;

		/*int I[2][2] = {{ 1,0 }, { 0,1 }};
		double yk_sk[2][2] = { {yk.x * sk.x,yk.x * sk.y},{yk.y * sk.x,yk.y * sk.y} };
		double sk_sk[2][2]= { {sk.x * sk.x,sk.x * sk.y},{sk.y * sk.x,sk.y * sk.y} };
		double k=1/(yk.x*sk.x+yk.y*sk.y);
		double sk_yk[2][2] = {{sk.x*yk.x,sk.x*yk.y},{sk.y * yk.x,sk.y * yk.y}};
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				(*H0)[i][j] = (I[i][j] - k * sk_yk[i][j]) * (I[i][j] - k * yk_sk[i][j])* (*H0)[i][j] + sk_sk[i][j];*/
	}
public:
	method_variable_metric(point point1)
	{
		FILE* File;
		fopen_s(&File, "output.txt", "w");
		fprintf(File, "%-5i%|-20s%|-20s%|-20s%|-20s%|-20s%|\n", "i", "(xi, yi)", "f(xi, yi)", "(s1, s2)", "Lambd", " |xi-xi-1|", "angle", "H");

		double check;
		double lambd = 0;
		point grad1 = { 0,0 }, grad2, pk, yk, sk, point2;
		double H0[2][2] = { 0 };
		H0[0][0] = H0[1][1] = 1;//единичная матрица
		grad1 = find_grad(point1);
		cout << function(point1.x, point1.y) << endl;
		check = sqrt(pow(grad1.x, 2) + pow(grad1.y, 2));
		for (int iter = 1; check > eps_f; iter++)
		{
			
			pk.x = -(H0[0][0] * grad1.x + H0[0][1] * grad1.y);
			pk.y = -(H0[1][0] * grad1.x + H0[1][1] * grad1.y);
			lambd = search_lambd(lambd, point1, pk);
			point2.x = point1.x + lambd * pk.x;
			point2.y = point1.y + lambd * pk.y;
			grad2 = find_grad(point2);
			sk.x = point2.x - point1.x;
			sk.y = point2.y - point1.y;
			yk.x = grad2.x - grad1.x;
			yk.y = grad2.y - grad1.y;
			cout << abs(point2.x - point1.x) << " " << abs(point2.y - point1.y) << " " << abs(function(point2) - function(point1)) << endl;

			Gessian(&H0, sk, yk);//проблема тут
			point1 = point2;
			grad1 = grad2;
			check = sqrt(pow(grad1.x, 2) + pow(grad1.y, 2));
			//fprintf(File, "%-5i%-10.7f%-10.7f%-10.7f%-10.7f%-10.7f%-10.7f%-20.7f%-20.7f%-20.7f%-20.7f%\n", i, point1.x, point1.y, function(point1.x,point1.y), -pk.x, -pk.y, lambd, sqrt(pow(point1.x, 2) + pow(point1.y, 2)), ;
			cout << iter << endl;
			cout << point1.x << " " << point1.y << endl;
			cout << function(point1.x, point1.y) << endl;
			cout << -pk.x << " " << -pk.y << endl;
			cout << lambd << endl;
			
			cout << point1 * (-1 * pk) / (norm(point1) * norm(pk))<<endl;
			cout << H0[0][0] << " " <<H0[0][1] << " " << endl;
			cout << H0[1][0] << " " << H0[1][1] << " " << endl;
			cout <<  endl;
		}
	}
};