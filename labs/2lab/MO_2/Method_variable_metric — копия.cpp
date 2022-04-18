#pragma once
#include <iostream>
#include <cmath>
#include "lib.cpp"
using namespace std;

class method_variable_metric
{
	double delta = 1E-8;
	double diff_x(point point1) //маасив аргументов/по какой переменной/в какую функцию подставляем/шаг
	{
		double h = 10E-7;
		double d0 = point1.x;
		point1.x += h;
		double f_right = function(point1);
		point1.x = d0 - h;
		double f_left = function(point1);
		count_f += 2;
		point1.x = d0;
		return (f_right - f_left) / (2 * h);
	}
	double diff_y(point point1) //маасив аргументов/по какой переменной/в какую функцию подставляем/шаг
	{
		double h = 10E-7;
		double d0 = point1.y;
		point1.y += h;
		double f_right = function(point1);
		point1.y = d0 - h;
		double f_left = function(point1);
		count_f += 2;
		point1.y = d0;
		return (f_right - f_left) / (2 * h);
	}
	point find_grad(point p)
	{
		point grad;
		/*grad.x = 2 * p.x - p.y + 9;
		grad.y = -p.x + 2 * p.y - 6;
		grad.f = function(grad);*/
		grad.x = diff_x(p);
		grad.y = diff_y(p);
		grad.f = function(grad);
		return grad;
	}

	void IntervalSearch(double lymbda, double *a, double *b, point point1, point p0)
	{
		double h = 0;
		double lymbda_previous = lymbda;
		//шаг 1. определяем направление поиска.
		point x1 = point1 + lymbda * p0;
		point x2 = point1 + (lymbda + delta) * p0;
		double f1 = function(x1);
		if ( f1 > function(x2))
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
		double f2 = function(point1 + lymbda * p0);
		count_f += 3;
		while (f1 > f2)
		{
			h *= 2;

			lymbda_previous = lymbda;
			lymbda = lymbda + h;
			f1 = f2;
			f2 = function(point1 + lymbda * p0);
			count_f++;
		}

		*a = lymbda_previous;
		*b = lymbda;
	}
	double GoldenSectionMethod(double a, double b, point point1, point pk)
	{
		int i;
		const double A_COEFF = ((3 - sqrt(5.0)) / 2);
		const double B_COEFF = ((sqrt(5.0) - 1) / 2);
		double x1 = a + A_COEFF * (b - a);
		double x2 = a + B_COEFF * (b - a);
		point p1 = point1 + x1 * pk;
		point p2 = point1 + x2 * pk;
		double  f1 = function(p1);
		double  f2 = function(p2);
		count_f += 2;

		for (i = 0; abs(b - a) > 1e-8; i++)
		{
			if (f1 > f2)
			{
				a = x1;
				x1 = x2;
				f1 = f2;
				x2 = a + B_COEFF * (b - a);
				p2 = point1 + x2 * pk;
				f2 = function(p2);
			}

			else
			{
				b = x2;
				x2 = x1;
				f2 = f1;
				x1 = a + A_COEFF * (b - a);
				p1 = point1 + x1 * pk;
				f1 = function(p1);
			}
			count_f++;
			if (i > 100) break;
		}
		return (a + b) / 2.0;
	}

	double search_lambd(double lambd,point point1, point pk)
	{
		double a = -1, b = 1;
		IntervalSearch(lambd,&a, &b, point1, pk);
		lambd = GoldenSectionMethod(a, b, point1, pk);
		return lambd;
	}
	/*
	_sk - delta xk
	_yk - delta grad fk
	
	*/
	void Gessian(double (&H)[2][2], point dxk, point gk)
	{
		point z = dxk - H * gk;
		double div = z * gk;
		H[0][0] += z.x * z.x/div;
		H[0][1] += z.x * z.y/div;
		H[1][0] += z.y * z.x/div;
		H[1][1] += z.y * z.y/div;
	}

public:
	method_variable_metric(point x)
	{
		int iterations = 0;
		double lambd = 0;
		int n = 2;
		count_f = 0;
		point grad1 = { 0,0 }, grad2, pk, gk, dxk, xk1;

		double H0[2][2] = { 0 };
		
		H0[0][0] = H0[1][1] = 1;//единичная матрица
		grad1 = find_grad(x);

		cout << "gradk1 " << grad1.x << " " << grad1.y << "\n";
		cout << "norm " << norm(grad1) << "\n";
		while (norm(grad1) > eps)
		{
			cout << "\n";
			cout << "iter " << iterations << "\n";
			pk = -1 * (H0 * grad1);
			cout << "pk " << pk.x << " " << pk.y <<"\n";
			cout << "angle " << x*pk / (norm(x)*norm(pk)) << "\n";
			lambd = search_lambd(lambd,x, pk);
			cout << "lambd " << lambd << "\n";
			xk1 = x + lambd * pk;
			cout << "xk1 " << xk1.x << " " << xk1.y << "\n";
 			grad2 = find_grad(xk1);
			cout << "gradk2 " << grad2.x << " " << grad2.y << "\n";
			dxk = xk1 - x;
			cout << "xk1-x " << dxk.x << " " << dxk.y << "\n";
			gk = grad2 - grad1;
			cout << "gk " << gk.x << " " << gk.y << "\n";
			Gessian(H0, dxk, gk);
			cout << "H " << "\n";
			cout << H0[0][0] << " "<< H0[0][1]<<"\n";
			cout << H0[1][0] << " " << H0[1][1] << "\n";
			x = xk1; 
			grad1 = grad2;
			iterations++;
			cout << "norm " << norm(grad1) << "\n";
			if (iterations%n == 0) //обновление алгоритма после n итераций
			{
				H0[0][0] = H0[1][1] = 1;//единичная матрица
				grad1 = find_grad(x);
			}
			if (iterations>50) 
			{
				break;
			}
			cout << x.x << " " << x.y << endl;
			
		}
		cout << "Func: " << function(x) << endl;
		cout << "Iterations: " << iterations << endl;
		cout << "Fcount: " << count_f << endl;
	}
};