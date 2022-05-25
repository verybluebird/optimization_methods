#pragma once
#include "lib.h"


class method_variable_metric
{
	OneDimensionFunction function;
	double delta = 1E-8;
	double diff_x(point point1) //маcсив аргументов/по какой переменной/в какую функцию подставляем/шаг
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
	double diff_y(point point1) //маcсив аргументов/по какой переменной/в какую функцию подставляем/шаг
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
		grad.x = diff_x(p);
		grad.y = diff_y(p);
		grad.f = function(grad);
		return grad;
	}

	void IntervalSearch(double lymbda, double* a, double* b, point point1, point p0)
	{
		double h = 0;
		double lymbda_previous = lymbda;
		//шаг 1. определяем направление поиска.
		point x1 = point1 + lymbda * p0;
		point x2 = point1 + (lymbda + delta) * p0;
		double f1 = function(x1);
		if (f1 > function(x2))
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

	double search_lambd(double lambd, point point1, point pk)
	{
		double a = -1, b = 1;
		IntervalSearch(lambd, &a, &b, point1, pk);
		lambd = GoldenSectionMethod(a, b, point1, pk);
		return lambd;
	}
	
	void Gessian(double (&H)[2][2], point dxk, point gk)
	{
		point z = dxk - H * gk;
		double div = z * gk;
		/*if (div < 1e-19)
			return;*/
		H[0][0] += z.x * z.x / div;
		H[0][1] += z.x * z.y / div;
		H[1][0] += z.y * z.x / div;
		H[1][1] += z.y * z.y / div;
	}

public:
	int count_f;
	method_variable_metric(point& x, double eps, OneDimensionFunction function, std::string fname)
	{
		this->function = function;
		FILE* out;
		FILE* out2;
		const char* c = fname.c_str();
		const char* c1 = "trajectory.txt";
		fopen_s(&out, c, "w");
		fopen_s(&out2, c1, "w");
		if (out)
		{
			fprintf(out, "| ITERS |       (Xi,Yi)        | f(Xi,Yi)|       (s1,s2)        |  Lambd  |/Xi-Xi-1/|  angle  |           H           |\n");
			fprintf(out, "|       |                      |         |                      |         |/Yi-Yi-1/|         |                       |\n");
			fprintf(out, "|       |                      |         |                      |         |/fi-fi-1/|         |                       |\n");
		}
		int iterations = 0;
		double lambd = 0;
		int n = 2;
		double angle;
		double df;
		count_f = 0;
		point grad1, grad2, pk, gk, dxk, xk1;
		double H0[2][2] = { 0 };
		H0[0][0] = H0[1][1] = 1;//единичная матрица
		grad1 = find_grad(x);
		if (out) {
			fprintf(out, "|%- 7i|(%- 7f, %- 7f)|%- 7f|                      |         |         |         |                       |\n", iterations, x.x, x.y, function(x));
		}
		if (out2) {
			fprintf(out2, "%- 7f	%- 7f\n", x.x, x.y);
		}
		while (norm(grad1) > eps)
		{
			
			pk = -1 * (H0 * grad1);
			lambd = search_lambd(lambd, x, pk);
			xk1 = x + lambd * pk;
			grad2 = find_grad(xk1);
			dxk = xk1 - x;
			gk = grad2 - grad1;
			Gessian(H0, dxk, gk);
			df = function(xk1) - function(x);
			x = xk1;
			angle = x * pk / (norm(x) * norm(pk));
			iterations++;
			if (out) {
				fprintf(out, "|%- 7i|(%- 7f, %- 7f)|%- 7f|(%- 7f, %- 7f)|%- 7f|%- 7f|%- 7f|[%- 7f   %- 7f]|\n", iterations, xk1.x, xk1.y, function(xk1), -pk.x, -pk.y, lambd, abs(dxk.x), angle, H0[0][0], H0[0][1]);
				fprintf(out, "|       |                      |         |                      |         |%- 7f|         |[%- 7f   %- 7f]|\n",  abs(dxk.y), H0[1][0], H0[1][1]);
				fprintf(out, "|       |                      |         |                      |         |%- 7f|         |                       |\n", abs(df));
			}
			if (out2) {
				fprintf(out2, "%- 7f	%- 7f\n", xk1.x, xk1.y);
			}

			grad1 = grad2;
			
			//if (iterations % n == 0) //обновление алгоритма после n итераций
			//{
			//	H0[0][0] = H0[1][1] = 1;//единичная матрица
			//	grad1 = find_grad(x);
			//}
			if (iterations > 50)
				break;
			std::cout << x.x << " " << x.y << std::endl;
		}
		if (out)
			fclose(out);
		if (out2)
			fclose(out2);
		std::cout << "Func: " << function(x) << std::endl;
		std::cout << "Iterations: " << iterations << std::endl;
		std::cout << "Fcount: " << count_f << std::endl;
		std::cout << std::endl;
	}
};