#pragma once
#include "lib.h"


class HookJ
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

	

public:
	int count_f = 0;
	int iterations = 0;
	void HookeJeeves(point& x0, OneDimensionFunction function, double eps) {
		this->function = function;
		point x = x0;
		point x1 = x0;
		point s = x0;// direction to 1D minimization

		double f0 = function(x0), flast = f0;// start f value
		double f1val = 0;// value in next finding point

		// optimization:
		while (true) {
			// examining search
			bool succesfulStep = false;
			double dx = 1e-2;		// start dx value

			int localIterations = 0;
			do {
				//шаг по x
				double fp = 0;			// value in x + dx point
				double fm = 0;			// value in x - dx point


				double temp = x1.x;
				x1.x += dx;
				fp = function(x1);

				if (fp > f0) {
					x1.x = temp - dx;
					fm = function(x1);

					if (fm > f0) x1.x = temp;	// x1[i] not changed
					else f0 = fm;
				}
				else f0 = fp;
				//шаг по y
				fp = 0;			// value in x + dx point
				fm = 0;			// value in x - dx point


				temp = x1.y;
				x1.y += dx;
				fp = function(x1);

				if (fp > f0) {
					x1.y = temp - dx;
					fm = function(x1);

					if (fm > f0) x1.y = temp;	// x1[i] not changed
					else f0 = fm;
				}
				else f0 = fp;

				if (norm(x1 - x) < 1e-13) dx /= 2;
				else succesfulStep = true;

				localIterations++;
				if (localIterations > 100) break;

			} while (!succesfulStep);

			// minimization in finded direction
			s = x1 - x;		// direction

			double lambda = 0;
			lambda = search_lambd(lambda, x, s);
			x1 = x + lambda * s; 

			f1val = function(x1);

			iterations++;

			
			double sub = norm(x - x1);
			double fsub = fabs(f1val - flast);
			double gg = fsub / sub;
			if (gg < eps || iterations > 100) break;
			else {	// prepare next iteration:
				f0 = f1val;
				flast = f0;
				x = x1;
			}

		}

		x0 = x1;

	}
};