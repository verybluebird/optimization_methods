#include"Method_variable_metric.cpp"
#include"Method_0.cpp"
#include"HookJ.cpp"
#include <random>
#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

double bestSolution = 0;
struct RandMethodsResult {
	int fcount;
	double value;
	point p;
	int N; // for simple Random Search
};

double func1(point v) {
	vector<int> c = { 2,1,7,2,8,4 };
	vector<int> a = { 5,2,-9,0,-3,-3 };
	vector<int> b = { 4,0,-6,-3,7,3 };
	

#define sqr(a) ((a)*(a))
	double sum = 0;
	for (int i = 0; i < 6; i++) {
		sum += c[i] / (1.0 + sqr(v.x - a[i]) + sqr(v.y - b[i]));
	}
	return -sum;
}

RandMethodsResult simpleRandomSearch(const double eps, const double p, const double x0, const double x1, const double y0, const double y1) {
	std::default_random_engine generator(0);
	std::uniform_real_distribution<double> xaxis(x0, x1);
	std::uniform_real_distribution<double> yaxis(y0, y1);

	int fcount = 0;
	// Calculated necessary amount of points:
	double p_eps = fabs(eps * eps / ((x1 - x0) * (y1 - y0)));

	int64_t N = log(1.0 - p) / log(1.0 - p_eps) + 1; 
	cout << N << endl;

	double valuemin = std::numeric_limits<double>::infinity();
	point pointmin;

	for (int i = 0; i < N; i++) {
		point p;
		p.x = xaxis(generator);
		p.y = yaxis(generator);

		double value = func1(p);
		fcount++;

		if (value < valuemin) {
			valuemin = value;
			pointmin = p;
		}
	}
	
	RandMethodsResult result = { fcount, valuemin, pointmin, N };
	return result;
}

void GetSolveByHookJ(point& p_0, const double eps, const double x0, const double x1, const double y0, const double y1, int& fcount)
{
	std::default_random_engine generator(0);
	std::uniform_real_distribution<double> xaxis(x0, x1);
	std::uniform_real_distribution<double> yaxis(y0, y1);
	p_0.x = xaxis(generator);
	p_0.y = yaxis(generator);
	//point p = p_0;
	HookJ h;
	h.HookeJeeves(p_0, func1, eps, fcount);
}

RandMethodsResult Alg1(int m, const double eps, const double x0, const double x1, const double y0, const double y1)
{
	int fcount = 0;
	point p_0, pointmin;
	double value;
	GetSolveByHookJ(pointmin, eps, x0, x1, y0, y1, fcount);
	double valuemin = func1(pointmin);
	fcount++;
	
	
	for (int i = 0; i < m; i++)
	{
		GetSolveByHookJ(p_0, eps, x0, x1, y0, y1, fcount);
		value = func1(p_0);
		fcount++;
		if (value < valuemin)
		{
			valuemin = value;
			pointmin = p_0;
		}
	}
	RandMethodsResult result = { fcount, valuemin, pointmin, m };
	return result;
}




RandMethodsResult Alg2(int m, const double eps, const double x0, const double x1, const double y0, const double y1)
{
	int fcount = 0;
	point p_0, pointmin;
	double value;
	// find first local minimum point
	GetSolveByHookJ(pointmin, eps, x0, x1, y0, y1, fcount);
	double valuemin = func1(pointmin);
	fcount++;
	
	for (int i = 0; i < m; i++)
	{
		// new random point
		std::default_random_engine generator(0);
		std::uniform_real_distribution<double> xaxis(x0, x1);
		std::uniform_real_distribution<double> yaxis(y0, y1);
		p_0.x = xaxis(generator);
		p_0.y = yaxis(generator);
		value = func1(p_0);
		fcount++;

		if (value < valuemin)
		{
			pointmin = p_0;
			GetSolveByHookJ(pointmin, eps, x0, x1, y0, y1, fcount);
			valuemin = func1(pointmin);
			i = 0;//сбросили число попыток
			fcount++;
		}
	}
	RandMethodsResult result = { fcount, valuemin, pointmin, m };
	return result;
}

void testAlgs() {
	std::ofstream fout("research2.txt");
	std::vector<int> ms = { 1, 3, 5, 7, 10, 50, 100, 500, 1000, 5000, 10000 };

	const double x0 = -10, x1 = 10;
	const double y0 = -10, y1 = 10;
	double eps = 0.00001;
	fout << "m\ta1_fcount\ta2_fcount\ta1_value\ta2_value\ta1_prec\ta2_prec\n";
	for (int maxMissSearch : ms) {
		RandMethodsResult a1 = Alg1(maxMissSearch,eps, x0, x1, y0, y1);
		RandMethodsResult a2 = Alg2(maxMissSearch,eps, x0, x1, y0, y1);

		fout << maxMissSearch << "\t" << a1.fcount << "\t" << a2.fcount
			<< "\t" << a1.value << "\t" << a2.value
			<< "\t" << a1.value - bestSolution << "\t" << a2.value - bestSolution;
		// precision of found result!?
		fout << std::endl;

	}
	fout.close();
}

void testSimpleRandom() {
	std::ofstream fout("research1.txt");
	std::vector<double> epss = { 1, 0.5, 0.1, 0.05, 0.01, 0.005 };
	std::vector<double> ps = { 0.1, 0.3, 0.5, 0.7, 0.9 };

	const double x0 = -10, x1 = 10;
	const double y0 = -10, y1 = 10;

	fout << "eps\tp\tN\t(x*,y*)\tf(x*,y*)\n";

	for (double eps : epss) {
		fout << eps << "\t";
		for (double p : ps) {
			fout << p << "\t";
			RandMethodsResult r = simpleRandomSearch(eps, p, x0, x1, y0, y1);
			fout << r.N << "\t" << r.p.x << "\t" << r.p.y << "\t" << r.value;
			fout << std::endl;
			fout << "\t";
		}
		
	}
	fout.close();
}

int main() {
	testSimpleRandom();
	const double x0 = -10, x1 = 10;
	const double y0 = -10, y1 = 10;
	int fcount = 0;
	point pointmin;
	GetSolveByHookJ(pointmin, 1e-14, x0, x1, y0, y1, fcount);
	bestSolution = func1(pointmin);
	//testAlgs();
	return 0;
}