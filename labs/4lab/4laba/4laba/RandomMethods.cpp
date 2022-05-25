#include"Method_variable_metric.cpp"
#include"Method_0.cpp"
#include"HookJ.cpp"
#include <random>
#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

struct RandMethodsResult {
	int fcount;
	double value;
	point p;
	int N; // for simple Random Search
};

double func1(point v) {
	vector<int> a = { 5,2,-9,0,-3,-3 };
	vector<int> b = { 4,0,-6,-3,7,3 };
	vector<int> c = { 2,1,7,2,8,4 };

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

		double value = f1(p);
		fcount++;

		if (value < valuemin) {
			valuemin = value;
			pointmin = p;
		}
	}
	
	RandMethodsResult result = { fcount, valuemin, pointmin, N };
	return result;
}

RandMethodsResult Alg1(const double eps, const double p, const double x0, const double x1, const double y0, const double y1)
{
	std::default_random_engine generator(0);
	std::uniform_real_distribution<double> xaxis(x0, x1);
	std::uniform_real_distribution<double> yaxis(y0, y1);

	for (int i = 0; i < N; i++) {
		point p_0;
		p_0.x = xaxis(generator);
		p_0.y = yaxis(generator);

		point p = p_0;
		HookJ h;
		h.HookeJeeves(p, f1, eps);
		double value = f1(p);
		fcount++;

		if (value < valuemin) {
			valuemin = value;
			pointmin = p;
		}
	}

}
void testSimpleRandom() {
	std::ofstream fout("research1.txt");
	std::vector<double> epss = { 1, 0.5, 0.1, 0.05, 0.01, 0.005 };
	std::vector<double> ps = { 0.1, 0.3, 0.5, 0.7, 0.9 };

	const double x0 = -10, x1 = 10;
	const double y0 = -10, y1 = 10;

	fout << "\t";
	for (double p : ps) fout << p << "\t";
	fout << std::endl;

	for (double eps : epss) {
		fout << eps << "\t";
		for (double p : ps) {
			RandMethodsResult r = simpleRandomSearch(eps, p, x0, x1, y0, y1);
			fout << r.N << " " << r.value << " " << r.p.x << " " << r.p.y << "\t";
		}
		fout << std::endl;
	}
	fout.close();
}

int main() {
	testSimpleRandom();
	return 0;
}