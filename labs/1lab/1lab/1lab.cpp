#include <iostream>
#include<vector>
#include <fstream>
#include <iomanip> 
using namespace std;

double a = -2;
double b = 20;

double eps = 0.0001;

double f(double x) {
    return (x - 8) * (x - 8);
}
//-----------------------------------------------------------------------------
/* Метод дихотомии для поиска минимума. */
void calcDichotomy(double a, double b) {
    std::ofstream fout("Dichotomy.txt");
    fout << std::setprecision(11);
    cout << std::setprecision(11);
    double x1, x2, f1, f2;
    double delta = eps / 10;
    int fcount = 0; //сколько раз посчиталась ф-я
    for (int i = 0; abs(b - a) > eps; i++)
    {
        x1 = (b + a - delta) / 2;
        x2 = (b + a + delta) / 2;
        f1 = f(x1);
        f2 = f(x2);
        if (f1 > f2)
            a = x1;
        else
            b = x2;
        fcount += 2;
        fout << "Iteration: " << i << "  a: " << a << "  b: " << b << "  X_min: " << (a + b) / 2 << "  fcount: " << fcount << endl;
        cout << "Iteration: " << i << "  a: " << a << "  b: " << b << "  X_min: " << (a + b) / 2 << "  fcount: " << fcount << endl;
    }
}
//-----------------------------------------------------------------------------
/* Метод золотого сечения для поиска минимума. */
void calcGoldenRatio(double a, double b) {
    std::ofstream fout("GoldenRatio.txt");

    int i;
    int fcount = 0;
    const double A_COEFF((3 - sqrt(5.0)) / 2);
    const double B_COEFF((sqrt(5.0) - 1) / 2);
    double x1 = a + A_COEFF * (b - a);
    double x2 = a + B_COEFF * (b - a);
    double  f1 = f(x1);
    double  f2 = f(x2);
    fcount += 2;

    for (i = 0; abs(b - a) > eps; i++)
    {
        if (f1 > f2)
        {
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + B_COEFF * (b - a);
            f2 = f(x2);
        }

        else {
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = a + A_COEFF * (b - a);
            f1 = f(x1);
        }
        fcount++;
        fout << "Iteration: " << i << "  a: " << a << "  b: " << b << "  X_min: " << (a + b) / 2 << "  fcount: " << fcount << endl;
        cout << "Iteration: " << i << "  a: " << a << "  b: " << b << "  X_min: " << (a + b) / 2 << "  fcount: " << fcount << endl;
    }
}
//-----------------------------------------------------------------------------
/* Функция для нахождения интервала, содержащего минимум для унимодальной функции. */
void findInterval(double a, double b, double x0)
{
    std::ofstream fout("Interval.txt");
    int fcount;
    double delta = eps / 10;
    double f0 = f(x0);
    double f1 = f(x0 + delta);
    double x1, h, x00;
    int k = 0;
    if (f0 > f1) //Если ф-я на участке убывает
    {
        k = 1;
        h = delta;
        
    }
    else //Если ф-я на участке возрастает
        h = -delta;
    
    x1 = x0 + h;
    do {
        h *= 2;
        x00 = x0;
        x0 = x1;
        x1 = x0 + h;
        f0 = f1;
        f1 = f(x1);
        k++;
        /*if (x00 > x1)
        {
            fout << x1 << " " << x00 << " iters " << k << endl;
            cout << x1 << " " << x00 << " iters " << k << endl;
        }
        else
        {
            fout << x00 << " " << x1 << " iters " << k << endl;
            cout << x00 << " " << x1 << " iters " << k << endl;
        }*/
    } while (f1 < f0);

    a = x00;
    b = x1;
    if (x00 > x1)
    {
        fout << x1<< " " << x00 << " iters " << k << endl;
        cout << x1 << " " << x00 << " iters " << k << endl;
    }
    else
    {
        fout << x00 << " " << x1 << " iters " << k << endl;
        cout << x00 << " " << x1 << " iters " << k << endl;
    }
}
//-----------------------------------------------------------------------------
/** Метод Фибоначчи для поиска минимума. */
void calcFibonacci(double a, double b)
{
    std::ofstream fout("Fibonacci.txt");

    int fcount = 0;
    int i;
    double x1, x2, f1, f2;
    double n = 2, max = (b - a) / eps, new_number = 0;//?
    vector<int> fibonacci_num;//массив чисел фибоначи 

    fibonacci_num.push_back(1);
    fibonacci_num.push_back(1);
    for (; max > new_number; n++)//заполняем массив числами фибоначи
    {
        new_number = fibonacci_num[n - 1] + fibonacci_num[n - 2];
        fibonacci_num.push_back(new_number);
    }
    n = fibonacci_num.size() - 3; //1 число привысившее максимум, 2 числа для использования формулы n+2

    x1 = a + fibonacci_num[n] * (b - a) / fibonacci_num[n + 2];
    x2 = a + fibonacci_num[n + 1] * (b - a) / fibonacci_num[n + 2];

    f1 = f(x1);
    f2 = f(x2);
    fcount += 2;
    for (i = 0; i < n - 2; i++) {
        if (f1 > f2)
        {
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + fibonacci_num[n - i - 1] * (b - a) / fibonacci_num[n - i];
            f2 = f(x2);
        }
        else
        {
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = a + fibonacci_num[n - i - 2] * (b - a) / fibonacci_num[n - i];
            f1 = f(x1);
        }
        fcount++;
        fout << "Iteration: " << i << "  a: " << a << "  b: " << b << "  X_min: " << (a + b) / 2 << "  fcount: " << fcount << endl;
        cout << "Iteration: " << i << "  a: " << a << "  b: " << b << "  X_min: " << (a + b) / 2 << "  fcount: " << fcount << endl;
    }
    fout << "Iteration: " << i << "  a: " << a << "  b: " << b << "  X_min: " << (a + b) / 2 << "  fcount: " << fcount << endl;
    cout << "Iteration: " << i << "  a: " << a << "  b: " << b << "  X_min: " << (a + b) / 2 << "  fcount: " << fcount << endl;

}
int main()
{
    float x0;
    cout << "Enter the initial value: " << endl << "x0 = ";
    cin >> x0;
    cout << "findInterval: " << endl;
    findInterval(a, b, x0);
    cout << "calcDichotomy: " << endl;
    calcDichotomy(a, b);
    cout << "calcGoldenRatio: " << endl;
    calcGoldenRatio(a, b);
    cout << "calcFibonacci: " << endl;
    calcFibonacci(a, b);
}