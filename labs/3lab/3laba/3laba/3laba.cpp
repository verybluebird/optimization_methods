#include"Method_variable_metric.cpp"
#include"Method_0.cpp"
#include"HookJ.cpp"

double r; //коэфф штрафа
double k; //коэфф увеличения штрафа

double g_restriction(point p)
{
    double x = p.x;
    double y = p.y;
    return -x - y + 1;
}
double h_restriction(point p)
{
    double x = p.x;
    double y = p.y;
    return x - y;
    //return (2 * x * x + y * y - 5);
}
double f(point p)
{
    double x = p.x;
    double y = p.y;

    return 5 * (x + y) * (x + y) + (x - 2) * (x - 2);
    //return x + y;
}
double H(double h)
{
    return pow(abs(h),1);
}

double G(double g)
{
    return pow((g + abs(g)) / 2, 1);
}

double Q_b(point p)
{
    return f(p) + r * H(h_restriction(p));
}
double Q_a(point p)
{
    return f(p) + r * G(g_restriction(p));
}

void Penalty_h() {
    HookJ h;
    FILE* out;
    fopen_s(&out, "problem_b.txt", "w");
    point p = point(1,1);
    r = 1;
    k = 2;
    if (out)
    {
        fprintf(out, "|   r     |       (Xi,Yi)        | f(Xi,Yi)| h(Xi,Yi)|\n");
        fprintf(out, "|%- 7f|(%- 7f, %- 7f)|%- 7f|%- 7f|\n", r, p.x, p.y, f(p), h_restriction(p));
    }
    //method_variable_metric create_metod_variable(p, 1e-6, Q_b, "test_b_1e-6.txt");
    //method_0 m = method_0(p, 1e-6, Q_b, "test_b_1e-6.txt");
    h.HookeJeeves(p, Q_b, 1e-3);
    if (out)
    {
        fprintf(out, "|%- 7f|(%- 7f, %- 7f)|%- 7f|%- 7f|\n", r, p.x, p.y, f(p), h_restriction(p));
    }
    for (int i = 0; i < 20; ++i) {
        r *= k;
        //method_variable_metric create_metod_variable(p, 1e-6, Q_b, "test_b_1e-6.txt");
        //method_0 m = method_0(p, 1e-6, Q_b, "test_b_1e-6.txt");
        h.HookeJeeves(p, Q_b, 1e-6);
        
        if (out)
        {
            fprintf(out, "|%- 7f|(%- 7f, %- 7f)|%- 7f|%- 7f|\n", r, p.x, p.y, f(p), h_restriction(p));
        }
        if ( k*abs(h_restriction(p)) < 1e-14)
            break;

    }
    if (out)
        fclose(out);
}

void Penalty_g() {
    HookJ h;
    FILE* out;
    fopen_s(&out, "problem_a.txt", "w");
    point p = point(1,1);
    r = 1;
    k = 2;
    if (out)
    {
        fprintf(out, "|   r     |       (Xi,Yi)        | f(Xi,Yi)| g(Xi,Yi)|\n");
        fprintf(out, "|%- 7f|(%- 7f, %- 7f)|%- 7f|%- 7f|\n", r, p.x, p.y, f(p), g_restriction(p));
    }
    //method_variable_metric create_metod_variable(p, 1e-6, Q_b, "test_b_1e-6.txt");
    //method_0 m = method_0(p, 1e-6, Q_b, "test_b_1e-6.txt");
    h.HookeJeeves(p, Q_a, 1e-3);
    if (out)
    {
        fprintf(out, "|%- 7f|(%- 7f, %- 7f)|%- 7f|%- 7f|\n", r, p.x, p.y, f(p), g_restriction(p));
    }
    for (int i = 0; i < 60; ++i) {
        r *= k;
        //method_variable_metric create_metod_variable(p, 1e-6, Q_b, "test_b_1e-6.txt");
        //method_0 m = method_0(p, 1e-6, Q_b, "test_b_1e-6.txt");
        h.HookeJeeves(p, Q_a, 1e-6);

        if (out)
        {
            fprintf(out, "|%- 7f|(%- 7f, %- 7f)|%- 7f|%- 7f|\n", r, p.x, p.y, f(p), g_restriction(p));
        }
        if (k*g_restriction(p) < 1e-14)
            break;

    }
    if (out)
        fclose(out);
}

int main()
{
    Penalty_h();
    Penalty_g();
}

