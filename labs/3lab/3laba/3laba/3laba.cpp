#include"Method_variable_metric.cpp"

double r; //коэфф штрафа
double k = 2; //коэфф увеличения штрафа
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
    //return 2 * x * x + y * y - 5;
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
    return abs(h);
}

double G(double g)
{
    return (g + abs(g)) / 2;
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
    FILE* out;
    fopen_s(&out, "problem_b.txt", "w");
    point p = point();
    r = 0.01;
    method_variable_metric create_metod_variable(p, 1e-6, Q_b, "test_b_1e-6.txt");
    if (out)
    {
        fprintf(out, "|   r     |       (Xi,Yi)        | f(Xi,Yi)| h(Xi,Yi)|\n");  
        fprintf(out, "|%- 7f|(%- 7f, %- 7f)|%- 7f|%- 7f|\n", r, p.x, p.y, f(p), h_restriction(p));
    }
    while (abs(h_restriction(p)) > 1e-14) {
        r *= k;
        method_variable_metric create_metod_variable(p, 1e-6, Q_b, "test_b_1e-6.txt");
        if (out)
        {
            fprintf(out, "|%- 7f|(%- 7f, %- 7f)|%- 7f|%- 7f|\n", r, p.x, p.y, f(p), h_restriction(p));
        }
    }
    if (out)
        fclose(out);
}

void Penalty_g() {
    FILE* out;
    fopen_s(&out, "problem_a.txt", "w");
    point p = point();
    r = 1;
    method_variable_metric create_metod_variable(p, 1e-6, Q_a, "test_a_1e-6.txt");
    if (out)
    {
        fprintf(out, "|   r     |       (Xi,Yi)        | f(Xi,Yi)| h(Xi,Yi)|\n");
        fprintf(out, "|%- 7f|(%- 7f, %- 7f)|%- 7f|%- 7f|\n", r, p.x, p.y, f(p), h_restriction(p));
    }
    while (g_restriction(p) > 1e-14) {
        r *= k;
        method_variable_metric create_metod_variable(p, 1e-6, Q_a, "test_a_1e-6.txt");
        if (out)
        {
            fprintf(out, "|%- 7f|(%- 7f, %- 7f)|%- 7f|%- 7f|\n", r, p.x, p.y, f(p), h_restriction(p));
        }
    }
    if (out)
        fclose(out);
}

int main()
{
    Penalty_h();
    Penalty_g();
}

