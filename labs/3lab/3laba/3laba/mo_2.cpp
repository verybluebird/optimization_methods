#include"Method_0.cpp"
#include"Method_variable_metric.cpp"
#include"lib.cpp"

int main()
{
    point point;
    std::cout << "Enter x0: ";
    std::cin >> point.x ;
    std::cout << "Enter y0: ";
    std::cin >> point.y;

    std::vector<test> Tests = {
    {f1, "table1_1e-3.txt",1e-3},
    {f2, "table2_1e-3.txt",1e-3},
    {f3, "table3_1e-3.txt",1e-3},
    {f1, "table1_1e-6.txt",1e-6},
    {f2, "table2_1e-6.txt",1e-6},
    {f3, "table3_1e-6.txt",1e-6},
    };

    for (test t :  Tests) {
        std::cout << "Broyden's\n";
        method_variable_metric create_metod_variable(point, t.eps, t.func, t.fname);
        std::cout << "Simplex\n";
        method_0 create_method0 (point, t.eps, t.func, t.fname);
    }
    
    std::cout << std::endl << count_f << std::endl;
}
