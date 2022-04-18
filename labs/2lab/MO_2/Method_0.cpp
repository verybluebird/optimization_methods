#pragma once
#include "lib.cpp"

class method_0
{
    int n = 2;
    OneDimensionFunction function;
    void sort(point* point_array)
    {
        point tmp;
        for (int i = 0; i < 3; i++) {
            for (int j = 2; j >= (i + 1); j--) {
                if (point_array[j].f < point_array[j - 1].f) {
                    tmp = point_array[j];
                    point_array[j] = point_array[j - 1];
                    point_array[j - 1] = tmp;
                }
            }
        }
    }

    double d1(double t)
    {
        return (t / (n * sqrt(2)) * (sqrt(n + 1) + n - 1));
    }

    double d2(double t)
    {
        return (t / (n * sqrt(2) * (sqrt(n + 1) - 1)));
    }

    point create_point(point point0, point vector)//ńîçäŕíčĺ äîď 2 ňî÷ĺę
    {
        point0.x = point0.x + 0.05 * vector.x;
        point0.y = point0.y + 0.05 * vector.y;
        return point0;
    }

    point elongation(point new_point, point mid, double gamma)//đŕńň˙ćĺíčĺ
    {
        new_point.x = mid.x + gamma * (new_point.x - mid.x);
        new_point.y = mid.y + gamma * (new_point.y - mid.y);
        new_point.f = function(new_point);
        return new_point;
    }

    void reduction(point* point_array)
    {
        for (int i = 1; i < 3; i++)
        {
            point_array[i].x = point_array[0].x + 0.5 * (point_array[i].x - point_array[0].x);
            point_array[i].y = point_array[0].y + 0.5 * (point_array[i].y - point_array[0].y);
        }
        point_array[1].f = function(point_array[1]);
        point_array[2].f = function(point_array[2]);
    }
    point compression(point worst, point mid, double betta)
    {
        worst.x = mid.x + betta * (worst.x - mid.x);
        worst.y = mid.y + betta * (worst.y - mid.y);
        worst.f = function(worst);
        count_f++;
        return worst;
    }

    void reflection(double alpha, double betta, double gamma, point mid, point* point_array) //îňđŕćĺíčĺ őóäřĺé ňî÷ęč îňíîńčňĺëüíî mid
    {
        point new_point, point;//íŕäî ëč point âîîáůĺ?
        new_point.x = mid.x + alpha * (mid.x - point_array[2].x);
        new_point.y = mid.y + alpha * (mid.y - point_array[2].y);
        new_point.f = function(new_point);
        count_f++;
        while (true)
        {
            if (new_point.f >= point_array[2].f)
            {
                reduction(point_array);//đĺäóęöč˙
                break;
            }
            if (new_point.f <= point_array[0].f)//ĺńëč ěĺíüřĺ ěčíčěŕëüíîăî íŕ k-ýňŕďĺ
            {
                point = elongation(new_point, mid, gamma);//đŕńň˙ćĺíčĺ
                if (point.f < point_array[0].f)
                {
                    point_array[2] = point_array[1];
                    point_array[1] = point_array[0];
                    point_array[0] = point;
                }
                else
                {
                    point_array[2] = point_array[1];
                    point_array[1] = point_array[0];
                    point_array[0] = new_point;
                }
                break;
            }
            if (new_point.f > point_array[1].f && new_point.f <= point_array[2].f)//ěĺćäó äâóě˙ őóäřčěč ňî÷ęŕěč
            {
                point_array[2] = compression(point_array[2], mid, betta);//ńćŕňčĺ íŕäî ëč âîçâđŕůňü âîîáůĺ??
                break;
            }
            if (new_point.f > point_array[0].f && new_point.f <= point_array[1].f)
            {
                point_array[2] = new_point;
                break;
            }
        }
    }

public:
    method_0(point point1, double eps, OneDimensionFunction function, std::string fname)
    {
        this->function = function;
        point point_array[3], mid, point_prew_array[3];
        int iter;
        double flag1;
        double alpha = 1;
        double betta = 0.5;
        double gamma = 2;
        bool flag_i = true, flag_f = true;
        int k = 0;
        double t = 0.1;
        point_array[0] = point1;
        point_array[1] = { point_array[0].x + d1(t),point_array[0].y + d2(t) };
        point_array[2] = { point_array[0].x + d2(t),point_array[0].y + d1(t) };
        point_array[0].f = function(point_array[0]);
        point_array[1].f = function(point_array[1]);
        point_array[2].f = function(point_array[2]);
       
        count_f = 3;
        for (iter = 1; flag_f; iter++)
        {

            sort(point_array);
            for (int i = 0; i < 3; i++) point_prew_array[i] = point_array[i];
            mid = { (point_array[0].x + point_array[1].x) / 2, (point_array[0].y + point_array[1].y) / 2 };
            mid.f = function(mid);
            count_f++;
            reflection(alpha, betta, gamma, mid, point_array);
            flag1 = sqrt(pow(point_array[0].f, 2) - pow(mid.f, 2)) / 3;
            if (flag1 < eps_f) flag_f = false;
            std::cout << iter << ": Direction: (" << mid.x - point_array[2].x << " , " << mid.x - point_array[2].x << ") " << std::endl;
            for (int i = 0; i < 3; i++) std::cout << std::setprecision(11) << " (" << abs(point_array[i].x - point_prew_array[i].x) << ", " << abs(point_array[i].y - point_prew_array[i].y) << " , " << abs(point_array[i].f - point_prew_array[i].f) << ")  " << std::endl;
            std::cout << std::setprecision(11) << "(" << point_array[0].x << ", " << point_array[0].y << ") " << point_array[0].f << std::endl;

        }
        std::cout << std::setprecision(11) << point_array[0].x << " " << point_array[0].y << " " << point_array[0].f << std::endl;
        std::cout << "Iter: " << iter;
    
    }
};



