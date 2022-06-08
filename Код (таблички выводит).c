#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define Y0 1.0
#define X0 1.0
#define N 5

// h - шаг дифференцирования. Я рассматриваю 3 случая (h=0,5; h=0,1; h=0,01).
// Первая строка таблицы - значения х. Оно возрастает с указанным выше шагом.
// Я исследую методы численного дифференцирования на промежутке [1;6]

// Со второй по пятую строки таблицы - значения y по 
//Методу Эйлера, Модифицированному методу Эйлера, Усовершенствованному методу Эйлера, Методу Рунге-Кутта каждая соответсвено.
// Если что, названия методов указанны прямо под строками их значений y.


double F (double x, double y){
    return (x*y-2)/pow(x,2);
}
double F1(double x, double y){
    return (x*y-2*(y*x-2))/pow(x,3);
}

void Euler_original(double h){
    FILE*file;
    file = fopen("file.txt", "a");
    
    double x = X0;
    fprintf(file, "\nh = %lf\n", h);
    while (x <= X0 + N){
        fprintf(file, "%20.5lf ", x);
        x += h;
    }

    fprintf(file, "\n");
    double y = Y0;
    x = X0;
    while (x <= X0 + N){
        fprintf(file, "%20.5lf ", y);
        y = y + h * F(x,y);
        x += h;
    }
    fprintf(file, "%20.5lf ", y);
    fprintf(file, "\n");
    fprintf(file, "Euler_original");
    fprintf(file, "\n");
    fclose(file);
}


void method_euler_mod(double h){
    FILE*file;
    file = fopen("file.txt", "a");
    double y = Y0, x = X0,xi, yi;
    while (x <= X0 + N){
        fprintf(file, "%20.5lf ", y);
        xi = x + h/2;
        yi =y + h/2*F1(x,y);
        y = y + h * F(xi,yi);
        x += h;
    }
    fprintf(file, "%20.5lf ", y);
    fprintf(file, "\n");
    fprintf(file, "euler_mod");
    fprintf(file, "\n");
    fclose(file);
}


void method_euler_improve(double h){
    FILE*file;
    file = fopen("file.txt", "a");
    double y = Y0, x = X0, yi;
    while (x <= X0 + N){
        fprintf(file, "%20.5lf ", y);
        yi = y + h*F1(x,y);
        x += h;
        y = y + h/2 * (F(x-h,y)+ F(x,yi));
        }
    fprintf(file, "%20.5lf ", y);
    fprintf(file, "\n");
    fprintf(file, "euler_improve");
    fprintf(file, "\n");
    fclose(file);
}


void method_runge_kutt(double h ){
    FILE*file;
    file = fopen("file.txt", "a");
    double y = Y0, x = X0, k1, k2, k3, k4;
    while (x <= X0 + N){
        fprintf(file, "%20.5lf ", y);
        k1 = h*F(x,y);
        k2 = h*F(x+h/2, y+k1/2*h);
        k3 = h*F(x+h/2, y+k2/2*h);
        k4 = h*F(x+h, y + k3*h);
        y = y + (k1 + 2*k2 + 2*k3+ k4)/6;
        x = x + h;
    }
    fprintf(file, "%20.5lf ", y);
    fprintf(file, "\n");
    fprintf(file, "runge_kutt");
    fprintf(file, "\n");
    fclose(file);
}




int main(void) {
    FILE*file;
    file = fopen("file.txt", "w");
    fclose(file);
    double h, step[3] = {0.5, 0.1, 0.01};
    int i;
    for(i = 0; i < 3; i++) {
        h = step[i];
        Euler_original(h);
        method_euler_mod(h);
        method_euler_improve(h);
        method_runge_kutt(h);
    }
    return EXIT_SUCCESS;
}
