#ifndef POLYNOMIALS_H
#define POLYNOMIAL_H

#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <complex>
#include <numbers>
#include <algorithm>
#include <fstream>
#include <queue>

using namespace std;

/*
    Имплементация метода решения линейного уравнения
*/
template<typename fp_t>
unsigned int solveLinear(fp_t a, fp_t b, vector<fp_t>& roots);

/*
    Имплементация метода решения линейного уравнения (версия с комплексной перегрузкой, для вещественных коэффициентов)
*/
template<typename fp_t>
unsigned int solveLinear(fp_t a, fp_t b, vector<complex<fp_t>>& roots);

/*
    Имплементация метода решения линейного уравнения (версия с комплексной перегрузкой, для комплексных коэффициентов)
*/
template<typename fp_t>
unsigned int solveLinear(complex<fp_t> a, complex<fp_t> b, vector<complex<fp_t>>& roots);

/*
    Имплементация метода решения квадратного уравнения
*/
template<typename fp_t>
unsigned int solveQuadratic(fp_t n, fp_t a, fp_t b, vector<fp_t>& roots);

/*
    Имплементация метода решения квадратного уравнения (версия с комплексной перегрузкой, для вещественных коэффициентов)
*/
template<typename fp_t>
unsigned int solveQuadratic(fp_t n, fp_t a, fp_t b, vector<complex<fp_t>>& roots);

/*
    Имплементация метода решения квадратного уравнения (версия с комплексной перегрузкой, для комплексных коэффициентов)
*/
template<typename fp_t>
unsigned int solveQuadratic(complex<fp_t> n, complex<fp_t> a, complex<fp_t> b, vector<complex<fp_t>>& roots);

/*
    Имплементация метода решения кубического уравнения — Numerical Recipes / Viète
    Информация о методе — https://quarticequations.com/Tutorial.pdf (стр. 1)
    Работу выполнил — Погосов Даниэль (https://github.com/DarklleS)
*/
template<typename fp_t>
unsigned int solveCubic(fp_t n, fp_t a, fp_t b, fp_t c, vector<fp_t>& roots);

/*
    Имплементация метода решения кубического уравнения — Numerical Recipes / Viète (версия с комплексной перегрузкой)
    Информация о методе — https://quarticequations.com/Tutorial.pdf (стр. 1)
    Работу выполнил — Погосов Даниэль (https://github.com/DarklleS)
*/
template<typename fp_t>
unsigned int solveCubic(fp_t n, fp_t a, fp_t b, fp_t c, vector<complex<fp_t>>& roots);

/*
    Имплементация метода решения уравнения четвертой степени — Ferrari's Method
    Информация о методе — https://quarticequations.com/Quartic2.pdf (стр. 3)
    Работу выполнил — Погосов Даниэль (https://github.com/DarklleS)
*/
template<typename fp_t>
unsigned int ferrari(fp_t n, fp_t a, fp_t b, fp_t c, fp_t d, vector<fp_t>& roots);

/*
    Имплементация метода решения уравнения четвертой степени — Descartes’ Method
    Информация о методе — https://quarticequations.com/Quartic2.pdf (стр. 4)
    Работу выполнил — Погосов Даниэль (https://github.com/DarklleS)
*/
template<typename fp_t>
unsigned int descartes(fp_t n, fp_t a, fp_t b, fp_t c, fp_t d, vector<fp_t>& roots);

/*
    Имплементация метода решения уравнения четвертой степени — NBS Method
    Информация о методе — https://quarticequations.com/Quartic2.pdf (стр. 5)
    Работу выполнил — Погосов Даниэль (https://github.com/DarklleS)
*/
template<typename fp_t>
unsigned int nbs(fp_t n, fp_t a, fp_t b, fp_t c, fp_t d, vector<fp_t>& roots);

/*
    Имплементация метода решения уравнения четвертой степени — Euler’s Method
    Информация о методе — https://quarticequations.com/Quartic2.pdf (стр. 6)
    Работу выполнил — Погосов Даниэль (https://github.com/DarklleS)
*/
template<typename fp_t>
unsigned int euler(fp_t n, fp_t a, fp_t b, fp_t c, fp_t d, vector<fp_t>& roots);

/*
    Имплементация метода решения уравнения четвертой степени — Van der Waerden’s Method
    Информация о методе — https://quarticequations.com/Quartic2.pdf (стр. 7)
    Работу выполнил — Погосов Даниэль (https://github.com/DarklleS)
*/
template<typename fp_t>
unsigned int vanDerWaerden(fp_t n, fp_t a, fp_t b, fp_t c, fp_t d, vector<fp_t>& roots);

/*
    Имплементация метода решения уравнения четвертой степени — Van der Waerden’s Method
    Информация о методе — https://dialnet.unirioja.es/descarga/articulo/6523979.pdf
    Работу выполнил — Погосов Даниэль (https://github.com/DarklleS)
*/
template<typename fp_t>
unsigned int tschirnhaus(fp_t n, fp_t a, fp_t b, fp_t c, fp_t d, vector<fp_t>& roots);


// Комплексная версия алгоритма Феррари
template<typename fp_t>
unsigned int ferrari_complex(fp_t n, fp_t a, fp_t b, fp_t c, fp_t d, vector<complex<fp_t>>& roots);

/*
    Имплементация метода решения уравнения четвертой степени — FQS method
    Информация о методе — https://www.academia.edu/es/27122162/The_fast_quartic_solver
    Работу выполнил — Погосов Даниэль (https://github.com/DarklleS) при содействии Алексея Федорова (https://github.com/Alexis777777)
*/
template<typename fp_t>
unsigned int fqs(fp_t n, fp_t a, fp_t b, fp_t c, fp_t d, vector<fp_t>& roots);

/*
    Имплементация метода решения уравнения четвертой степени — Mansfield Merriman
    Информация о методе — https://www.jstor.org/stable/2369666
    Работу выполнил — Погосов Даниэль (https://github.com/DarklleS) при содействии Артёма Успенского (https://github.com/MadBunny0)
*/
template<typename fp_t>
unsigned int merriman(fp_t n, fp_t a_, fp_t b_, fp_t c_, fp_t d_, vector<fp_t>& roots);

// Функция поиска корней в случае когда (p^2 - a0) = 0, т. е. p = p1 = sqrt(a0)
// Коэффициенты q, q1 могут быть как действительными, так и комплексносопряженными
// Функция получает на вход коэффициенты a3 и a2, коэффициент p, ссылку на вектор корней,
// возвращает число найденных действительных корней
template<typename fp_t>
unsigned int findRoots(fp_t a3, fp_t a2, fp_t p, vector<fp_t>& roots);

/*
    Имплементация метода решения уравнения четвёртой степени — Solution of quartic equations by William Squire
    Информация о методе — https://www.tandfonline.com/doi/abs/10.1080/0020739790100223
    Работу выполнил — Погосов Даниэль (https://github.com/DarklleS) при содействии Азизхона Ишаннова (https://github.com/ZeekFeed07)
*/
template <typename fp_t>
unsigned int squire(fp_t a, fp_t b, fp_t c, fp_t d, fp_t e, std::vector<fp_t>& roots);

/*
    Имплементация метода решения уравнения четвертой степени — Ungar Method
    Информация о методе — https://core.ac.uk/download/pdf/82351384.pdf
    Работу выполнил — Погосов Даниэль (https://github.com/DarklleS) при содействии Погоса Анесяна (https://github.com/Pogosito)
*/
template<typename fp_t>
unsigned int ungar(fp_t n, fp_t a, fp_t b, fp_t c, fp_t d, vector<fp_t>& roots);

/*
    Имплементация метода решения уравнения четвертой степени —A_Note_on_the_Solution_of_Quartic_Equations-Salzer-1960
    Информация о методе — https://www.ams.org/journals/mcom/1960-14-071/S0025-5718-1960-0117882-6/S0025-5718-1960-0117882-6.pdf
    Работу выполнил — Погосов Даниэль (https://github.com/DarklleS)
*/
template<typename fp_t>
unsigned int salzer(fp_t N, fp_t A, fp_t B, fp_t C, fp_t D, vector<fp_t>& roots);

#endif