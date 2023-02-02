#define _USE_MATH_DEFINES 
#define PRINT true

#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <complex>
#include <numbers>
#include <algorithm>
#include "excerpt.h"

typedef float fp_t;

using namespace std;

// ================================= ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ ================================= //

template<typename fp_t>
inline bool isZero(const fp_t& x)
{
    return FP_ZERO == fpclassify(x);
}

template <typename fp_t>
inline fp_t nearZero(fp_t x)
{
    return x = abs(x) < numeric_limits<fp_t>::epsilon() ? static_cast<fp_t>(0.0L) : x;
}

template<typename fp_t>
inline fp_t fms(fp_t a, fp_t b, fp_t c, fp_t d) // a * b - c * c
{
    fp_t cd = -c * d;

    return fma(a, b, cd) - fma(c, d, cd);
}

template <typename fp_t>
inline int sgn(fp_t x)
{
    return (fp_t(0) < x) - (x < fp_t(0));
}

template<typename fp_t>
complex<fp_t> epsilonC(complex<fp_t> x)
{
    return abs(x) * numeric_limits<fp_t>::epsilon() > abs(x.imag()) ? complex<fp_t>(x.real(), 0) : x;
}

// ======================== ФУНКЦИИ ПОИСКА КОРНЕЙ СТЕПЕННЫХ ПОЛИНОМОВ ======================== //

/*
    Имплементация метода решения кубического уравнения — Numerical Recipes / Viète
    Информация о методе — https://quarticequations.com/Tutorial.pdf (стр. 1)
    Работу выполнил — Погосов Даниэль (https://github.com/DarklleS)
*/
template<typename fp_t>
unsigned int solveCubic(fp_t n, fp_t a, fp_t b, fp_t c, vector<fp_t>& roots)
{
    // Нормировка коэффициентов
    if (isZero(n) || isinf(a /= n))
        return 0;
    if (isinf(b /= n))
        return 0;
    if (isinf(c /= n))
        return 0;

    // Объявление констант
    const fp_t PI = static_cast<fp_t>(M_PI);
    const fp_t ONE_HALF = static_cast<fp_t>(0.5L);
    const fp_t ONE_THIRD = static_cast<fp_t>(1.0L / 3.0L);
    const fp_t ONE_NINTH = static_cast<fp_t>(1.0L / 9.0L);
    const fp_t ONE_54TH = static_cast<fp_t>(1.0L / 54.0L);
    const fp_t A_THIRD = -a * ONE_THIRD;
    constexpr fp_t EPS = numeric_limits<fp_t>::epsilon();

    // Основные расчетные коэффициенты
    fp_t Q = fms(static_cast<fp_t>(3.0L), b, a, a) * ONE_NINTH;
    fp_t R = fms(fms(static_cast<fp_t>(9.0L), b, static_cast<fp_t>(2.0L) * a, a), a, static_cast<fp_t>(27.0L), c) * ONE_54TH;
    
    // Дискриминант
    fp_t D = fms(R, R, -Q, Q * Q);
    
    if (D <= 0) // Обрабатываем вещественные корни
    {
        /* Костыль, который помогает увеличить количество удачных вычислений и случаи когда theta = nand (в основном для NBS метода) */
        if (isZero(D) && abs(Q) <= EPS * EPS && abs(R) <= EPS * EPS * EPS)
        {
            Q = static_cast<fp_t>(0.0L);
            R = static_cast<fp_t>(0.0L);
        }

        if (isZero(Q)) // Все все корни вещественные и кратны трем
        {
            roots =
            {
                A_THIRD,
                A_THIRD,
                A_THIRD
            };
        }
        else // Все все корни вещественные, один из которых кратен двум
        {
            // Дополнительные расчетные коэффициенты
            fp_t theta = acos(R / pow(-Q, static_cast<fp_t>(1.5L)));
            fp_t phi = theta * ONE_THIRD;
            fp_t sqrtQ = static_cast<fp_t>(2.0L) * sqrt(-Q);

            roots =
            {
                fms(sqrtQ, cos(phi), a, ONE_THIRD),
                fms(sqrtQ, cos(fma(static_cast<fp_t>(-2.0L) * PI, ONE_THIRD, phi)), a, ONE_THIRD),
                fms(sqrtQ, cos(fma(static_cast<fp_t>(2.0L) * PI, ONE_THIRD, phi)), a, ONE_THIRD)
            };

            sort(roots.begin(), roots.end());
        }

        return 3;
    }
    else // Обрабатываем два комплексных корня, определяя их вещественность
    {
        // Дополнительные расчетные коэффициенты
        fp_t A = cbrt(abs(R) + sqrt(D));
        fp_t T = R >= 0 || abs(R) <= EPS ? A - Q / A :
                                           Q / A - A;
        // Комплексные корни представлены как: z = x +- iy, где мнимая часть: y = sqrt(3) / 2 * (A + Q / A)
        fp_t imagPart = fms(static_cast<fp_t>(sqrt(3.0L) * 0.5L), A, -static_cast<fp_t>(sqrt(3.0L) * 0.5L), Q / A);

        if (abs(imagPart) <= EPS) // Все корни вещественные, один из которых кратен двум
        {   
            roots =
            {
                fma(-ONE_THIRD, a, T),
                fms(-ONE_THIRD, a, ONE_HALF, T),
                fms(-ONE_THIRD, a, ONE_HALF, T)
            };

            sort(roots.begin(), roots.end());

            return 3;
        }
        else // Только один корень вещественный, остальные комплексные
        {
            roots =
            {
                 fma(-ONE_THIRD, a, T)
            };

            return 1;
        }
    }
}

/*
    Имплементация метода решения кубического уравнения — Numerical Recipes / Viète (версия с комплексной перегрузкой)
    Информация о методе — https://quarticequations.com/Tutorial.pdf (стр. 1)
    Работу выполнил — Погосов Даниэль (https://github.com/DarklleS)
*/
template<typename fp_t>
unsigned int solveCubic(fp_t n, fp_t a, fp_t b, fp_t c, vector<complex<fp_t>>& roots)
{
    // Нормировка коэффициентов
    if (isZero(n) || isinf(a /= n))
        return 0;
    if (isinf(b /= n))
        return 0;
    if (isinf(c /= n))
        return 0;

    // Объявление констант
    const fp_t PI = static_cast<fp_t>(M_PI);
    const fp_t ONE_HALF = static_cast<fp_t>(0.5L);
    const fp_t ONE_THIRD = static_cast<fp_t>(1.0L / 3.0L);
    const fp_t ONE_NINTH = static_cast<fp_t>(1.0L / 9.0L);
    const fp_t ONE_54TH = static_cast<fp_t>(1.0L / 54.0L);
    const fp_t A_THIRD = -a * ONE_THIRD;
    constexpr fp_t EPS = numeric_limits<fp_t>::epsilon();

    // Основные расчетные коэффициенты
    fp_t Q = fms(static_cast<fp_t>(3.0L), b, a, a) * ONE_NINTH;
    fp_t R = fms(fms(static_cast<fp_t>(9.0L), b, static_cast<fp_t>(2.0L) * a, a), a, static_cast<fp_t>(27.0L), c) * ONE_54TH;

    // Дискриминант
    fp_t D = fms(R, R, -Q, Q * Q);

    if (D <= 0) // Обрабатываем вещественные корни
    {
        if (isZero(Q)) // Все все корни вещественные и кратны трем
        {
            roots =
            {
                complex<fp_t>(A_THIRD, static_cast<fp_t>(0.0L)),
                complex<fp_t>(A_THIRD, static_cast<fp_t>(0.0L)),
                complex<fp_t>(A_THIRD, static_cast<fp_t>(0.0L))
            };
        }
        else // Все все корни вещественные, один из которых кратен двум
        {
            // Дополнительные расчетные коэффициенты
            fp_t theta = acos(R / pow(-Q, static_cast<fp_t>(1.5L)));
            fp_t phi = theta * ONE_THIRD;
            fp_t sqrtQ = static_cast<fp_t>(2.0L) * sqrt(-Q);

            roots =
            {
                complex<fp_t>(fms(sqrtQ, cos(phi), a, ONE_THIRD), static_cast<fp_t>(0.0L)),
                complex<fp_t>(fms(sqrtQ, cos(fma(static_cast<fp_t>(-2.0L) * PI, ONE_THIRD, phi)), a, ONE_THIRD), static_cast<fp_t>(0.0L)),
                complex<fp_t>(fms(sqrtQ, cos(fma(static_cast<fp_t>(2.0L) * PI, ONE_THIRD, phi)), a, ONE_THIRD), static_cast<fp_t>(0.0L))
            };

            sort(roots.begin(), roots.end(), 
                [](complex<fp_t>& a, complex<fp_t>& b) 
                {
                    return a.real() < b.real(); 
                });
        }

        return 3;
    }
    else // Обрабатываем два комплексных корня, определяя их вещественность
    {
        // Дополнительные расчетные коэффициенты
        fp_t A = cbrt(abs(R) + sqrt(D));

        fp_t T = R >= 0 || abs(R) <= EPS ? A - Q / A :
                                           Q / A - A;

        // Комплексные корни представлены как: z = x +- iy, где мнимая часть: y = sqrt(3) / 2 * (A + Q / A)
        fp_t imagPart = fms(static_cast<fp_t>(sqrt(3.0L) / 2.0L), A, -static_cast<fp_t>(sqrt(3.0L) / 2.0L), Q / A);

        if (abs(imagPart) <= EPS) // Все корни вещественные, один из которых кратен двум
        {
            roots =
            {
                complex<fp_t>(fma(-ONE_THIRD, a, T), static_cast<fp_t>(0.0L)),
                complex<fp_t>(fms(-ONE_THIRD, a, ONE_HALF, T), static_cast<fp_t>(0.0L)),
                complex<fp_t>(fms(-ONE_THIRD, a, ONE_HALF, T), static_cast<fp_t>(0.0L))
            };

            sort(roots.begin(), roots.end(),
                [](complex<fp_t>& a, complex<fp_t>& b)
                {
                    return a.real() < b.real();
                });

            return 3;
        }
        else // Только один корень вещественный, остальные комплексные
        {
            roots =
            {
                complex<fp_t>(fma(-ONE_THIRD, a, T), static_cast<fp_t>(0.0L)),
                complex<fp_t>(fms(-ONE_THIRD, a, ONE_HALF, T), imagPart),
                complex<fp_t>(fms(-ONE_THIRD, a, ONE_HALF, T), -imagPart)
            };

            return 1;
        }
    }
}

/*
    Имплементация метода решения уравнения четвертой степени — Ferrari's Method
    Информация о методе — https://quarticequations.com/Quartic2.pdf (стр. 3)
    Работу выполнил — Погосов Даниэль (https://github.com/DarklleS)
*/
template<typename fp_t>
unsigned int ferrari(fp_t n, fp_t a, fp_t b, fp_t c, fp_t d, vector<fp_t>& roots)
{
    // Нормировка коэффициентов
    if (isZero(n) || isinf(a /= n))
        return solveCubic(a, b, c, d, roots);
    if (isinf(b /= n))
        return 0;
    if (isinf(c /= n))
        return 0;
    if (isinf(d /= n))
        return 0;

    // Объявление констант
    const fp_t ONE_HALF = static_cast<fp_t>(0.5L);
    const fp_t ONE_QUARTER = static_cast<fp_t>(0.25L);
    constexpr fp_t EPS = numeric_limits<fp_t>::epsilon();

    // Количество вещественных корней
    unsigned numberOfRoots = 0;

    // Вычисляем расчетные коэффициенты
    fp_t C = a * ONE_QUARTER;
    fp_t a_ = fma(static_cast<fp_t>(-6.0L) * C, C, b);
    fp_t b_ = fma(fms(static_cast<fp_t>(8.0L) * C, C, static_cast<fp_t>(2.0L), b), C, c);
    fp_t c_ = fma(fma(fma(static_cast<fp_t>(-3.0L) * C, C, b), C, -c), C, d);

    // Решаем резольвентное кубическое уравнение вида: m^3 + a_ * m^2 + (a_^2 / 4 - c_) * m - b_^2 / 8 = 0
    vector<fp_t> cubicRoots(3);
    unsigned numberOfCubicRoots = solveCubic(static_cast<fp_t>(1.0L), a_, fma(ONE_QUARTER * a_, a_, -c_), -b_ * b_ * static_cast<fp_t>(0.125L), cubicRoots);
    
    // Выбираем корень, который удовлетворяет условию: m > 0 и является вещественным корнем, иначе m = 0
    fp_t m = cubicRoots[numberOfCubicRoots - 1] > 0 ? cubicRoots[numberOfCubicRoots - 1] : static_cast<fp_t>(0.0L);

    // Определяем знак радикала R
    fp_t sigma = b_ > 0 ? static_cast<fp_t>(1.0L) : static_cast<fp_t>(-1.0L);

    fp_t subR = fma(m, m, fma(fma(ONE_QUARTER, a_, m), a_, -c_));
    // Если отрицательно, то имеем комплексный радикал, отсюда уравнение будет иметь исключительно комплексное решение
    if (subR < 0)
        return 0;

    // Радикал R
    fp_t R = sigma * sqrt(subR);

    // Вычисляем радиканды будущего решения
    fp_t radicand = fms(-ONE_HALF, m, ONE_HALF, a_) - R;
    fp_t radicand_ = fms(-ONE_HALF, m, ONE_HALF, a_) + R;

    fp_t sqrtM = sqrt(m * ONE_HALF);

    // Если полученный радиканд (radicand) > 0 либо близок к 0, то уравнение имеет два либо больше вещественных решений
    if (radicand >= 0)
    {
        roots =
        {
            sqrtM - C + sqrt(radicand),
            sqrtM - C - sqrt(radicand)
        };

        numberOfRoots += 2;
    }
    else if (abs(radicand) <= EPS)
    {
        roots =
        {
            sqrtM - C,
            sqrtM - C
        };

        numberOfRoots += 2;
    }

    // Если полученный радиканд (radicand_) > 0 либо близок к 0, то уравнение имеет два либо больше вещественных решений
    if (radicand_ >= 0)
    {
        roots[numberOfRoots] = -sqrtM - C + sqrt(radicand_);
        roots[numberOfRoots + 1] = -sqrtM - C - sqrt(radicand_);

        numberOfRoots += 2;
    }
    else if (abs(radicand_) <= EPS)
    {
        roots[numberOfRoots] = -sqrtM - C;
        roots[numberOfRoots + 1] = -sqrtM - C;

        numberOfRoots += 2;
    }

    return numberOfRoots;
}

/*
    Имплементация метода решения уравнения четвертой степени — Descartes’ Method 
    Информация о методе — https://quarticequations.com/Quartic2.pdf (стр. 4)
    Работу выполнил — Погосов Даниэль (https://github.com/DarklleS)
*/
template<typename fp_t>
unsigned int descartes(fp_t n, fp_t a, fp_t b, fp_t c, fp_t d, vector<fp_t>& roots)
{
    // Нормировка коэффициентов
    if (isZero(n) || isinf(a /= n))
        return solveCubic(a, b, c, d, roots);
    if (isinf(b /= n))
        return 0;
    if (isinf(c /= n))
        return 0;
    if (isinf(d /= n))
        return 0;

    // Объявление констант
    const fp_t ONE_HALF = static_cast<fp_t>(0.5L);
    const fp_t ONE_QUARTER = static_cast<fp_t>(0.25L);
    constexpr fp_t EPS = numeric_limits<fp_t>::epsilon();

    // Количество вещественных корней
    unsigned numberOfRoots = 0;

    // Вычисляем расчетные коэффициенты
    fp_t C = a * ONE_QUARTER;
    fp_t a_ = fma(static_cast<fp_t>(-6.0L) * C, C, b);
    fp_t b_ = fma(fms(static_cast<fp_t>(8.0L) * C, C, static_cast<fp_t>(2.0L), b), C, c);
    fp_t c_ = fma(fma(fma(static_cast<fp_t>(-3.0L) * C, C, b), C, -c), C, d);

    // Решаем резольвентное бикубическое уравнение вида: y^6 + 2 * a_ * y^4 + (a_^2 - 4 * c_) * y^2 - b_^2 = 0
    vector<fp_t> cubicRoots(3);
    unsigned numberOfCubicRoots = solveCubic(static_cast<fp_t>(1.0L), static_cast<fp_t>(2.0L) * a_, fms(a_, a_, static_cast<fp_t>(4.0L), c_), -b_ * b_, cubicRoots);

    // Выбираем корень, который удовлетворяет условию: y^2 > 0 и является вещественным корнем, иначе y^2 = 0
    fp_t yy = cubicRoots[numberOfCubicRoots - 1] > 0 ? cubicRoots[numberOfCubicRoots - 1] : static_cast<fp_t>(0.0L);

    // Положительный вещественный корень бикубического уравнения
    fp_t y = sqrt(yy);

    // Определяем знак радикала R
    fp_t sigma = b_ > 0 ? static_cast<fp_t>(1.0L) : static_cast<fp_t>(-1.0L);

    fp_t subR = fma(ONE_QUARTER * a_, a_, fma(fms(ONE_QUARTER, yy, -ONE_HALF, a_), yy, -c_));
    // Если отрицательно, то имеем комплексный радикал, отсюда уравнение будет иметь исключительно комплексное решение
    if (subR < 0)
        return 0;

    // Радикал R
    fp_t R = sigma * sqrt(subR);

    // Вычисляем радиканды будущего решения
    fp_t radicand = fms(-ONE_QUARTER, yy, ONE_HALF, a_) - R;
    fp_t radicand_ = fms(-ONE_QUARTER, yy, ONE_HALF, a_) + R;

    // Если полученный радиканд (radicand) > 0 либо близок к 0, то уравнение имеет два либо больше вещественных решений
    if (radicand >= 0)
    {
        roots =
        {
            fma(ONE_HALF, y, sqrt(radicand)) - C,
            fma(ONE_HALF, y, -sqrt(radicand)) - C
        };

        numberOfRoots += 2;
    }
    else if (abs(radicand) <= EPS)
    {
        roots =
        {
            fma(ONE_HALF, y, -C),
            fma(ONE_HALF, y, -C)
        };

        numberOfRoots += 2;
    }

    // Если полученный радиканд (radicand) > 0 либо близок к 0, то уравнение имеет два либо больше вещественных решений
    if (radicand_ >= 0)
    {
        roots[numberOfRoots] = fma(-ONE_HALF, y, sqrt(radicand_)) - C;
        roots[numberOfRoots + 1] = fma(-ONE_HALF, y, -sqrt(radicand_)) - C;

        numberOfRoots += 2;
    }
    else if (abs(radicand_) <= EPS)
    {
        roots[numberOfRoots] = fma(-ONE_HALF, y, -C);
        roots[numberOfRoots + 1] = fma(-ONE_HALF, y, -C);

        numberOfRoots += 2;
    }

    return numberOfRoots;
}

/*
    Имплементация метода решения уравнения четвертой степени — NBS Method
    Информация о методе — https://quarticequations.com/Quartic2.pdf (стр. 5)
    Работу выполнил — Погосов Даниэль (https://github.com/DarklleS)
*/
template<typename fp_t>
unsigned int nbs(fp_t n, fp_t a, fp_t b, fp_t c, fp_t d, vector<fp_t>& roots)
{
    // Нормировка коэффициентов
    if (isZero(n) || isinf(a /= n))
        return solveCubic(a, b, c, d, roots);
    if (isinf(b /= n))
        return 0;
    if (isinf(c /= n))
        return 0;
    if (isinf(d /= n))
        return 0;

    // Объявление констант
    const fp_t ONE_HALF = static_cast<fp_t>(0.5L);
    const fp_t ONE_QUARTER = static_cast<fp_t>(0.25L);
    constexpr fp_t EPS = numeric_limits<fp_t>::epsilon();

    // Количество вещественных корней
    unsigned numberOfRoots = 0;

    // Вычисляем расчетные коэффициенты
    fp_t a_ = -b;
    fp_t b_ = fms(a, c, static_cast<fp_t>(4.0L), d);
    fp_t c_ = fms(fms(static_cast<fp_t>(4.0L), b, a, a), d, c, c);

    // Решаем резольвентное кубическое уравнение вида: u^3 + a_ * u^2 + b_ * u + c_ = 0
    vector<fp_t> cubicRoots(3);
    unsigned numberOfCubicRoots = solveCubic(static_cast<fp_t>(1.0L), a_, b_, c_, cubicRoots);

    // Выбираем наибольший вещественный корень
    fp_t u = cubicRoots[numberOfCubicRoots - 1];

    // Определяем знак радикандов
    fp_t sigma = fma(-ONE_HALF * a, u, c) > 0 ? static_cast<fp_t>(1.0L) : static_cast<fp_t>(-1.0L);

    // Вычисляем рандиканды переменных p и q
    fp_t radicandP = fma(ONE_QUARTER * a, a, u) - b;
    fp_t radicandQ = fma(ONE_QUARTER * u, u, -d);

    // Если хоть один из радикандов отрицателен, то имеем исключительно комплексное решение
    if (radicandP < 0 || radicandQ < 0)
        return 0;

    // Вычисляем итоговые расчетные коэффициенты p_n и q_n (n = 1, 2)
    fp_t p1 = fma(ONE_HALF, a, -sqrt(radicandP));
    fp_t p2 = fma(ONE_HALF, a, sqrt(radicandP));
    fp_t q1 = fms(ONE_HALF, u, -sigma, sqrt(radicandQ));
    fp_t q2 = fms(ONE_HALF, u, sigma, sqrt(radicandQ));

    // Вычисляем радиканды будущего решения
    fp_t radicand = fma(ONE_QUARTER * p1, p1, -q1);
    fp_t radicand_ = fma(ONE_QUARTER * p2, p2, -q2);
    
    // Если полученный радиканд (radicand) > 0 либо близок к 0, то уравнение имеет два либо больше вещественных решений
    if (radicand >= 0)
    {
        roots =
        {
            fma(-ONE_HALF, p1, sqrt(radicand)),
            fma(-ONE_HALF, p1, -sqrt(radicand))
        };

        numberOfRoots += 2;
    }
    else if (abs(radicand) <= EPS)
    {
        roots =
        {
            -p1 * ONE_HALF,
            -p1 * ONE_HALF
        };

        numberOfRoots += 2;
    }

    // Если полученный радиканд (radicand_) > 0 либо близок к 0, то уравнение имеет два либо больше вещественных решений
    if (radicand_ >= 0)
    {

        roots[numberOfRoots] = fma(-ONE_HALF, p2, sqrt(radicand_));
        roots[numberOfRoots + 1] = fma(-ONE_HALF, p2, -sqrt(radicand_));;

        numberOfRoots += 2;
    }
    else if (abs(radicand_) <= EPS)
    {
        roots[numberOfRoots] = -p2 * ONE_HALF;
        roots[numberOfRoots + 1] = -p2 * ONE_HALF;

        numberOfRoots += 2;
    }

    
    return numberOfRoots;
}

/*
    Имплементация метода решения уравнения четвертой степени — Euler’s Method 
    Информация о методе — https://quarticequations.com/Quartic2.pdf (стр. 6)
    Работу выполнил — Погосов Даниэль (https://github.com/DarklleS)
*/
template<typename fp_t>
unsigned int euler(fp_t n, fp_t a, fp_t b, fp_t c, fp_t d, vector<fp_t>& roots)
{
    // Нормировка коэффициентов
    if (isZero(n) || isinf(a /= n))
        return solveCubic(a, b, c, d, roots);
    if (isinf(b /= n))
        return 0;
    if (isinf(c /= n))
        return 0;
    if (isinf(d /= n))
        return 0;

    // Объявление констант
    const fp_t ONE_HALF = static_cast<fp_t>(0.5L);
    const fp_t ONE_QUARTER = static_cast<fp_t>(0.25L);
    constexpr fp_t EPS = numeric_limits<fp_t>::epsilon();

    // Количество вещественных корней
    unsigned numberOfRoots = 0;

    // Вычисляем расчетные коэффициенты
    fp_t C = a * ONE_QUARTER;
    fp_t a_ = fma(static_cast<fp_t>(-6.0L) * C, C, b);
    fp_t b_ = fma(fms(static_cast<fp_t>(8.0L) * C, C, static_cast<fp_t>(2.0L), b), C, c);
    fp_t c_ = fma(fma(fma(static_cast<fp_t>(-3.0L) * C, C, b), C, -c), C, d);

    // Решаем резольвентное кубическое уравнение вида: r^3 + a_ / 2 * r^2 + (a_^2 - 4 * c_) / 16 * r - b_^2 / 64 = 0
    vector<complex<fp_t>> cubicRoots(3);
    unsigned numberOfCubicRoots = solveCubic(static_cast<fp_t>(1.0L), a_ * ONE_HALF, fms(static_cast<fp_t>(0.0625L) * a_, a_, ONE_QUARTER, c_), -b_ * b_ * static_cast<fp_t>(0.015625L), cubicRoots);

    // Определяем знак радикандов
    fp_t sigma = b_ > 0 ? static_cast<fp_t>(1.0L) : static_cast<fp_t>(-1.0L);

    // Объявление вычислительный переменных для поиска радикандов
    fp_t r;
    fp_t x2, x3, y;
    fp_t subradicand;
    fp_t radicand, radicand_;

    // Вычисление радикандов. Разбераем 2 случая дабы оптимизировать вычисления
    // - Если решение резольвентного кубического уравнения имеет единственный вещественный корень
    if (numberOfCubicRoots == 1)
    {
        r = cubicRoots[0].real();
        if (r < 0)
            return 0;

        x2 = cubicRoots[1].real();
        x3 = x2;
        y = cubicRoots[1].imag();

        subradicand = fms(x2, x3, -y, y);
        // Если отрицательно, то имеем комплексные радиканды, отсюда уравнение будет иметь исключительно комплексное решение
        if (subradicand < 0)
            return 0;

        radicand = fms(static_cast<fp_t>(2.0L), x2, static_cast<fp_t>(2.0L) * sigma, sqrt(subradicand));
        radicand_ = fms(static_cast<fp_t>(2.0L), x2, static_cast<fp_t>(-2.0L) * sigma, sqrt(subradicand));
    }
    // - Иначе все корни резольвентного кубического уравнения вещественные
    else
    {
        r = cubicRoots[2].real();
        if (r < 0)
            return 0;

        x2 = cubicRoots[0].real();
        x3 = cubicRoots[1].real();
        y = static_cast<fp_t>(0.0L);

        subradicand = x2 * x3;
        // Если отрицательно, то имеем комплексные радиканды, отсюда уравнение будет иметь исключительно комплексное решение
        if (subradicand < 0)
            return 0;

        radicand = fma(static_cast<fp_t>(-2.0L) * sigma, sqrt(subradicand), x2) + x3;
        radicand_ = fma(static_cast<fp_t>(2.0L) * sigma, sqrt(subradicand), x2) + x3;
    }

    fp_t sqrtR = sqrt(r);

    // Если полученный радиканд (radicand) > 0 либо близок к 0, то уравнение имеет два либо больше вещественных решений
    if (radicand >= 0)
    {
        roots =
        {
            sqrtR + sqrt(radicand) - C,
            sqrtR - sqrt(radicand) - C
        };

        numberOfRoots += 2;
    }
    else if (abs(radicand) <= EPS)
    {
        roots =
        {
            sqrtR - C,
            sqrtR - C
        };

        numberOfRoots += 2;
    }

    // Если полученный радиканд (radicand_) > 0 либо близок к 0, то уравнение имеет два либо больше вещественных решений
    if (radicand_ >= 0)
    {

        roots[numberOfRoots] = -sqrtR + sqrt(radicand_) - C;
        roots[numberOfRoots + 1] = -sqrtR - sqrt(radicand_) - C;

        numberOfRoots += 2;
    }
    else if (abs(radicand_) <= EPS)
    {
        roots[numberOfRoots] = -sqrtR - C;
        roots[numberOfRoots + 1] = -sqrtR - C;

        numberOfRoots += 2;
    }

    return numberOfRoots;
}

/*
    Имплементация метода решения уравнения четвертой степени — Van der Waerden’s Method 
    Информация о методе — https://quarticequations.com/Quartic2.pdf (стр. 7)
    Работу выполнил — Погосов Даниэль (https://github.com/DarklleS)
*/
template<typename fp_t>
unsigned int vanDerWaerden(fp_t n, fp_t a, fp_t b, fp_t c, fp_t d, vector<fp_t>& roots)
{
    // Нормировка коэффициентов
    if (isZero(n) || isinf(a /= n))
        return solveCubic(a, b, c, d, roots);
    if (isinf(b /= n))
        return 0;
    if (isinf(c /= n))
        return 0;
    if (isinf(d /= n))
        return 0;

    // Объявление констант
    const fp_t ONE_HALF = static_cast<fp_t>(0.5L);
    const fp_t ONE_QUARTER = static_cast<fp_t>(0.25L);
    constexpr fp_t EPS = numeric_limits<fp_t>::epsilon();

    // Количество вещественных корней
    unsigned numberOfRoots = 0;

    // Вычисляем расчетные коэффициенты
    fp_t C = a * ONE_QUARTER;
    fp_t a_ = fma(static_cast<fp_t>(-6.0L) * C, C, b);
    fp_t b_ = fma(fms(static_cast<fp_t>(8.0L) * C, C, static_cast<fp_t>(2.0L), b), C, c);
    fp_t c_ = fma(fma(fma(static_cast<fp_t>(-3.0L) * C, C, b), C, -c), C, d);

    // Решаем резольвентное кубическое уравнение вида: theta^3 - 2 * a_ * theta^2 + (a_^2 - 4 * c_) * theta + b_^2 = 0
    vector<complex<fp_t>> cubicRoots(3);
    unsigned numberOfCubicRoots = solveCubic(static_cast<fp_t>(1.0L), static_cast<fp_t>(-2.0L) * a_, fms(a_, a_, static_cast<fp_t>(4.0L), c_), b_ * b_, cubicRoots);

    // Определяем знак радикандов
    fp_t sigma = b_ > 0 ? static_cast<fp_t>(1.0L) : static_cast<fp_t>(-1.0L);

    // Объявление вычислительный переменных для поиска радикандов
    fp_t theta;
    fp_t x2, x3, y;
    fp_t subradicand;
    fp_t radicand, radicand_;

    // Вычисление радикандов. Разбераем 2 случая дабы оптимизировать вычисления
    // - Если решение резольвентного кубического уравнения имеет единственный вещественный корень
    if (numberOfCubicRoots == 1)
    {
        theta = cubicRoots[0].real();
        if (theta > 0)
            return 0;

        x2 = cubicRoots[1].real();
        x3 = x2;
        y = cubicRoots[1].imag();

        subradicand = fms(x2, x3, -y, y);
        // Если отрицательно, то имеем комплексные радиканды, отсюда уравнение будет иметь исключительно комплексное решение
        if (subradicand < 0)
            return 0;

        radicand = fms(static_cast<fp_t>(-2.0L), x2, static_cast<fp_t>(2.0L) * sigma, sqrt(subradicand));
        radicand_ = fms(static_cast<fp_t>(-2.0L), x2, static_cast<fp_t>(-2.0L) * sigma, sqrt(subradicand));
    }
    // - Иначе все корни резольвентного кубического уравнения вещественные
    else
    {
        theta = cubicRoots[0].real();
        if (theta > 0)
            return 0;

        x2 = cubicRoots[1].real();
        x3 = cubicRoots[2].real();
        y = static_cast<fp_t>(0.0L);

        subradicand = x2 * x3;
        // Если отрицательно, то имеем комплексные радиканды, отсюда уравнение будет иметь исключительно комплексное решение
        if (subradicand < 0)
            return 0;

        radicand = fma(static_cast<fp_t>(-2.0L) * sigma, sqrt(subradicand), -x2) - x3;
        radicand_ = fma(static_cast<fp_t>(2.0L) * sigma, sqrt(subradicand), -x2) - x3;
    }

    fp_t sqrtTheta = sqrt(-theta);

    // Если полученный радиканд (radicand) > 0 либо близок к 0, то уравнение имеет два либо больше вещественных решений
    if (radicand >= 0)
    {
        roots =
        {
            fma(ONE_HALF, sqrt(radicand), fma(ONE_HALF, sqrtTheta, -C)),
            fma(-ONE_HALF, sqrt(radicand), fma(ONE_HALF, sqrtTheta, -C))
        };

        numberOfRoots += 2;
    }
    else if (abs(radicand) <= EPS)
    {
        roots =
        {
            fma(ONE_HALF, sqrtTheta, -C),
            fma(ONE_HALF, sqrtTheta, -C)
        };

        numberOfRoots += 2;
    }

    // Если полученный радиканд (radicand_) > 0 либо близок к 0, то уравнение имеет два либо больше вещественных решений
    if (radicand_ >= 0)
    {

        roots[numberOfRoots] = fma(ONE_HALF, sqrt(radicand_), fma(-ONE_HALF, sqrtTheta, -C));
        roots[numberOfRoots + 1] = fma(-ONE_HALF, sqrt(radicand_), fma(-ONE_HALF, sqrtTheta, -C));

        numberOfRoots += 2;
    }
    else if (abs(radicand_) <= EPS)
    {
        roots[numberOfRoots] = fma(-ONE_HALF, sqrtTheta, -C);
        roots[numberOfRoots + 1] = fma(-ONE_HALF, sqrtTheta, -C);

        numberOfRoots += 2;
    }

    return numberOfRoots;
}

// ======================== ФУНКЦИИ ТЕСТИРОВАНИЯ МЕТОДОВ ======================== //

template <typename fp_t>
void testCubicAdv(const int testCount, const fp_t dist)
{
    int P = 3; // power, total number of tests
    fp_t low = -1, high = 1; // [low, high], max distance between clustered roots
    fp_t absMaxError, relMaxError; // variables for each test Errors
    int numOfFoundRoots, cantFind = 0;
    fp_t maxAbsAllofTest = -1, maxRelAllofTest = -1; // maximum from maxAbsoluteError and maxRelError from all [testCount] tests

    long double absErrors = 0;
    long double relError = 0;
    int count = 0;
    fp_t count1 = 0;

    vector<fp_t> coefficients(P + 1);
    vector<fp_t> trueRoots(P);

    for (size_t i = 0; i < testCount; ++i)
    {
        vector<fp_t> foundRoots(P);

        generate_polynomial<fp_t>(P, 0, P, 0, dist, low, high, trueRoots, coefficients);
        numOfFoundRoots = solveCubic<fp_t>(coefficients[3], coefficients[2], coefficients[1], coefficients[0], foundRoots);

        compare_roots<fp_t>(numOfFoundRoots, P, foundRoots, trueRoots, absMaxError, relMaxError);

        if (isinf(absMaxError))
            cantFind += 1;
        else
        {
            maxAbsAllofTest = absMaxError > maxAbsAllofTest ? absMaxError : maxAbsAllofTest;
            absErrors += absMaxError;
            maxRelAllofTest = relMaxError > maxRelAllofTest ? relMaxError : maxRelAllofTest;
            relError += relMaxError;

            count += relMaxError > 1 ? 1 : 0;
        }
    }

    if (PRINT)
    {
        cout << "\n\n\t\t\tCUBIC TEST RESULTS\n\n";
        cout << "Max distance: " << dist << endl;
        cout << "Total count of tests: " << testCount << endl;
        cout << "Couldn't find roots: " << cantFind << " times " << endl;
        cout << "Mean absMaxError = " << absErrors / (testCount - cantFind) << endl;
        cout << "Max {absMaxError_i | i = 0, ..., 1e6} from all of the tests: " << maxAbsAllofTest << endl;
        cout << "Mean RelMaxError = " << relError / (testCount - cantFind) << endl;
        cout << "Max {RelMaxError_i | i = 0, ..., 1e6} all of the tests: " << maxRelAllofTest << endl;
        cout << "RelMaxError > 1: " << count << " times" << endl;
    }
}

template<typename fp_t>
void testQuarticAdv(const int testCount, const fp_t dist) {
    int P = 4; // power
    fp_t low = -1, high = 1; // [low, high]
    fp_t absMaxError, relMaxError; // variables for each test Errors
    int numOfFoundRoots, cantFind = 0;
    fp_t maxAbsAllofTest = -1, maxRelAllofTest = -1; // maximum from maxAbsoluteError and maxRelError from all [testCount] tests

    long double absErrors = 0;
    long double relError = 0;
    int count = 0;

    vector<fp_t> coefficients(P + 1);

    for (size_t i = 0; i < testCount; ++i) 
    {
        vector<fp_t> foundRoots(P);
        vector<fp_t> trueRoots(P);
        fp_t rad;
        fp_t rad_;
        generate_polynomial<fp_t>(P, 0, P, 0, dist,
            low, high, trueRoots, coefficients);
        numOfFoundRoots = nbs<fp_t>(coefficients[4], coefficients[3], coefficients[2], coefficients[1], coefficients[0], foundRoots);

        compare_roots<fp_t>(numOfFoundRoots, 4, foundRoots, trueRoots, absMaxError, relMaxError);

        if (isinf(absMaxError))
            cantFind += 1;
        else {
            maxAbsAllofTest = absMaxError > maxAbsAllofTest ? absMaxError : maxAbsAllofTest;
            absErrors += absMaxError;
            maxRelAllofTest = relMaxError > maxRelAllofTest ? relMaxError : maxRelAllofTest;
            relError += relMaxError;

            count += relMaxError > 1 ? 1 : 0;
        }
    }

    if (PRINT)
    {
        cout << "\n\n\t\t\tQUARTIC TEST RESULTS\n\n";
        cout << "Max distance: " << dist << endl;
        cout << "Total count of tests: " << testCount << endl;
        cout << "Couldn't find roots: " << cantFind << " times " << endl;
        cout << "Mean absMaxError = " << absErrors / (testCount - cantFind) << endl;
        cout << "Max {absMaxError_i | i = 0, ..., 1e6} from all of the tests: " << maxAbsAllofTest << endl;
        cout << "Mean RelMaxError = " << relError / (testCount - cantFind) << endl;
        cout << "Max {RelMaxError_i | i = 0, ..., 1e6} all of the tests: " << maxRelAllofTest << endl;
        cout << "RelMaxError > 1: " << count << " times" << endl;
    }
}

// ============================================================================= //

int main()
{
    //testCubicAdv(1000000, static_cast<fp_t>(1e-05L));
    testQuarticAdv(10000000, static_cast<fp_t>(1e-05L));
}