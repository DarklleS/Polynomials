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

// Функция проверки элемента на ноль (x == 0)
template<typename fp_t>
inline bool isZero(static const fp_t& x)
{
    return FP_ZERO == fpclassify(x);
}

// Функция проверки равенства двух переменных (A == B)
template<typename fp_t>
inline bool isEqual(const fp_t& A, const fp_t& B)
{
    return abs(A - B) < numeric_limits<fp_t>::epsilon() * (abs(B) + abs(A)) * static_cast<fp_t>(0.5L);
}

// Функция вычисления разности произведений с оптимизацией погрешности (a * b - c * d)
template<typename fp_t>
inline fp_t fms(fp_t a, fp_t b, fp_t c, fp_t d) 
{
    fp_t cd = -c * d;

    return fma(a, b, cd) - fma(c, d, cd);
}

template<typename fp_t>
complex<fp_t> fmac(complex<fp_t> x, complex<fp_t> y, complex<fp_t> z)
{
    fp_t r = fma(-x.imag(), y.imag(), fma(x.real(), y.real(), z.real()));
    fp_t i = fma(y.real(), x.imag(), fma(x.real(), y.imag(), z.imag()));

    return complex<fp_t>(r, i);
}

template<typename fp_t>
complex<fp_t> fmsc(complex<fp_t> a, complex<fp_t> b, complex<fp_t> c, complex<fp_t> d)
{
    complex<fp_t> cd = -c * d;

    return fmac(a, b, cd) - fmac(c, d, cd);
}

// Сигнатурная функция
template <typename fp_t>
inline int sgn(static const fp_t& x)
{
    return (fp_t(0) < x) - (x < fp_t(0));
}

// Функция определяющая является ли число вещественным, путем определения того, насколько мала мнимая часть
template<typename fp_t>
inline complex<fp_t> epsilonComplex(complex<fp_t> x)
{
    return abs(x) * numeric_limits<fp_t>::epsilon() > abs(x.imag()) ? complex<fp_t>(x.real(), 0) : x;
}

template<typename fp_t>
inline bool isComplex(static const complex<fp_t>& x)
{
    return abs(x) * numeric_limits<fp_t>::epsilon() <= abs(x.imag());
}

template<typename fp_t>
fp_t dscrmt(fp_t A, fp_t B, fp_t C)
{
    fp_t p = B * B;
    fp_t q = A * C;
    //Use the hardware's FMA
    fp_t dp = std::fma(B, B, -p);
    fp_t dq = std::fma(A, C, -q);
    fp_t d = (p - q) + (dp - dq);
    return d;
}

template<typename fp_t>
complex<fp_t> dscrmt(complex<fp_t> A, complex<fp_t> B, complex<fp_t> C) 
{
    complex<fp_t> p = B * B;
    complex<fp_t> q = A * C;
    //Use the hardware's FMA
    complex<fp_t> dp = fmac(B, B, -p);
    complex<fp_t> dq = fmac(A, C, -q);
    complex<fp_t> d = (p - q) + (dp - dq);
    return d;
}

// ======================== ФУНКЦИИ ПОИСКА КОРНЕЙ СТЕПЕННЫХ ПОЛИНОМОВ ======================== //

/*
    Имплементация метода решения линейного уравнения
*/
template<typename fp_t>
unsigned int solveLinear(fp_t a, fp_t b, vector<fp_t>& roots)
{
    roots[0] = -b / a;

    return isnan(roots[0]) || isinf(roots[0]) ? 0 : 1;
}

/*
    Имплементация метода решения квадратного уравнения
*/
template<typename fp_t>
unsigned int solveQuadratic(fp_t n, fp_t a, fp_t b, vector<fp_t>& roots)
{
    if (isZero(n) || isinf(a /= n))
        return 0;
    if (isinf(b /= n))
        return 0;
     
    a /= static_cast<fp_t>(-2.0L);

    fp_t d = dscrmt(n, a, b);

    if (d < 0)
        return 0;

    fp_t S = a;
    S = fma(sqrt(d), (sgn(S) + (S == 0)), S);

    fp_t Z1 = S;
    fp_t Z2 = (b) / S;

    roots[0] = Z1;
    roots[1] = Z2;

    return 2;
}

template<typename fp_t>
unsigned int solveQuadratic(fp_t n, fp_t a, fp_t b, vector<complex<fp_t>>& roots)
{
    if (isZero(n) || isinf(a /= n))
        return 0;
    if (isinf(b /= n))
        return 0;

    a /= static_cast<fp_t>(-2.0L);

    fp_t d = dscrmt(n, a, b);
    complex<fp_t> sqrtD = sqrt(complex<fp_t>(d, static_cast<fp_t>(0.0L)));

    complex<fp_t> S(a, static_cast<fp_t>(0.0L));
    complex<fp_t> partS(sgn(a) + (isZero(a), static_cast<fp_t>(0.0L)));

    S = fmac(sqrtD, partS, S);

    complex<fp_t> Z1 = S;
    complex<fp_t> Z2 = (b) / S;

    roots[0] = Z1;
    roots[1] = Z2;
}

template<typename fp_t>
unsigned int solveQuadratic(complex<fp_t> n, complex<fp_t> a, complex<fp_t> b, vector<complex<fp_t>>& roots)
{
    a /= static_cast<fp_t>(-2.0L);

    complex<fp_t> d = dscrmt(n, a, b);
    complex<fp_t> sqrtD = sqrt(complex<fp_t>(d));

    complex<fp_t> S(a);
    complex<fp_t> partS(sgn(a) + (isZero(a.real()) && isZero(a.imag()), static_cast<fp_t>(0.0L)));

    S = fmac(sqrtD, partS, S);

    complex<fp_t> Z1 = S;
    complex<fp_t> Z2 = (b) / S;

    roots[0] = Z1;
    roots[1] = Z2;

    return 0;
}

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
        return solveQuadratic(a, b, c, roots);
    if (isinf(b /= n))
        return 0;
    if (isinf(c /= n))
        return 0;

    // Объявление констант
    static const fp_t PI = static_cast<fp_t>(numbers::pi);
    static const fp_t ONE_HALF = static_cast<fp_t>(0.5L);
    static const fp_t ONE_THIRD = static_cast<fp_t>(1.0L / 3.0L);
    static const fp_t ONE_NINTH = static_cast<fp_t>(1.0L / 9.0L);
    static const fp_t ONE_54TH = static_cast<fp_t>(1.0L / 54.0L);
    static const fp_t A_THIRD = -a * ONE_THIRD;
    static constexpr fp_t EPS = numeric_limits<fp_t>::epsilon();

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
        fp_t T = R >= -EPS ? A - Q / A :
                             Q / A - A;
        
        fp_t sqrt3Half = static_cast<fp_t>(sqrt(3.0L) * 0.5L);

        // Комплексные корни представлены как: z = x +- iy, где мнимая часть: y = sqrt(3) / 2 * (A + Q / A)
        fp_t imagPart = fms(sqrt3Half, A, -sqrt3Half, Q / A);

        if (abs(imagPart) <= EPS) // Все корни вещественные, один из которых кратен двум
        {
            fp_t multipleRoot = fms(-ONE_THIRD, a, ONE_HALF, T);

            roots =
            {
                fma(-ONE_THIRD, a, T),
                multipleRoot,
                multipleRoot
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
    static const fp_t PI = static_cast<fp_t>(numbers::pi);
    static const fp_t ONE_HALF = static_cast<fp_t>(0.5L);
    static const fp_t ONE_THIRD = static_cast<fp_t>(1.0L / 3.0L);
    static const fp_t ONE_NINTH = static_cast<fp_t>(1.0L / 9.0L);
    static const fp_t ONE_54TH = static_cast<fp_t>(1.0L / 54.0L);
    static const fp_t A_THIRD = -a * ONE_THIRD;
    static constexpr fp_t EPS = numeric_limits<fp_t>::epsilon();

    // Основные расчетные коэффициенты
    fp_t Q = fms(static_cast<fp_t>(3.0L), b, a, a) * ONE_NINTH;
    fp_t R = fms(fms(static_cast<fp_t>(9.0L), b, static_cast<fp_t>(2.0L) * a, a), a, static_cast<fp_t>(27.0L), c) * ONE_54TH;

    // Дискриминант
    fp_t D = fms(R, R, -Q, Q * Q);

    if (D <= 0) // Обрабатываем вещественные корни
    {
        if (isZero(Q)) // Все все корни вещественные и кратны трем
        {
            complex<fp_t> root(A_THIRD, static_cast<fp_t>(0.0L));

            roots =
            {
                root,
                root,
                root
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
        fp_t T = R >= -EPS ? A - Q / A :
                             Q / A - A;

        fp_t sqrt3Half = static_cast<fp_t>(sqrt(3.0L) * 0.5L);

        // Комплексные корни представлены как: z = x +- iy, где мнимая часть: y = sqrt(3) / 2 * (A + Q / A)
        fp_t imagPart = fms(sqrt3Half, A, -sqrt3Half, Q / A);

        if (abs(imagPart) <= EPS) // Все корни вещественные, один из которых кратен двум
        {
            complex<fp_t> multipleRoot(fms(-ONE_THIRD, a, ONE_HALF, T), static_cast<fp_t>(0.0L));

            roots =
            {
                complex<fp_t>(fma(-ONE_THIRD, a, T), static_cast<fp_t>(0.0L)),
                multipleRoot,
                multipleRoot
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
            fp_t rootRealPart = fms(-ONE_THIRD, a, ONE_HALF, T);

            roots =
            {
                complex<fp_t>(fma(-ONE_THIRD, a, T), static_cast<fp_t>(0.0L)),
                complex<fp_t>(rootRealPart, imagPart),
                complex<fp_t>(rootRealPart, -imagPart)
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
    static const fp_t ONE_HALF = static_cast<fp_t>(0.5L);
    static const fp_t ONE_QUARTER = static_cast<fp_t>(0.25L);
    static constexpr fp_t EPS = numeric_limits<fp_t>::epsilon();

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

    fp_t radicandPart = fms(-ONE_HALF, m, ONE_HALF, a_);

    // Вычисляем радиканды будущего решения
    fp_t radicand = radicandPart - R;
    fp_t radicand_ = radicandPart + R;

    fp_t sqrtM = sqrt(m * ONE_HALF);

    fp_t rootPart = sqrtM - C;

    // Если полученный радиканд (radicand) > 0 либо близок к 0, то уравнение имеет два либо больше вещественных решений
    if (radicand >= 0)
    {
        fp_t radical = sqrt(radicand);

        roots =
        {
            rootPart + radical,
            rootPart - radical
        };

        numberOfRoots += 2;
    }
    else if (abs(radicand) <= EPS)
    {
        roots =
        {
            rootPart,
            rootPart
        };

        numberOfRoots += 2;
    }

    rootPart = -sqrtM - C;

    // Если полученный радиканд (radicand_) > 0 либо близок к 0, то уравнение имеет два либо больше вещественных решений
    if (radicand_ >= 0)
    {
        fp_t radical = sqrt(radicand_);

        roots[numberOfRoots] = rootPart + radical;
        roots[numberOfRoots + 1] = rootPart - radical;

        numberOfRoots += 2;
    }
    else if (abs(radicand_) <= EPS)
    {
        roots[numberOfRoots] = rootPart;
        roots[numberOfRoots + 1] = rootPart;

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
    static const fp_t ONE_HALF = static_cast<fp_t>(0.5L);
    static const fp_t ONE_QUARTER = static_cast<fp_t>(0.25L);
    static constexpr fp_t EPS = numeric_limits<fp_t>::epsilon();

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

    fp_t radicandPart = fms(-ONE_QUARTER, yy, ONE_HALF, a_);

    // Вычисляем радиканды будущего решения
    fp_t radicand = radicandPart - R;
    fp_t radicand_ = radicandPart + R;

    fp_t rootPart = fma(ONE_HALF, y, -C);

    // Если полученный радиканд (radicand) > 0 либо близок к 0, то уравнение имеет два либо больше вещественных решений
    if (radicand >= 0)
    {
        fp_t radical = sqrt(radicand);

        roots =
        {
            rootPart + radical,
            rootPart - radical
        };

        numberOfRoots += 2;
    }
    else if (abs(radicand) <= EPS)
    {
        roots =
        {
            rootPart,
            rootPart
        };

        numberOfRoots += 2;
    }

    rootPart = fma(-ONE_HALF, y, -C);

    // Если полученный радиканд (radicand) > 0 либо близок к 0, то уравнение имеет два либо больше вещественных решений
    if (radicand_ >= 0)
    {
        fp_t radical = sqrt(radicand_);

        roots[numberOfRoots] = rootPart + radical;
        roots[numberOfRoots + 1] = rootPart - radical;

        numberOfRoots += 2;
    }
    else if (abs(radicand_) <= EPS)
    {
        roots[numberOfRoots] = rootPart;
        roots[numberOfRoots + 1] = rootPart;

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
    static const fp_t ONE_HALF = static_cast<fp_t>(0.5L);
    static const fp_t ONE_QUARTER = static_cast<fp_t>(0.25L);
    static constexpr fp_t EPS = numeric_limits<fp_t>::epsilon();

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

    fp_t radicalP = sqrt(radicandP);
    fp_t radicalQ = sqrt(radicandQ);

    // Вычисляем итоговые расчетные коэффициенты p_n и q_n (n = 1, 2)
    fp_t p1 = fma(ONE_HALF, a, -radicalP);
    fp_t p2 = fma(ONE_HALF, a, radicalP);
    fp_t q1 = fms(ONE_HALF, u, -sigma, radicalQ);
    fp_t q2 = fms(ONE_HALF, u, sigma, radicalQ);

    // Вычисляем радиканды будущего решения
    fp_t radicand = fma(ONE_QUARTER * p1, p1, -q1);
    fp_t radicand_ = fma(ONE_QUARTER * p2, p2, -q2);

    // Если полученный радиканд (radicand) > 0 либо близок к 0, то уравнение имеет два либо больше вещественных решений
    if (radicand >= 0)
    {
        fp_t radical = sqrt(radicand);

        roots =
        {
            fma(-ONE_HALF, p1, radical),
            fma(-ONE_HALF, p1, -radical)
        };

        numberOfRoots += 2;
    }
    else if (abs(radicand) <= EPS)
    {
        fp_t root = -p1 * ONE_HALF;

        roots =
        {
            root,
            root
        };

        numberOfRoots += 2;
    }

    // Если полученный радиканд (radicand_) > 0 либо близок к 0, то уравнение имеет два либо больше вещественных решений
    if (radicand_ >= 0)
    {
        fp_t radical = sqrt(radicand_);

        roots[numberOfRoots] = fma(-ONE_HALF, p2, radical);
        roots[numberOfRoots + 1] = fma(-ONE_HALF, p2, -radical);;

        numberOfRoots += 2;
    }
    else if (abs(radicand_) <= EPS)
    {
        fp_t root = -p2 * ONE_HALF;

        roots[numberOfRoots] = root;
        roots[numberOfRoots + 1] = root;

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
    static const fp_t ONE_HALF = static_cast<fp_t>(0.5L);
    static const fp_t ONE_QUARTER = static_cast<fp_t>(0.25L);
    static constexpr fp_t EPS = numeric_limits<fp_t>::epsilon();

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

        fp_t subRadical = sqrt(subradicand);

        radicand = fms(static_cast<fp_t>(2.0L), x2, static_cast<fp_t>(2.0L) * sigma, subRadical);
        radicand_ = fms(static_cast<fp_t>(2.0L), x2, static_cast<fp_t>(-2.0L) * sigma, subRadical);
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

        fp_t subRadical = sqrt(subradicand);

        radicand = fma(static_cast<fp_t>(-2.0L) * sigma, subRadical, x2) + x3;
        radicand_ = fma(static_cast<fp_t>(2.0L) * sigma, subRadical, x2) + x3;
    }

    fp_t rootPart = sqrt(r) - C;

    // Если полученный радиканд (radicand) > 0 либо близок к 0, то уравнение имеет два либо больше вещественных решений
    if (radicand >= 0)
    {
        fp_t radical = sqrt(radicand);

        roots =
        {
            rootPart + radical,
            rootPart - radical
        };

        numberOfRoots += 2;
    }
    else if (abs(radicand) <= EPS)
    {
        roots =
        {
            rootPart,
            rootPart
        };

        numberOfRoots += 2;
    }

    rootPart = -sqrt(r) - C;

    // Если полученный радиканд (radicand_) > 0 либо близок к 0, то уравнение имеет два либо больше вещественных решений
    if (radicand_ >= 0)
    {
        fp_t radical = sqrt(radicand_);

        roots[numberOfRoots] = rootPart + radical;
        roots[numberOfRoots + 1] = rootPart - radical;

        numberOfRoots += 2;
    }
    else if (abs(radicand_) <= EPS)
    {
        roots[numberOfRoots] = rootPart;
        roots[numberOfRoots + 1] = rootPart;

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
    static const fp_t ONE_HALF = static_cast<fp_t>(0.5L);
    static const fp_t ONE_QUARTER = static_cast<fp_t>(0.25L);
    static constexpr fp_t EPS = numeric_limits<fp_t>::epsilon();

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

        fp_t subRadical = sqrt(subradicand);

        radicand = fms(static_cast<fp_t>(-2.0L), x2, static_cast<fp_t>(2.0L) * sigma, subRadical);
        radicand_ = fms(static_cast<fp_t>(-2.0L), x2, static_cast<fp_t>(-2.0L) * sigma, subRadical);
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

        fp_t subRadical = sqrt(subradicand);

        radicand = fma(static_cast<fp_t>(-2.0L) * sigma, subRadical, -x2) - x3;
        radicand_ = fma(static_cast<fp_t>(2.0L) * sigma, subRadical, -x2) - x3;
    }

    fp_t sqrtTheta = sqrt(-theta);

    fp_t rootPart = fma(ONE_HALF, sqrtTheta, -C);

    // Если полученный радиканд (radicand) > 0 либо близок к 0, то уравнение имеет два либо больше вещественных решений
    if (radicand >= 0)
    {
        fp_t radical = sqrt(radicand);

        roots =
        {
            fma(ONE_HALF, radical, rootPart),
            fma(-ONE_HALF, radical, rootPart)
        };

        numberOfRoots += 2;
    }
    else if (abs(radicand) <= EPS)
    {
        roots =
        {
            rootPart,
            rootPart
        };

        numberOfRoots += 2;
    }

    rootPart = fma(-ONE_HALF, sqrtTheta, -C);

    // Если полученный радиканд (radicand_) > 0 либо близок к 0, то уравнение имеет два либо больше вещественных решений
    if (radicand_ >= 0)
    {
        fp_t radical = sqrt(radicand_);

        roots[numberOfRoots] = fma(ONE_HALF, radical, rootPart);
        roots[numberOfRoots + 1] = fma(-ONE_HALF, radical, rootPart);

        numberOfRoots += 2;
    }
    else if (abs(radicand_) <= EPS)
    {
        roots[numberOfRoots] = rootPart;
        roots[numberOfRoots + 1] = rootPart;

        numberOfRoots += 2;
    }

    return numberOfRoots;
}

// ======================== ФУНКЦИИ ТЕСТИРОВАНИЯ МЕТОДОВ ======================== //

template<typename fp_t>
void testCubicPolynomial(int testCount, long double maxDistance)
{
    unsigned P = 3; // Степень исходного полинома
    fp_t low = -1, high = 1; // Интервал на котором заданы корни полинома
    fp_t absMaxError, relMaxError; // Абсолютная и относительная погрешность по итогам пройденного теста
    fp_t absMaxErrorTotal = -1, relMaxErrorTotal = -1; // Итоговая максимальная абсолютная и относительная погрешность по итогам всех тестов
    long double absErrorAvg = 0, relErrorAvg = 0; // Средняя абсолютная и относительная погрешность по итогам всех тестов
    unsigned numberOfFoundRoots; // Количество найденных корней
    unsigned cantFind = 0; // Счетчик количества ситуаций, когда методу не удалось найти корни (numberOfFoundRoots == 0)
    vector<fp_t> coefficients(P + 1); // Вектор коэффициентов полинома
    unsigned count = 0; // Счетчик количества ситуаций, когда относительная погрешность больше определенного числа (relMaxError > n)

    for (size_t i = 0; i < testCount; ++i)
    {
        vector<fp_t> foundRoots(P);
        vector<fp_t> trueRoots(P);

        generate_polynomial<fp_t>(P, 0, P, 0, static_cast<fp_t>(maxDistance), low, high, trueRoots, coefficients);

        numberOfFoundRoots = solveCubic<fp_t>(coefficients[3], coefficients[2], coefficients[1], coefficients[0], foundRoots);

        if (numberOfFoundRoots > 0)
        {
            compare_roots<fp_t>(numberOfFoundRoots, P, foundRoots, trueRoots, absMaxError, relMaxError);

            absMaxErrorTotal = absMaxError > absMaxErrorTotal ? absMaxError : absMaxErrorTotal;
            absErrorAvg += absMaxError;

            relMaxErrorTotal = relMaxError > relMaxErrorTotal ? relMaxError : relMaxErrorTotal;
            relErrorAvg += relMaxError;

            count += relMaxError > 1 ? 1 : 0;
        }
        else
            cantFind += 1;
    }

    absErrorAvg /= (testCount - cantFind);
    relErrorAvg /= (testCount - cantFind);

    if (PRINT)
    {
        cout << "CUBIC TEST RESULTS" << endl;
        cout << "========================================" << endl;
        cout << "Max distance: " << maxDistance << endl;
        cout << "Total count of tests: " << testCount << endl;
        cout << "Couldn't find roots: " << cantFind << " times " << endl;
        cout << "Average absolute error: " << absErrorAvg << endl;
        cout << "Total maximum absolute error: " << absMaxErrorTotal << endl;
        cout << "Average relative error: " << relErrorAvg << endl;
        cout << "Total maximum relative error: " << relMaxErrorTotal << endl;
        cout << "RelMaxError > 1: " << count << " times" << endl;
        cout << "========================================" << endl;
    }
}

template<typename fp_t>
void testQuarticPolynomial(int testCount, long double maxDistance)
{
    unsigned P = 4; // Степень исходного полинома
    fp_t low = -1, high = 1; // Интервал на котором заданы корни полинома
    fp_t absMaxError, relMaxError; // Абсолютная и относительная погрешность по итогам пройденного теста
    fp_t absMaxErrorTotal = -1, relMaxErrorTotal = -1; // Итоговая максимальная абсолютная и относительная погрешность по итогам всех тестов
    long double absErrorAvg = 0, relErrorAvg = 0; // Средняя абсолютная и относительная погрешность по итогам всех тестов
    unsigned numberOfFoundRoots; // Количество найденных корней
    unsigned cantFind = 0; // Счетчик количества ситуаций, когда методу не удалось найти корни (numberOfFoundRoots == 0)
    vector<fp_t> coefficients(P + 1); // Вектор коэффициентов полинома
    unsigned count = 0; // Счетчик количества ситуаций, когда относительная погрешность больше определенного числа (relMaxError > n)
    int countExcessRoots = 0;
    int countLostRoots = 0;

    for (size_t i = 0; i < testCount; ++i)
    {
        vector<fp_t> foundRoots(P);
        vector<fp_t> trueRoots(P);
        int excessRoots = 0;
        int lostRoots = 0;

        generate_polynomial<fp_t>(P, 0, 0, 0, static_cast<fp_t>(maxDistance), low, high, trueRoots, coefficients);

        numberOfFoundRoots = euler<fp_t>(coefficients[4], coefficients[3], coefficients[2], coefficients[1], coefficients[0], foundRoots);

        if (numberOfFoundRoots > 0)
        {
            compare_roots<fp_t>(numberOfFoundRoots, P, foundRoots, trueRoots, absMaxError, relMaxError, excessRoots, lostRoots);

            absMaxErrorTotal = absMaxError > absMaxErrorTotal ? absMaxError : absMaxErrorTotal;
            absErrorAvg += absMaxError;

            relMaxErrorTotal = relMaxError > relMaxErrorTotal ? relMaxError : relMaxErrorTotal;
            relErrorAvg += relMaxError;

            countExcessRoots += excessRoots;
            countLostRoots += lostRoots;

            count += relMaxError > 1 ? 1 : 0;
        }
        else
        {
            countLostRoots += 4;
            cantFind += 1;
        }
    }

    absErrorAvg /= (testCount - cantFind);
    relErrorAvg /= (testCount - cantFind);

    if (PRINT)
    {
        cout << "QUARTIC TEST RESULTS" << endl;
        cout << "========================================" << endl;
        cout << "Max distance: " << maxDistance << endl;
        cout << "Total count of tests: " << testCount << endl;
        cout << "Couldn't find roots: " << cantFind << " times " << endl;
        cout << "----------------------------------------" << endl;
        cout << "Average absolute error: " << absErrorAvg << endl;
        cout << "Total maximum absolute error: " << absMaxErrorTotal << endl;
        cout << "Average relative error: " << relErrorAvg << endl;
        cout << "Total maximum relative error: " << relMaxErrorTotal << endl;
        cout << "----------------------------------------" << endl;
        cout << "Total count of lost roots: " << countLostRoots << endl;
        cout << "Total count of excess roots: " << countExcessRoots << endl;
        cout << "----------------------------------------" << endl;
        cout << "relMaxError > 1: " << count << " times" << endl;
        cout << "========================================" << endl;
    }
}

int main()
{
    testQuarticPolynomial<fp_t>(10000000, 1e-5);
}
