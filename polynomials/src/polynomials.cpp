#include "../inc/excerpt.h"
#include "../inc/modules.h"

/*
    Имплементация метода решения линейного уравнения
*/
template<typename fp_t>
unsigned int solveLinear(fp_t a, fp_t b, vector<fp_t>& roots)
{
    roots = { -b / a };

    return isnan(roots[0]) || isinf(roots[0]) ? 0 : 1;
}

template<typename fp_t>
unsigned int solveLinear(fp_t a, fp_t b, vector<complex<fp_t>>& roots)
{
    fp_t root = -b / a;

    if (isfinite(root))
        return 0;
    else
    {
        roots = { complex<fp_t>(root, static_cast<fp_t>(0.0L)) };
        return 1;
    }
}

template<typename fp_t>
unsigned int solveLinear(complex<fp_t> a, complex<fp_t> b, vector<complex<fp_t>>& roots)
{
    complex<fp_t> root(-b / a);

    if (isfinite(root.real()) || isfinite(root.imag()))
        return 0;
    else
    {
        roots = { root };
        return 1;
    }
}

template<typename fp_t>
unsigned int solveQuadratic(fp_t n, fp_t a, fp_t b, vector<fp_t>& roots)
{
    if (isZero(n) || isinf(a /= n))
        return solveLinear(a, b, roots);
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
        return solveLinear(a, b, roots);
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
    if (isZero(n) || isinfc(a /= n))
        return solveLinear(a, b, roots);
    if (isinfc(b /= n))
        return 0;

    a /= static_cast<fp_t>(-2.0L);

    complex<fp_t> d = dscrmt(n, a, b);
    complex<fp_t> sqrtD = sqrt(complex<fp_t>(d));

    complex<fp_t> S(a);
    complex<fp_t> partS(sgn(a) + (isZero(a.real()) && isZero(a.imag())), static_cast<fp_t>(0.0L));

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
    //cout << a << " " << b << " " << c << endl;
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
    //a_ = b - 6 * C * C;
    //b_ = c - 2 * b * C + 8 * C * C * C;
    //c_ = d - c * C + b * C * C - 3 * C * C * C * C;
    //cout << "COEFFS\n" << setprecision(1000) << a_ << "\n" << fma(ONE_QUARTER * a_, a_, -c_) << "\n" << -b_ * b_ * static_cast<fp_t>(0.125L) << endl;
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

    if (abs(radicandP) <= EPS)
    {
        radicandP = static_cast<fp_t>(0.0L);
    }

    if (abs(radicandQ) <= EPS)
    {
        radicandQ = static_cast<fp_t>(0.0L);
    }

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

    //cout << "x2 = " << x2 << endl;
    //cout << "x3 = " << x3 << endl;
    //cout << "y2 = " << y << endl;
    //cout << "theta1 = " << theta << endl;
    //cout << "sigma = " << sigma << endl;

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

/*
    Имплементация метода решения уравнения четвертой степени — Van der Waerden’s Method
    Информация о методе — https://dialnet.unirioja.es/descarga/articulo/6523979.pdf
    Работу выполнил — Погосов Даниэль (https://github.com/DarklleS)
*/
template<typename fp_t>
unsigned int tschirnhaus(fp_t n, fp_t a, fp_t b, fp_t c, fp_t d, vector<fp_t>& roots)
{
    // Нормировка коэффициентов
    if (isZero(n))
        return solveCubic(a, b, c, d, roots);
    if (isinf(a /= n) || isinf(b /= n) || isinf(c /= n) || isinf(d /= n))
        return 0;

    // Объявление констант
    static const fp_t ONE_HALF = static_cast<fp_t>(0.5L);
    static const fp_t ONE_QUARTER = static_cast<fp_t>(0.25L);
    static constexpr fp_t EPS = numeric_limits<fp_t>::epsilon();
    static const complex<fp_t> TWO_C(static_cast<fp_t>(2.0L), static_cast<fp_t>(0.0L));
    static const complex<fp_t> ONE_HALF_C(static_cast<fp_t>(0.5L), static_cast<fp_t>(0.0L));

    // Вычисляем расчетные коэффициенты
    fp_t C = a * ONE_QUARTER;
    fp_t a_ = fma(static_cast<fp_t>(-6.0L) * C, C, b);
    fp_t b_ = fma(fms(static_cast<fp_t>(8.0L) * C, C, static_cast<fp_t>(2.0L), b), C, c);
    fp_t c_ = fma(fma(fma(static_cast<fp_t>(-3.0L) * C, C, b), C, -c), C, d);

    // Счетчик количества вещественных корней
    unsigned numberOfRoots = 0;

    // Вычисляем коэффициенты для бикубического резольвентного уравнения
    fp_t A_ = a_ * ONE_HALF;
    fp_t B_ = -c_;
    fp_t C_ = fms(b_, b_, static_cast<fp_t>(4.0L) * a_, c_) * static_cast<fp_t>(0.125L);

    // Решаем резольвентное бикубическое уравнение вида: f^6 + a_/2 * f^4 - c_ * f^2 + (b_^2 - 4 * a_ * c_) / 8 = 0
    vector<complex<fp_t>> cubicRoots(3);
    unsigned numberOfCubicRoots = solveCubic(static_cast<fp_t>(1.0L), A_, B_, C_, cubicRoots);

    // Определяем решение (наибольшее) бикубического уравнения
    complex<fp_t> ff = cubicRoots[numberOfCubicRoots - 1];
    complex<fp_t> f = sqrt(ff);

    if (isnan(ff.real()) || isnan(ff.imag()))
        return 0;

    // Определение вспомогательных коэффициентов
    complex<fp_t> k = fmac(TWO_C, ff, complex<fp_t>(a_, static_cast<fp_t>(0.0L)));
    complex<fp_t> sqrtK = sqrt(-k);

    // Вычисляем коэффициенты для квадратных вспомогательных уравнений
    complex<fp_t> A = fmac(TWO_C, f, -sqrtK);
    complex<fp_t> B = fmsc(f, k, ONE_HALF_C, complex<fp_t>(-b_, static_cast<fp_t>(0.0L))) * sqrtK / k;

    if (isnan(B.real()))
    {
        return 0;
    }

    // Решаем резольвентное бикубическое уравнение вида: d^2 + (2 * f - sqrt(-k)) * f - sqrt(-k) * ((2 * f * k + b) / (2 * k)) = 0
    vector<complex<fp_t>> quadraticRoots(2);
    unsigned numberOfQuadraticRoots = solveQuadratic(complex<fp_t>(static_cast<fp_t>(1.0L), static_cast<fp_t>(0.0L)), A, -B, quadraticRoots);

    vector<complex<fp_t>> rootsComplex =
    {
        quadraticRoots[0] + f - C,
        quadraticRoots[1] + f - C
    };

    // Отбрасываем комплексные решения
    selectRealRoots(rootsComplex, numberOfRoots, roots);

    // Переопределяем коэффициент при d для нового квадратного уравнения
    A = fmac(TWO_C, f, sqrtK);

    // Решаем резольвентное бикубическое уравнение вида: d^2 + (2 * f + sqrt(-k)) * f + sqrt(-k) * ((2 * f * k + b) / (2 * k)) = 0
    numberOfQuadraticRoots = solveQuadratic(complex<fp_t>(static_cast<fp_t>(1.0L), static_cast<fp_t>(0.0L)), A, B, quadraticRoots);


    rootsComplex =
    {
        quadraticRoots[0] + f - C,
        quadraticRoots[1] + f - C
    };

    // Отбрасываем комплексные решения
    selectRealRoots(rootsComplex, numberOfRoots, roots);

    return numberOfRoots;
}

// Комплексная версия алгоритма Феррари
template<typename fp_t>
unsigned int ferrari_complex(fp_t n, fp_t a, fp_t b, fp_t c, fp_t d, vector<complex<fp_t>>& roots)
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

    // Радикал R
    complex<fp_t> R = sigma * sqrt(complex<fp_t>(subR, static_cast<fp_t>(0.0L)));

    complex<fp_t> radicandPart(fms(-ONE_HALF, m, ONE_HALF, a_), static_cast<fp_t>(0.0L));

    // Вычисляем радиканды будущего решения
    complex<fp_t> radicand = radicandPart - R;
    complex<fp_t> radicand_ = radicandPart + R;

    fp_t sqrtM = sqrt(m * ONE_HALF);

    complex<fp_t> rootPart(sqrtM - C, static_cast<fp_t>(0.0L));
    complex<fp_t> rootPart_(-sqrtM - C, static_cast<fp_t>(0.0L));

    complex<fp_t> radical(sqrt(radicand));
    complex<fp_t> radical_(sqrt(radicand_));

    roots =
    {
        rootPart + radical,
        rootPart - radical,
        rootPart_ + radical_,
        rootPart_ - radical_
    };

    // Определяем кол-во вещественных корней
    if (radicand.real() >= 0 && isZero(radicand.imag()))
        numberOfRoots += 2;
    if (radicand_.real() >= 0 && isZero(radicand_.imag()))
        numberOfRoots += 2;

    return numberOfRoots;
}

/*
    Имплементация метода решения уравнения четвертой степени — FQS method
    Информация о методе — https://www.academia.edu/es/27122162/The_fast_quartic_solver
    Работу выполнил — Погосов Даниэль (https://github.com/DarklleS) при содействии Алексея Федорова (https://github.com/Alexis777777)
*/
template<typename fp_t>
unsigned int fqs(fp_t n, fp_t a, fp_t b, fp_t c, fp_t d, vector<fp_t>& roots) {

    // Нормировка коэффициентов
    if (isZero(n) || isinf(a /= n))
        return solveCubic(a, b, c, d, roots);
    if (isinf(b /= n))
        return 0;
    if (isinf(c /= n))
        return 0;
    if (isinf(d /= n))
        return 0;

    vector<complex<fp_t>> roots1(2);
    vector<complex<fp_t>> roots2(2);

    unsigned numberOfRoots = 0;

    // Объявление констант
    static const fp_t ONE_HALF = static_cast<fp_t>(0.5L);
    static const fp_t ONE_6TH = static_cast<fp_t>(1.0L / 6.0L);

    //Case 1. Два двойных корня / Один четвертичный корень(4 - кратный корень).

    // Определяем коэффициенты вспомогательного уравнения
    fp_t alpha_ = a * ONE_HALF;
    fp_t beta_ = fma(-alpha_, alpha_, b) * ONE_HALF; // (b - alpha_^2) / 2

    // Определяющие коэффициенты
    fp_t E1 = fma(static_cast<fp_t>(-2.0L) * alpha_, beta_, c); // c - 2 * alpha_ * beta
    fp_t E2 = fma(-beta_, beta_, d); // d - beta_^2

    if (isZero(E1) and isZero(E2)) {
        // Поиск корней вспомогательного уравнения
        solveQuadratic<fp_t>(static_cast<fp_t>(1.0L), alpha_, beta_, roots1);

        //Определим действительные корни
        if (isZero(roots1[0].imag())) {
            roots[numberOfRoots] = roots1[0].real();
            roots[numberOfRoots + 1] = roots1[0].real();
            numberOfRoots += 2;
        }
        if (isZero(roots1[1].imag())) {
            roots[numberOfRoots] = roots1[1].real();
            roots[numberOfRoots + 1] = roots1[1].real();
            numberOfRoots += 2;
        }

        return numberOfRoots;
    }

    //Case 2. Один тройной корень и один простой корень.

    // Определяем коэффициенты вспомогательного уравнения
    alpha_ = a * ONE_HALF; beta_ = b * ONE_6TH;

    // Поиск корней вспомогательного уравнения
    vector<fp_t> roots_(2);
    solveQuadratic<fp_t>(static_cast<fp_t>(1.0L), alpha_, beta_, roots_);

    // Определяем корни исходного уравнения
    fp_t x1 = roots_[0], x2 = -fma(static_cast<fp_t>(3.0L), x1, a); // x2 = -3 * x1 - a
    // Определяющие коэффициенты
    E1 = fma(x1 * x1, fma(static_cast<fp_t>(3.0L), x2, x1), c); // c + x1 * x1 * (x1 + 3 * x2);
    E2 = fma(-x1 * x1, x1 * x2, d); // d - x^3 * x2 

    if (isZero(E1) and isZero(E2)) {
        roots = { x1, x2, x2, x2 }; return numberOfRoots += 4;
    }

    // Определяем корни исходного уравнения
    x1 = roots_[1]; x2 = -fma(static_cast<fp_t>(3.0L), x1, a); // -3 * x1 - a
    // Определяющие коэффициенты
    E1 = fma(x1 * x1, fma(static_cast<fp_t>(3.0L), x2, x1), c); // c + x1 * x1 * (x1 + 3 * x2);
    E2 = fma(-x1 * x1, x1 * x2, d); // d -  x1^3 * x2

    if (isZero(E1) and isZero(E2)) {
        roots = { x1, x2, x2, x2 }; return numberOfRoots += 4;
    }

    //Case 3. Один двойной корень и два простых корня / 4 простых корня. 
    vector<fp_t> coeffs(4);
    // Определяем коэффициенты вспомогательных уравнений с помощью итеративного метода
    FQScoeffs(n, a, b, c, d, coeffs);
    fp_t alpha = coeffs[0], beta = coeffs[1], gamma = coeffs[2], delta = coeffs[3];

    // Поиск корней вспомогательного уравнения
    solveQuadratic<fp_t>(static_cast<fp_t>(1.0L), alpha, beta, roots1);
    solveQuadratic<fp_t>(static_cast<fp_t>(1.0L), gamma, delta, roots2);

    //Определим действительные корни
    if (isZero(roots1[0].imag())) {
        roots[numberOfRoots] = roots1[0].real();    numberOfRoots++;
    }
    if (isZero(roots1[1].imag())) {
        roots[numberOfRoots] = roots1[1].real();    numberOfRoots++;
    }
    if (isZero(roots2[0].imag())) {
        roots[numberOfRoots] = roots2[0].real();    numberOfRoots++;
    }
    if (isZero(roots2[1].imag())) {
        roots[numberOfRoots] = roots2[1].real();    numberOfRoots++;
    }

    return numberOfRoots;
}

/*
    Имплементация метода решения уравнения четвертой степени — Mansfield Merriman
    Информация о методе — https://www.jstor.org/stable/2369666
    Работу выполнил — Погосов Даниэль (https://github.com/DarklleS) при содействии Артёма Успенского (https://github.com/MadBunny0)
*/
template<typename fp_t>
unsigned int merriman(fp_t n, fp_t a_, fp_t b_, fp_t c_, fp_t d_, vector<fp_t>& roots) {
    // Нормировка коэффициентов
    if (isZero(n) || isinf(a_ /= n))
        return solveCubic(a_, b_, c_, d_, roots);
    if (isinf(b_ /= n))
        return 0;
    if (isinf(c_ /= n))
        return 0;
    if (isinf(d_ /= n))
        return 0;

    // Объявление констант
    static const fp_t ONE_HALF = static_cast<fp_t>(0.5L);
    static const fp_t ONE_QUARTER = static_cast<fp_t>(0.25L);
    static const fp_t NINE_QUARTER = static_cast<fp_t>(2.25L);
    static const fp_t ONE_SIXTH = static_cast<fp_t>(1.0L / 6.0L);
    static const fp_t ONE_THIRD = static_cast<fp_t>(1.0L / 3.0L);
    static const fp_t FOUR_THIRD = static_cast<fp_t>(4.0L / 3.0L);
    static const fp_t TWO = static_cast<fp_t>(2.0L);
    static const fp_t THREE = static_cast<fp_t>(3.0L);
    static const fp_t SIX = static_cast<fp_t>(6.0L);
    static const fp_t NINE = static_cast<fp_t>(9.0L);
    static const fp_t EPS = static_cast<fp_t>(1e-6);

    // Пересчёт переменных
    fp_t a = a_ * ONE_QUARTER;
    fp_t b = b_ * ONE_SIXTH;
    fp_t c = c_ * ONE_QUARTER;
    fp_t d = d_;

    // Количество вещественных корней
    unsigned numberOfRoots = 0;
    //vector<fp_t> roots(4);

    // Вычисляем расчетные коэффициенты
    fp_t m = fma(b * b, b, fms(d, fma(a, a, -b), -c, fma(-2 * a, b, c))); // a^2 * d + b^3 + c^2 - 2abc - bd
    fp_t nn = pow(fma(b, b, fms(ONE_THIRD, d, FOUR_THIRD * a, c)), THREE); // (b^2 + 1/3d - 4/3ac)^3

    // переменная для развитвления логики m^2 - n: приводит к разным u,v,w
    fp_t radical = fma(m, m, -nn);
    // переменная для развитвления логики 2a^3 - 3ab + c: определяет знаки в корнях
    fp_t cond_coef = fms(TWO, pow(a, THREE), THREE, a * b) + c; // 2 * a^3 - 3 * a * b + c

    if (radical > 0) {
        fp_t s1 = m + sqrt(radical);
        fp_t s = ONE_HALF * pow(s1, ONE_THIRD); // 1/2 (m + sqrt(m^2 - n))^1/3

        fp_t j = (m - sqrt(radical) < 0) ? -1 : 1; // определяем знак для t1
        fp_t t1 = j * (m - sqrt(radical));

        fp_t t = j * ONE_HALF * pow(t1, ONE_THIRD); // 1/2 (m - sqrt(m^2 - n))^1/3

        // Вычисляем коэффициенты решения
        fp_t u = fma(a, a, -b) + s + t; // a^2 - b + s + t
        fp_t v = fms(TWO, pow(a, TWO), TWO, b) - s - t; // 2a^2 - 2b - s - t
        fp_t w = fms(THREE, pow(s, TWO), SIX * s, t) +
            fma(v, v, THREE * pow(t, TWO)); // v^2 + 3s^2 - 6st + 3t^2

        // вычисляем слагаемые исходного решения
        fp_t roots_term_1 = sqrt(u);
        fp_t roots_term_2 = sqrt(v + sqrt(w));
        fp_t roots_term_3 = sqrt(v - sqrt(w));
        // все корни комплексные
        if (isnan(roots_term_1) || (isnan(roots_term_2) && isnan(roots_term_3))) {
            //cout << "All roots are complex CASE 1" << endl;
            return 0;
        }
        // вещественные корни, 2 корня кратные
        else if (isnan(roots_term_3)) {
            if (cond_coef < 0) {
                roots[0] = -a + roots_term_1 + roots_term_2;
                roots[1] = -a + roots_term_1 - roots_term_2;
                roots[2] = -a - roots_term_1; //УБРАЛИ МНИМУЮ ЧАСТЬ
                roots[3] = -a - roots_term_1; //УБРАЛИ МНИМУЮ ЧАСТЬ
            }
            else {
                roots[0] = -a - roots_term_1 - roots_term_2;
                roots[1] = -a - roots_term_1 + roots_term_2;
                roots[2] = -a + roots_term_1; //УБРАЛИ МНИМУЮ ЧАСТЬ
                roots[3] = -a + roots_term_1; //УБРАЛИ МНИМУЮ ЧАСТЬ
            }
        }
        // вещественные корни, 2 корня кратные
        else if (isnan(roots_term_2)) {
            if (cond_coef < 0) {
                roots[0] = -a + roots_term_1; //УБРАЛИ МНИМУЮ ЧАСТЬ
                roots[1] = -a + roots_term_1; //УБРАЛИ МНИМУЮ ЧАСТЬ
                roots[2] = -a - roots_term_1 + roots_term_3;
                roots[3] = -a - roots_term_1 - roots_term_3;
            }
            else {
                roots[0] = -a - roots_term_1; //УБРАЛИ МНИМУЮ ЧАСТЬ
                roots[1] = -a - roots_term_1; //УБРАЛИ МНИМУЮ ЧАСТЬ
                roots[2] = -a + roots_term_1 - roots_term_3;
                roots[3] = -a + roots_term_1 + roots_term_3;
            }
        }
        // вещественные корни, 4 различных корня
        else {
            if (cond_coef < 0) {
                roots[0] = -a + roots_term_1 + roots_term_2;
                roots[1] = -a + roots_term_1 - roots_term_2;
                roots[2] = -a - roots_term_1 + roots_term_3;
                roots[3] = -a - roots_term_1 - roots_term_3;
            }
            else {
                roots[0] = -a - roots_term_1 - roots_term_2;
                roots[1] = -a - roots_term_1 + roots_term_2;
                roots[2] = -a + roots_term_1 - roots_term_3;
                roots[3] = -a + roots_term_1 + roots_term_3;
            }
            //cout << "ibanutiy else";
        }
        numberOfRoots = 4;
        //cout << "m * m - nn > 0" << endl;
    }
    // если мнимая часть слишком мала (на уровне погрешности fp_t)
    else if (abs(radical) <= EPS) {
        // вычисляем расчетные коэффициенты для исходного решения
        fp_t u = fma(a, a, -b) + pow(nn, ONE_SIXTH); // a^2 - b + n^1/6
        fp_t v = fms(TWO, pow(a, TWO), TWO, b) - pow(nn, ONE_SIXTH); // 2a^2 - 2b - n^1/6

        if (isZero(v)) {
            // вычисляем слагаемое исходного решения
            fp_t roots_term = sqrt(fms(THREE, pow(a, TWO), THREE, b)); // sqrt(fms)

            if (cond_coef < 0) {
                roots[0] = -a + roots_term;
                roots[1] = roots[0];
                roots[2] = -a - roots_term;
                roots[3] = roots[2];
            }
            else {
                roots[0] = -a - roots_term;
                roots[1] = roots[0];
                roots[2] = -a + roots_term;
                roots[3] = roots[2];
            }
            //cout << "m * m - nn == 0, v == 0" << endl;
            numberOfRoots = 4;
        }
        else {
            // вычисляем слагаемые исходного решения
            fp_t roots_term_1 = sqrt(u);
            fp_t roots_term_2 = sqrt(TWO * v);


            if (cond_coef < 0) {
                roots[0] = -a + roots_term_1 + roots_term_2;
                roots[1] = -a + roots_term_1 - roots_term_2;
                roots[2] = -a - roots_term_1;
                roots[3] = roots[2];
            }
            else {
                roots[0] = -a - roots_term_1 - roots_term_2;
                roots[1] = -a - roots_term_1 + roots_term_2;
                roots[2] = -a + roots_term_1;
                roots[3] = roots[2];

            }
            //cout << "m * m - nn == 0, v != 0" << endl;
            numberOfRoots = 4;
        }
    }
    else if (isZero(m)) {
        // вычисляем расчетные коэффициенты для исходного решения
        fp_t u = fma(a, a, -b); // a^2 - b
        fp_t v = fms(TWO, pow(a, TWO), TWO, b); // 2a^2 - 2b
        fp_t w = fms(v, v, THREE, pow(b, TWO)) - d; // v^2 - 3b^2 - d

        // вычисляем слагаемые исходного решения
        fp_t roots_term_1 = sqrt(u);
        fp_t roots_term_2 = sqrt(v + sqrt(w));
        fp_t roots_term_3 = sqrt(v - sqrt(w));


        if (isnan(roots_term_1) || (isnan(roots_term_2) && isnan(roots_term_3))) {
            //cout << "ALL of roots are complex CASE 4" << endl;
            return 0;
        }
        else if (isnan(roots_term_2)) {
            if (cond_coef < 0) {
                roots[0] = -a - roots_term_1;
                roots[1] = -a - roots_term_1;
                roots[2] = -a + roots_term_1 + roots_term_3;
                roots[3] = -a + roots_term_1 - roots_term_3;
            }
            else {
                roots[0] = -a + roots_term_1;
                roots[1] = -a + roots_term_1;
                roots[2] = -a - roots_term_1 - roots_term_3;
                roots[3] = -a - roots_term_1 + roots_term_3;
            }
            //cout << "ALL of roots are complex term_2                  ";

        }
        else if (isnan(roots_term_3)) {
            if (cond_coef < 0) {
                roots[0] = -a - roots_term_1 + roots_term_2;
                roots[1] = -a - roots_term_1 - roots_term_2;
                roots[2] = -a + roots_term_1; //у
                roots[3] = -a + roots_term_1;
            }
            else {
                roots[0] = -a + roots_term_1 - roots_term_2;
                roots[1] = -a + roots_term_1 + roots_term_2;
                roots[2] = -a - roots_term_1;
                roots[3] = -a - roots_term_1;
            }
            //cout << "ALL of roots are complex term_3                   ";

        }
        //cout << "m == 0" << endl;
        numberOfRoots = 4;
    }
    // не попали ни в 1 кейс -> не удалось найти корни аналитически, возвращаем 0
    else {
        //cout << "Irreducible case" << endl;

        return 0;
    }
    //cout << roots[0] << " " << roots[1] << " " << roots[2] << " " << roots[3] << endl;
    return numberOfRoots;
}

// Функция поиска корней в случае когда (p^2 - a0) = 0, т. е. p = p1 = sqrt(a0)
// Коэффициенты q, q1 могут быть как действительными, так и комплексносопряженными
// Функция получает на вход коэффициенты a3 и a2, коэффициент p, ссылку на вектор корней,
// возвращает число найденных действительных корней
template<typename fp_t>
unsigned int findRoots(fp_t a3, fp_t a2, fp_t p, vector<fp_t>& roots)
{
    vector<complex<fp_t>> qroots(2);
    solveQuadratic<fp_t>(1, -a3, a2 - 2 * p, qroots);

    unsigned int numberOfRoots = 0;

    if (isZero(qroots[0].imag()) and isZero(qroots[1].imag())) {
        fp_t q = qroots[0].real(), q1 = qroots[1].real();

        vector<fp_t> quadraticRoots(2);
        unsigned int nRoots = 0;

        //fp_t q = fma(0.5, a3, -sr), q1 = fma(0.5, a3, sr);

        nRoots = solveQuadratic<fp_t>(1, q, p, quadraticRoots);
        roots[numberOfRoots] = quadraticRoots[0] * (nRoots > 0);
        numberOfRoots += (nRoots > 0);

        roots[numberOfRoots] = quadraticRoots[1] * (nRoots > 0);
        numberOfRoots += (nRoots > 0);

        nRoots = solveQuadratic<fp_t>(1, q1, p, quadraticRoots);
        roots[numberOfRoots] = quadraticRoots[0] * (nRoots > 0);
        numberOfRoots += (nRoots > 0);

        roots[numberOfRoots] = quadraticRoots[1] * (nRoots > 0);
        numberOfRoots += (nRoots > 0);

        return numberOfRoots;
    }
    else
    {
        complex<fp_t> q = qroots[0], q1 = qroots[1];
        vector<complex<fp_t>> rootsComplex(2);

        complex<fp_t> pc = p;
        complex<fp_t> one; one.real(1);

        solveQuadratic(one, q, pc, rootsComplex);
        roots[numberOfRoots] =
            rootsComplex[0].real() * (isZero(rootsComplex[0].imag()));
        numberOfRoots += (isZero(rootsComplex[0].imag()));

        roots[numberOfRoots] =
            rootsComplex[1].real() * (isZero(rootsComplex[1].imag()));
        numberOfRoots += (isZero(rootsComplex[1].imag()));

        solveQuadratic(one, q1, pc, rootsComplex);
        roots[numberOfRoots] =
            rootsComplex[0].real() * (isZero(rootsComplex[0].imag()));
        numberOfRoots += (isZero(rootsComplex[0].imag()));

        roots[numberOfRoots] =
            rootsComplex[1].real() * (isZero(rootsComplex[1].imag()));
        numberOfRoots += (isZero(rootsComplex[1].imag()));

        return numberOfRoots;
    }
}

/*
    Имплементация метода решения уравнения четвёртой степени — Solution of quartic equations by William Squire
    Информация о методе — https://www.tandfonline.com/doi/abs/10.1080/0020739790100223
    Работу выполнил — Погосов Даниэль (https://github.com/DarklleS) при содействии Азизхона Ишаннова (https://github.com/ZeekFeed07)
*/
template <typename fp_t>
unsigned int squire(fp_t a, fp_t b, fp_t c, fp_t d, fp_t e, std::vector<fp_t>& roots)
{
    // Нормировка коэффициентов
    if (isZero(a) || isinf(b /= a))
        return solveCubic(b, c, d, e, roots);
    if (isinf(c /= a))
        return 0;
    if (isinf(d /= a))
        return 0;
    if (isinf(e /= a))
        return 0;
    // переопределяем исходные коэффициенты
    fp_t a0 = e, a1 = d, a2 = c, a3 = b;
    unsigned int numberOfRoots = 0;

    //if (abs(a0) > 1 or abs(a1) > 1 or abs(a2) > 1 or abs(a3) > 1) return 0;

    fp_t p; // переменная полиноме 6 степени
    fp_t bd = 0; // заданная граница приближения
    fp_t coef[6];

    // Вычисление коэффициентов полинома 6 степени функции S(p)
    coef[0] = pow(a0, 3);
    if (coef[0] < bd) bd = coef[0];
    coef[1] = -a2 * a0 * a0;
    if (coef[1] < bd) bd = coef[1];
    coef[2] = a0 * fma(a1, a3, -a0); // a0*(a1*a3-a0)
    if (coef[2] < bd) bd = coef[2];
    coef[3] = fma(-a0, a3 * a3, fms<fp_t>(2 * a0, a2, a1, a1)); // 2*a0*a2-a1^2 +- a0*a3^2
    if (coef[3] < bd) bd = coef[3];
    coef[4] = fma(a1, a3, -a0); // a1*a3-a0
    if (coef[4] < bd) bd = coef[4];
    coef[5] = -a2;
    if (coef[5] < bd) bd = coef[5];

    // Функция S(p) = p^6+coef[5]*p^5+coef[5]*p^4+coef[4]*p^5+coef[3]*p^3+coef[2]*p^2+coef[1]*p+coef[0]
    auto fun = +[](fp_t p, fp_t* coef) {
        return static_cast<fp_t> (fma(fma(fma(fma(fma(
            coef[5] + p, p, coef[4]), p, coef[3]), p, coef[2]), p, coef[1]), p, coef[0]));
    };

    // переопределяем границу приближения
    bd = 1 - bd;

    fp_t eps1 = eps_temp<fp_t>;

    if (a0 > 0) {
        // определяем решение уравнения 6-й степени, если a0 > 0 и abs(fun(p, coef)) < eps1
        p = sqrt(a0);
        // подставляем решение p в уравнений 6-й степени и сравниваем с величиной близкой к 0
        if (abs(fun(p, coef)) < eps1)
        {
            // находим корни исходного уравнения
            return findRoots(a3, a2, p, roots);
        }

        // находим корень методом бисекции (a0 > 0)
        p = refine<fp_t>(p, bd, coef, fun);
    }
    else p = refine<fp_t>(0, bd, coef, fun); // находим корень методом бисекции (a0 <= 0)

    fp_t den = fma(p, p, -a0); // p*p-a0
    fp_t eps2 = eps_temp<fp_t>;

    if (abs(den) < eps2)
    {
        // находим корень методом бисекции
        return findRoots(a3, a2, p, roots);
    }

    // Нахождение коэффициентов при квадратных уравнениях
    fp_t p1 = a0 / p;
    fp_t q = p * fma(a3, p, -a1) / den; // p * (a3*p - a1) / den
    fp_t q1 = fms(a1, p, a3, a0) / den; // (a1*p - a3*a0) / den

    vector<fp_t> quadraticRoots(2);
    unsigned int nRoots = 0;

    // Нахождение корней квадратных уравнений
    nRoots = solveQuadratic<fp_t>(1, q, p, quadraticRoots);
    roots[numberOfRoots] = quadraticRoots[0] * (nRoots > 0);
    numberOfRoots += (nRoots > 0);
    roots[numberOfRoots] = quadraticRoots[1] * (nRoots > 0);
    numberOfRoots += (nRoots > 0);

    nRoots = solveQuadratic<fp_t>(1, q1, p1, quadraticRoots);
    roots[numberOfRoots] = quadraticRoots[0] * (nRoots > 0);
    numberOfRoots += (nRoots > 0);
    roots[numberOfRoots] = quadraticRoots[1] * (nRoots > 0);
    numberOfRoots += (nRoots > 0);

    return numberOfRoots;
}

/*
    Имплементация метода решения уравнения четвертой степени — Ungar Method
    Информация о методе — https://core.ac.uk/download/pdf/82351384.pdf
    Работу выполнил — Погосов Даниэль (https://github.com/DarklleS) при содействии Погоса Анесяна (https://github.com/Pogosito)
*/
template<typename fp_t>
unsigned int ungar(fp_t n, fp_t a, fp_t b, fp_t c, fp_t d, vector<fp_t>& roots)
{
    unsigned int numberOfRoots = 0;

    // MARK: - Объявление констант

    static const fp_t EIGHT_THIRDS = static_cast<fp_t>(8.0L / 3.0L);
    static const fp_t ONE_FOURTH = static_cast<fp_t>(0.25L);

    static const complex<fp_t> COMPLEX_ONE_FOURTH = complex<fp_t>(ONE_FOURTH, static_cast<fp_t>(0.0L));
    static const complex<fp_t> COMPLEX_FOUR_THIRDS = complex<fp_t>(static_cast<fp_t>(4.0L / 3.0L), static_cast<fp_t>(0.0L));
    static const complex<fp_t> COMPLEX_HALF = complex<fp_t>(static_cast<fp_t>(0.5L), static_cast<fp_t>(0.0L));

    static const complex<fp_t> q0 = complex<fp_t>(static_cast<fp_t>(1.0L), static_cast<fp_t>(0.0L));
    static const complex<fp_t> q1 = complex<fp_t>(static_cast<fp_t>(-0.5L), static_cast<fp_t>(0.5L) * sqrt(static_cast<fp_t>(3.0L)));
    static const complex<fp_t> q2 = complex<fp_t>(static_cast<fp_t>(-0.5L), static_cast<fp_t>(-0.5L) * sqrt(static_cast<fp_t>(3.0L)));

    // MARK: - Вычисляем расчетные коэффициенты P, Q, R, alpha0, betta0, gamma0

    // P
    fp_t P_helper = fms(static_cast<fp_t>(4.0L), -b, -a, a); // -4b + a^2
    fp_t P = fms(static_cast<fp_t>(8.0L), c, -a, P_helper); // 8c + a*P_helper

    // Q
    fp_t Q_helper = fms(static_cast<fp_t>(4.0L), d, a, c); // 4d - ac
    fp_t Q = fms(static_cast<fp_t>(3.0L), Q_helper, -b, b); // 3Q_helper + b^2

    // R
    fp_t R_helper0 = fms(static_cast<fp_t>(3.0L) * d, a, b, c); // 3ad-bc
    fp_t R_helper1 = fms(R_helper0, a, c * c, -static_cast<fp_t>(3.0L)); // R_helper0*a + 3c^2
    fp_t R_helper2 = fms(b, b, static_cast<fp_t>(36.0L), d); // b^2 - 36d
    fp_t R = fms(static_cast<fp_t>(2.0L) * R_helper2, b, -R_helper1, static_cast<fp_t>(9.0L)); // 2 * R_helper2 * b + 9 * R_helper1

    fp_t alpha0 = fms(a, a, EIGHT_THIRDS, b); // a^2 - 8 * b / 3

    fp_t betta_gamma_helper = fma(ONE_FOURTH, R * R, -pow(Q, static_cast<fp_t>(3.0L))); // R^2 / 4 - Q^3

    complex<fp_t> sqrt_betta_gamma_helper = (betta_gamma_helper > 0) ? sqrt(betta_gamma_helper) : complex<fp_t>(0, sqrt(abs(betta_gamma_helper)));
    complex<fp_t> betta_helper_1 = fmac(COMPLEX_HALF, complex<fp_t>(R, 0), sqrt_betta_gamma_helper); // R/2 + sqrt_betta_gamma_helper
    complex<fp_t> gamma_helper_1 = fmac(COMPLEX_HALF, complex<fp_t>(R, 0), -sqrt_betta_gamma_helper); // R/2 - sqrt_betta_gamma_helper

    // MARK: Находим все значения кубического корня (R + sqrt(R^2 - 4Q^3)) / 2. В статье ошибка (см. файл с опечатками)
    vector<complex<fp_t>> betta0CubeRoots = cubeRoot(betta_helper_1);
    vector<complex<fp_t>> gamma0CubeRoots = cubeRoot(gamma_helper_1);

    // MARK: - Согласно теореме 5.1 считаем количество комплексных и действительных корней.
    // 1) R*R - 4*pow(Q,3) > 0  --- > 2 действит. + 2 комплексных корня
    // 2) R*R - 4*pow(Q,3) = 0  --- > 2 действит. + (2 комплексных <=> T >= 0 для всех cuberoots(betta_helper_1))
    // 3) R*R - 4*pow(Q,3) < 0  --- > 3.1) 4 действит. <=> T >= 0 для всех cuberoots(betta_helper_1)
    int cnt_negative_T = 0;
    int cnt_pozitive_T = 0;

    std::vector<fp_t> T;

    for (int i = 0; i < 3; ++i)
    {
        fp_t _3a1_8b = fms(static_cast<fp_t>(3.0L) * a, a, static_cast<fp_t>(8.0L), b); // 3*a^2 - 8*b
        fp_t t = fma(static_cast<fp_t>(8.0L), real(betta0CubeRoots[i]), _3a1_8b); // 8 * Re(betta0CubeRoots[i]) + _3a1_8b
        // MARK: Используем round т.к никогда не получаем чистый 0. (Проверка с eps не проходит. При вычитании двух одинаковых чисел получали не 0, но при этом результат больше eps)
        if (round(abs(static_cast<fp_t>(8.0L) * real(betta0CubeRoots[i]) + _3a1_8b)) <= numeric_limits<fp_t>::epsilon()) { t = 0; }
        T.push_back(t);
        if (T[i] < 0) cnt_negative_T++;
        if (T[i] >= 0) cnt_pozitive_T++;
    }

    // MARK: Используем round т.к никогда не получаем чистый 0. (Проверка с eps не проходит)

    // Исправление
    numberOfRoots = (round(abs(betta_gamma_helper)) <= numeric_limits<fp_t>::epsilon() && cnt_pozitive_T == 3) ? 4 :
        (betta_gamma_helper > 0) ? 2 : (cnt_pozitive_T == 3) ? 4 : 0;

    // Провекра условия 4.10
    int beta0_I = 0;
    int gamma0_J = 0;
    check4_10Condition(betta0CubeRoots, gamma0CubeRoots, Q, beta0_I, gamma0_J);

    complex<fp_t> alpha0_complex = complex<fp_t>(alpha0, static_cast<fp_t>(0.0L));
    // +-sqrt(4/3*(q0*betta0CubeRoots[beta0_I] + q0*gamma0CubeRoots[gamma0_J])+alpha0_complex)
    vector<complex<fp_t>> bettaSqrts = { sqrt(fmac(COMPLEX_FOUR_THIRDS, fmsc(q0, betta0CubeRoots[beta0_I], -q0, gamma0CubeRoots[gamma0_J]), alpha0_complex)),
                                        -sqrt(fmac(COMPLEX_FOUR_THIRDS, fmsc(q0, betta0CubeRoots[beta0_I], -q0, gamma0CubeRoots[gamma0_J]), alpha0_complex)) };
    // +-sqrt(4/3*(q1*betta0CubeRoots[beta0_I] + q2*gamma0CubeRoots[gamma0_J])+alpha0_complex)
    vector<complex<fp_t>> gammaSqrts = { sqrt(fmac(COMPLEX_FOUR_THIRDS, fmsc(q1, betta0CubeRoots[beta0_I], -q2, gamma0CubeRoots[gamma0_J]), alpha0_complex)),
                                        -sqrt(fmac(COMPLEX_FOUR_THIRDS, fmsc(q1, betta0CubeRoots[beta0_I], -q2, gamma0CubeRoots[gamma0_J]), alpha0_complex)) };
    // +-sqrt(4/3*(q2*betta0CubeRoots[beta0_I] + q1*gamma0CubeRoots[gamma0_J])+alpha0_complex)
    vector<complex<fp_t>> deltaSqrts = { sqrt(fmac(COMPLEX_FOUR_THIRDS, fmsc(q2, betta0CubeRoots[beta0_I], -q1, gamma0CubeRoots[gamma0_J]), alpha0_complex)),
                                        -sqrt(fmac(COMPLEX_FOUR_THIRDS, fmsc(q2, betta0CubeRoots[beta0_I], -q1, gamma0CubeRoots[gamma0_J]), alpha0_complex)) };
    // Провекра условия 4.11
    int bettaI = 0;
    int gammaJ = 0;
    int deltaK = 0;
    check4_11Condition(bettaSqrts, gammaSqrts, deltaSqrts, P, bettaI, gammaJ, deltaK);

    // MARK: - Находим корни

    fp_t alpha = -ONE_FOURTH * a;
    complex<fp_t> alpha_complex = complex<fp_t>(alpha, static_cast<fp_t>(0.0L));
    // (deltaSqrts[deltaK]/4 + (gammaSqrts[gammaJ]/4+(bettaSqrts[bettaI]/4 + alpha_complex)))
    complex<fp_t> w0 = fmac(COMPLEX_ONE_FOURTH, deltaSqrts[deltaK], fmac(COMPLEX_ONE_FOURTH, gammaSqrts[gammaJ], fmac(COMPLEX_ONE_FOURTH, bettaSqrts[bettaI], alpha_complex)));
    // (-deltaSqrts[deltaK]/4 + (-gammaSqrts[gammaJ]/4+(bettaSqrts[bettaI]/4 + alpha_complex)))
    complex<fp_t> w1 = fmac(-COMPLEX_ONE_FOURTH, deltaSqrts[deltaK], fmac(-COMPLEX_ONE_FOURTH, gammaSqrts[gammaJ], fmac(COMPLEX_ONE_FOURTH, bettaSqrts[bettaI], alpha_complex)));
    // (-deltaSqrts[deltaK]/4 + (gammaSqrts[gammaJ]/4+(-bettaSqrts[bettaI]/4 + alpha_complex)))
    complex<fp_t> w2 = fmac(-COMPLEX_ONE_FOURTH, deltaSqrts[deltaK], fmac(COMPLEX_ONE_FOURTH, gammaSqrts[gammaJ], fmac(-COMPLEX_ONE_FOURTH, bettaSqrts[bettaI], alpha_complex)));
    // (deltaSqrts[deltaK]/4 + (-gammaSqrts[gammaJ]/4+(-bettaSqrts[bettaI]/4 + alpha_complex)))
    complex<fp_t> w3 = fmac(COMPLEX_ONE_FOURTH, deltaSqrts[deltaK], fmac(-COMPLEX_ONE_FOURTH, gammaSqrts[gammaJ], fmac(-COMPLEX_ONE_FOURTH, bettaSqrts[bettaI], alpha_complex)));

    if (numberOfRoots == 4)
    {
        roots[0] = w0.real();
        roots[1] = w1.real();
        roots[2] = w2.real();
        roots[3] = w3.real();
    }
    else if (numberOfRoots == 2) {
        if (!isComplex(w0)) roots[0] = w0.real();
        if (!isComplex(w1)) roots[1] = w1.real();
        if (!isComplex(w2)) roots[2] = w2.real();
        if (!isComplex(w3)) roots[3] = w3.real();
    }
    return numberOfRoots;
}

/*
    Имплементация метода решения уравнения четвертой степени —A_Note_on_the_Solution_of_Quartic_Equations-Salzer-1960
    Информация о методе — https://www.ams.org/journals/mcom/1960-14-071/S0025-5718-1960-0117882-6/S0025-5718-1960-0117882-6.pdf
    Работу выполнил — Погосов Даниэль (https://github.com/DarklleS)
*/
template<typename fp_t>
unsigned int salzer(fp_t N, fp_t A, fp_t B, fp_t C, fp_t D, vector<fp_t>& roots)
{
    // Нормировка исходных коэффициентов, где N - старший коэффициент уравнения
    if (isZero(N) || isinf(A /= N))
        return solveCubic(A, B, C, D, roots);
    if (isinf(B /= N))
        return 0;
    if (isinf(C /= N))
        return 0;
    if (isinf(D /= N))
        return 0;

    // Количество вещественных корней
    unsigned int numberOfRoots = 0;

    // Вычисление коэффициентов кубического уравнения
    fp_t a = -B;
    fp_t b = fms(A, C, static_cast<fp_t>(4), D);  // A * C - 4 * D
    fp_t c = fms(D, fms(static_cast<fp_t>(4), B, A, A), C, C); // D * (4 * B - A^2) -C ^ 2
    vector<fp_t> cubicRoots(3);

    // Решение кубического уравнения
    solveCubic(static_cast<fp_t>(1.L), a, b, c, cubicRoots);

    // Вещественный корень кубического уравнения
    fp_t x = cubicRoots[0];

    // Вычисление начальных расчетных коэффициентов
    fp_t m;
    fp_t n;
    fp_t mm = fma(static_cast<fp_t>(0.25L) * A, A, -B) + x; // A^2 / 4 - B + c

    if (mm > 0)
    {
        m = sqrt(mm);
        n = static_cast<fp_t>(0.25L) * fms(A, x, static_cast<fp_t>(2), C) / m; // (A * x - 2 * C) / (2 * m)
    }
    else if (isZero(mm))
    {
        m = 0;
        n = sqrt(fma(static_cast<fp_t>(0.25L) * x, x, -D));  // sqrt(x^2 / 4 - D)
    }
    else // m - комплексное, следовательно уравнение не будет иметь вещественных корней
    {
        return 0;
    }

    // Вычисление расчетных коэффициентов
    fp_t alpha = fma(static_cast<fp_t>(0.5L) * A, A, -x) - B; // A^2 / 2 - x - B
    fp_t beta = fms(static_cast<fp_t>(4), n, A, m); // 4 * n - A * m
    fp_t gamma = alpha + beta;
    fp_t delta = alpha - beta;

    if (gamma >= 0) // Если gamma >= 0, то уравнение имеет два либо больше вещественных корней
    {
        gamma = sqrt(gamma);

        roots[0] = fma(static_cast<fp_t>(0.5L), gamma, fms(static_cast<fp_t>(0.5L), m, static_cast<fp_t>(0.25L), A)); // gamma / 2 + (m / 2 - A / 4)
        roots[1] = fma(static_cast<fp_t>(-0.5L), gamma, fms(static_cast<fp_t>(0.5L), m, static_cast<fp_t>(0.25L), A)); // -gamma / 2 + (m / 2 - A / 4)

        numberOfRoots += 2;
    }
    if (delta >= 0) // Если delta >= 0, то уравнение имеет два либо больше вещественных корней
    {
        delta = sqrt(delta);

        roots[numberOfRoots] = fma(static_cast<fp_t>(0.5L), delta, fms(static_cast<fp_t>(-0.5L), m, static_cast<fp_t>(0.25L), A)); // delta / 2 + (-m / 2 - A / 4)
        roots[numberOfRoots + 1] = fma(static_cast<fp_t>(-0.5L), delta, fms(static_cast<fp_t>(-0.5L), m, static_cast<fp_t>(0.25L), A)); // delta / 2 + (-m / 2 - A / 4)

        numberOfRoots += 2;
    }

    return numberOfRoots;
}

template unsigned int solveLinear(float a, float b, vector<float>& roots);
template unsigned int solveLinear(double a, double b, vector<double>& roots);
template unsigned int solveLinear(long double a, long double b, vector<long double>& roots);

template unsigned int solveLinear(float a, float b, vector<complex<float>>& roots);
template unsigned int solveLinear(double a, double b, vector<complex<double>>& roots);
template unsigned int solveLinear(long double a, long double b, vector<complex<long double>>& roots);

template unsigned int solveLinear(complex<float> a, complex<float> b, vector<complex<float>>& roots);
template unsigned int solveLinear(complex<double> a, complex<double> b, vector<complex<double>>& roots);
template unsigned int solveLinear(complex<long double> a, complex<long double> b, vector<complex<long double>>& roots);

template unsigned int solveQuadratic(float n, float a, float b, vector<float>& roots);
template unsigned int solveQuadratic(double n, double a, double b, vector<double>& roots);
template unsigned int solveQuadratic(long double n, long double a, long double b, vector<long double>& roots);
 
template unsigned int solveQuadratic(float n, float a, float b, vector<complex<float>>& roots);
template unsigned int solveQuadratic(double n, double a, double b, vector<complex<double>>& roots);
template unsigned int solveQuadratic(long double n, long double a, long double b, vector<complex<long double>>& roots);

template unsigned int solveQuadratic(complex<float> n, complex<float> a, complex<float> b, vector<complex<float>>& roots);
template unsigned int solveQuadratic(complex<double> n, complex<double> a, complex<double> b, vector<complex<double>>& roots);
template unsigned int solveQuadratic(complex<long double> n, complex<long double> a, complex<long double> b, vector<complex<long double>>& roots);

template unsigned int solveCubic(float n, float a, float b, float c, vector<float>& roots);
template unsigned int solveCubic(double n, double a, double b, double c, vector<double>& roots);
template unsigned int solveCubic(long double n, long double a, long double b, long double c, vector<long double>& roots);

template unsigned int solveCubic(float n, float a, float b, float c, vector<complex<float>>& roots);
template unsigned int solveCubic(double n, double a, double b, double c, vector<complex<double>>& roots);
template unsigned int solveCubic(long double n, long double a, long double b, long double c, vector<complex<long double>>& roots);

template unsigned int ferrari(float n, float a, float b, float c, float d, vector<float>& roots);
template unsigned int ferrari(double n, double a, double b, double c, double d, vector<double>& roots);
template unsigned int ferrari(long double n, long double a, long double b, long double c, long double d, vector<long double>& roots);

template unsigned int descartes(float n, float a, float b, float c, float d, vector<float>& roots);
template unsigned int descartes(double n, double a, double b, double c, double d, vector<double>& roots);
template unsigned int descartes(long double n, long double a, long double b, long double c, long double d, vector<long double>& roots);

template unsigned int nbs(float n, float a, float b, float c, float d, vector<float>& roots);
template unsigned int nbs(double n, double a, double b, double c, double d, vector<double>& roots);
template unsigned int nbs(long double n, long double a, long double b, long double c, long double d, vector<long double>& roots);

template unsigned int euler(float n, float a, float b, float c, float d, vector<float>& roots);
template unsigned int euler(double n, double a, double b, double c, double d, vector<double>& roots);
template unsigned int euler(long double n, long double a, long double b, long double c, long double d, vector<long double>& roots);

template unsigned int vanDerWaerden(float n, float a, float b, float c, float d, vector<float>& roots);
template unsigned int vanDerWaerden(double n, double a, double b, double c, double d, vector<double>& roots);
template unsigned int vanDerWaerden(long double n, long double a, long double b, long double c, long double d, vector<long double>& roots);

template unsigned int tschirnhaus(float n, float a, float b, float c, float d, vector<float>& roots);
template unsigned int tschirnhaus(double n, double a, double b, double c, double d, vector<double>& roots);
template unsigned int tschirnhaus(long double n, long double a, long double b, long double c, long double d, vector<long double>& roots);

template unsigned int ferrari_complex(float n, float a, float b, float c, float d, vector<complex<float>>& roots);
template unsigned int ferrari_complex(double n, double a, double b, double c, double d, vector<complex<double>>& roots);
template unsigned int ferrari_complex(long double n, long double a, long double b, long double c, long double d, vector<complex<long double>>& roots);

template unsigned int fqs(float n, float a, float b, float c, float d, vector<float>& roots);
template unsigned int fqs(double n, double a, double b, double c, double d, vector<double>& roots);
template unsigned int fqs(long double n, long double a, long double b, long double c, long double d, vector<long double>& roots);

template unsigned int merriman(float n, float a, float b, float c, float d, vector<float>& roots);
template unsigned int merriman(double n, double a, double b, double c, double d, vector<double>& roots);
template unsigned int merriman(long double n, long double a, long double b, long double c, long double d, vector<long double>& roots);

template unsigned int findRoots(float a3, float a2, float p, vector<float>& roots);
template unsigned int findRoots(double a3, double a2, double p, vector<double>& roots);
template unsigned int findRoots(long double a3, long double a2, long double p, vector<long double>& roots);

template unsigned int squire(float n, float a, float b, float c, float d, vector<float>& roots);
template unsigned int squire(double n, double a, double b, double c, double d, vector<double>& roots);
template unsigned int squire(long double n, long double a, long double b, long double c, long double d, vector<long double>& roots);

template unsigned int ungar(float n, float a, float b, float c, float d, vector<float>& roots);
template unsigned int ungar(double n, double a, double b, double c, double d, vector<double>& roots);
template unsigned int ungar(long double n, long double a, long double b, long double c, long double d, vector<long double>& roots);

template unsigned int salzer(float N, float A, float B, float C, float D, vector<float>& roots);
template unsigned int salzer(double N, double A, double B, double C, double D, vector<double>& roots);
template unsigned int salzer(long double N, long double A, long double B, long double C, long double D, vector<long double>& roots);