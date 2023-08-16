#include "../inc/excerpt.h"
#include "../inc/polynomials.h"

// Функция проверки элемента на ноль (x == 0)
template<typename fp_t>
inline bool isZero(static const fp_t& x)
{
    return FP_ZERO == fpclassify(x);
}

// Функция проверки комплексного элемента на ноль (x == 0)
template<typename fp_t>
inline bool isZero(static const complex<fp_t>& x)
{
    return (FP_ZERO == fpclassify(x.real())) && (FP_ZERO == fpclassify(x.imag()));
}

// Комплексная версия метода isinf (проверка на бесконечность)
template<typename fp_t>
inline bool isinfc(static const complex<fp_t>& x)
{
    return (FP_INFINITE == fpclassify(x.real())) || (FP_INFINITE == fpclassify(x.imag()));
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

// Комплексная версия fma
template<typename fp_t>
inline complex<fp_t> fmac(complex<fp_t> x, complex<fp_t> y, complex<fp_t> z)
{
    fp_t r = fma(-x.imag(), y.imag(), fma(x.real(), y.real(), z.real()));
    fp_t i = fma(y.real(), x.imag(), fma(x.real(), y.imag(), z.imag()));

    return complex<fp_t>(r, i);
}

// Комплексная версия fms
template<typename fp_t>
inline complex<fp_t> fmsc(complex<fp_t> a, complex<fp_t> b, complex<fp_t> c, complex<fp_t> d)
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

// Сигнатурная функция (для комплексной величины)
template <typename fp_t>
inline int sgn(static const complex<fp_t>& x)
{
    return sgn(x.real()) != 0 ? sgn(x.real()) : sgn(x.imag());
}

// Функция определяющая является ли число вещественным, путем определения того, насколько мала мнимая часть
template<typename fp_t>
inline complex<fp_t> epsilonComplex(complex<fp_t> x)
{
    return abs(x) * numeric_limits<fp_t>::epsilon() > abs(x.imag()) ? complex<fp_t>(x.real(), 0) : x;
}

// Проверка значения на комплексность 
template<typename fp_t>
inline bool isComplex(static const complex<fp_t>& x)
{
    return abs(x) * numeric_limits<fp_t>::epsilon() <= abs(x.imag());
}

// Дискриминант 
template<typename fp_t>
fp_t dscrmt(fp_t A, fp_t B, fp_t C)
{
    fp_t p = B * B;
    fp_t q = A * C;

    fp_t dp = std::fma(B, B, -p);
    fp_t dq = std::fma(A, C, -q);
    fp_t d = (p - q) + (dp - dq);

    return d;
}

// Дискриминант (для комплексной величины)
template<typename fp_t>
complex<fp_t> dscrmt(complex<fp_t> A, complex<fp_t> B, complex<fp_t> C)
{
    complex<fp_t> p = B * B;
    complex<fp_t> q = A * C;

    complex<fp_t> dp = fmac(B, B, -p);
    complex<fp_t> dq = fmac(A, C, -q);
    complex<fp_t> d = (p - q) + (dp - dq);

    return d;
}

// Кубический корень от комплексного числа
template<typename fp_t>
vector<complex<fp_t>> cubeRoot(complex<fp_t> number)
{
    vector<complex<fp_t>> result;
    static const fp_t _2PI = static_cast<fp_t>(numbers::pi) * static_cast<fp_t>(2.0L);
    static const fp_t ONE_THIRDS = static_cast<fp_t>(1.0L / 3.0L);

    fp_t argValue = arg(number);
    fp_t modulus = abs(number);

    for (int k = 0; k < 3; ++k) {

        fp_t realOfCube = cos(fma(_2PI, k, argValue) * ONE_THIRDS);
        fp_t imagOfCube = sin(fma(_2PI, k, argValue) * ONE_THIRDS);

        if (abs(imagOfCube) < std::numeric_limits<fp_t>::epsilon()) { imagOfCube = 0; }

        complex<fp_t> root_k = pow(modulus, ONE_THIRDS) * complex<fp_t>(realOfCube, imagOfCube);

        result.push_back(root_k);
    }

    return result;
}

// Функция отбора вещественных корней
template<typename fp_t>
inline void selectRealRoots(vector<complex<fp_t>> foundRoots, unsigned int& numberOfRealRoots, vector<fp_t>& roots)
{
    for (int i = 0; i < foundRoots.size(); ++i)
    {
        if (!isComplex(foundRoots[i]))
        {
            roots[numberOfRealRoots] = foundRoots[i].real();
            ++numberOfRealRoots;
        }
    }

    return;
}

//Быстрый алгоритм наименьших квадратов для вычисления коэффициентов {γ0, δ0}
//для заданного набора{ α0, β0 }
template <typename fp_t>
unsigned int least_squares(fp_t a, fp_t b, fp_t c, fp_t d
    , fp_t alpha0, fp_t beta0, vector<fp_t>& coeff) {

    fp_t F1, F2;
    F1 = fma(beta0, beta0, fma(alpha0, alpha0, static_cast<fp_t>(1.0L))); //F1 = 1 + pow(alpha0, 2) + pow(beta0, 2);   
    F2 = fma(alpha0, beta0, alpha0);        //F2 = alpha0 * (1 + beta0);

    fp_t c1, c2;
    c1 = fma(alpha0, b - beta0, -alpha0) + fma(beta0, c, a);     //c1 = a - alpha0 + alpha0 * (b - beta0) + beta0 * c;
    c2 = fma(alpha0, c, b) + fma(beta0, d, -beta0);     //c2 = b - beta0 + alpha0 * c + beta0 * d;

    fp_t L1, L2, L3;
    L1 = sqrt(F1);
    L3 = F2 / L1;
    L2 = sqrt(fma(-F2 / F1, F2, F1));       //L2 = sqrt(F1 - (F2 / F1) * F2);

    fp_t y1, y2;
    y1 = c1 / L1;
    y2 = fma(-y1, L3, c2) / L2;     //y2 = (c2 - y1 * L3) / L2;

    //Находим коэффициенты γ0, δ0
    fp_t delta0, gamma0;
    delta0 = y2 / L2;
    gamma0 = fma(-delta0, L3, y1) / L1;     //gamma0 = (y1 - delta0 * L3) / L1;

    coeff[0] = delta0;
    coeff[1] = gamma0;

    return 0;
}

//Итеративный алгоритм уточнения коэффициентов {α, β, γ, δ} из заданного начального
//приближения{ α0, β0, γ0, δ0 }.
template <typename fp_t>
unsigned int backward_optimizer(fp_t a, fp_t b, fp_t c, fp_t d
    , fp_t alpha0, fp_t beta0, fp_t gamma0, fp_t delta0, unsigned maxIters, vector<fp_t>& coeff) {

    fp_t x1, x2, x3, x4;
    fp_t y1, y2, y3, y4;

    fp_t U23, U33, U44, L43;

    //Начальное приближение коэффициентов.
    fp_t alpha = alpha0, beta = beta0, gamma = gamma0, delta = delta0;

    //Инициализация величин ошибок на коэффициентах.
    fp_t e1 = a - alpha - gamma;
    fp_t e2 = b - fma(alpha, gamma, beta + delta); //b - beta - alpha * gamma - delta;
    fp_t e3 = fma(-alpha, delta, fma(-beta, gamma, c)); //c - beta * gamma - alpha * delta;
    fp_t e4 = fma(-beta, delta, d); //d - beta * delta;

    //Величина суммарной ошибки
    fp_t eps = 0;
    //Счетчик итераций
    unsigned iters = 0;

    //Очередь для хранения значений eps для 4-х предыдущих итераций
    queue <fp_t> q; q.push(0); q.push(0); q.push(0); q.push(0);
    fp_t tmp1, tmp2, tmp3, tmp4;

    for (int i = 0; i < maxIters; i++) {

        //LU-факторизация.
        U23 = alpha - gamma;
        U33 = -fma(gamma, U23, delta - beta);   //beta - delta - gamma * U23;

        L43 = -delta * U23 / U33;

        U44 = -fma(L43, U23, delta - beta); //beta - delta - L43 * U23;

        //Вычисление компонент {x1, x2, x3, x3} вспомогательного вектора x. Lx = e → x.
        x1 = e1;
        x2 = fma(-gamma, x1, e2);   //e2 - gamma * x1;
        x3 = fma(-gamma, x2, fma(-delta, x1, e3));  //e3 - delta * x1 - gamma * x2;
        x4 = fma(-L43, x3, fma(-delta, x2, e4));    //e4 - delta * x2 - L43 * x3;

        //Определения компонент {y1, y2, y3, y3} вектора обновления y. Uy = x → y.
        y4 = x4 / U44;
        y3 = fma(-U23, y4, x3) / U33; //(x3 - U23 * y4) / U33;
        y2 = fma(-U23, y3, x2 - y4);    //x2 - U23 * y3 - y4;
        y1 = x1 - y3;

        //Обновление коэффициентов.
        alpha = alpha + y1;
        beta = beta + y2;
        gamma = gamma + y3;
        delta = delta + y4;

        //Вычисление величин ошибок на коэффициентах.
        e1 = a - alpha - gamma;
        e2 = b - fma(alpha, gamma, beta + delta);  //b - beta - alpha * gamma - delta;
        e3 = fma(-alpha, delta, fma(-beta, gamma, c));  //c - beta * gamma - alpha * delta;
        e4 = fma(-beta, delta, d);  //d - beta * delta;

        iters = i;

        //Величина суммарной ошибки
        eps = abs(e1) + abs(e2) + abs(e3) + abs(e4);

        tmp1 = q.front(); q.pop();
        tmp2 = q.front(); q.pop();
        tmp3 = q.front(); q.pop();
        tmp4 = q.front(); q.pop();

        //Условие сходимости алгоритма.
        if (eps < 1e-07 or isEqual(eps, tmp1) or isEqual(eps, tmp2) or isEqual(eps, tmp3) or isEqual(eps, tmp4)) break;
        //eps < 1e-07 подобрано экспериментально

        //Обновление элементов очереди.
        q.push(tmp2); q.push(tmp3); q.push(tmp4); q.push(eps);
    }

    coeff[0] = alpha;
    coeff[1] = beta;
    coeff[2] = gamma;
    coeff[3] = delta;
    coeff[4] = eps;

    return iters;
}

//Функция вычисления коэффициентов квадратных уравнений
template<typename fp_t>
unsigned int FQScoeffs(fp_t n, fp_t a, fp_t b, fp_t c, fp_t d, vector<fp_t>& coeffs)
{
    int maxIters = 10000;
    vector<complex<fp_t>> x(4);
    ferrari_complex<fp_t>(n, a, b, c, d, x); //Получаем приближение корней методом Феррари

    //Сортируем корни в порядке убывания
    complex<fp_t> tmp = 0;
    for (int i = 0; i < 4; i++)
        for (int j = (4 - 1); j >= (i + 1); j--)
            if (abs(x[j]) > abs(x[j - 1])) {
                tmp = x[j]; x[j] = x[j - 1]; x[j - 1] = tmp;
            }

    //Вычисляем начальное приближение коэффициентов {α01, β01} и {α02, β02}
    fp_t alpha01, beta01, gamma01, delta01;
    fp_t alpha02, beta02, gamma02, delta02;

    alpha01 = -(x[0] + x[1]).real(); beta01 = (x[0] * x[1]).real();
    alpha02 = -(x[1] + x[2]).real(); beta02 = (x[1] * x[2]).real();

    //С помощью алгоритма наименьших квадратов для каждого из набора найдем коэффициенты { γ01, δ01 } и { γ02, δ02 } соответственно.
    vector<fp_t> coeff1(2);
    vector<fp_t> coeff2(2);

    least_squares(a, b, c, d, alpha01, beta01, coeff1);
    least_squares(a, b, c, d, alpha02, beta02, coeff2);

    delta01 = coeff1[0];
    gamma01 = coeff1[1];

    delta02 = coeff2[0];
    gamma02 = coeff2[1];

    //Найдем коэффициенты для каждого из двух начальных приближений
    vector<fp_t> res1(5);
    vector<fp_t> res2(5);

    unsigned itersCount1 = backward_optimizer(a, b, c, d, alpha01, beta01, gamma01, delta01, maxIters, res1);
    unsigned itersCount2 = backward_optimizer(a, b, c, d, alpha02, beta02, gamma02, delta02, maxIters, res2);

    fp_t eps1 = res1[4];
    fp_t eps2 = res2[4];

    //Определим какой набор коэффициентов правильный
    fp_t alpha, beta, gamma, delta;

    if (itersCount1 < itersCount2) {
        alpha = res1[0]; beta = res1[1]; gamma = res1[2]; delta = res1[3];
    }
    else if (itersCount1 > itersCount2) {
        alpha = res2[0]; beta = res2[1]; gamma = res2[2]; delta = res2[3];
    }
    else if (eps1 < eps2) {
        alpha = res1[0]; beta = res1[1]; gamma = res1[2]; delta = res1[3];
    }
    else {
        alpha = res2[0]; beta = res2[1]; gamma = res2[2]; delta = res2[3];
    }

    coeffs[0] = alpha;
    coeffs[1] = beta;
    coeffs[2] = gamma;
    coeffs[3] = delta;
}

// Итеративная функция нахождения корней уравнения методом бисекции
template<typename fp_t>
fp_t refine(fp_t left, fp_t right, fp_t* coef, fp_t(*fun)(fp_t, fp_t*))
{
    fp_t xl = left;
    fp_t xr = right;

    fp_t eps3 = eps_temp<fp_t>;
    fp_t p = fma((xr - xl), static_cast<fp_t>(0.5L), xl);

    while (xr - xl > eps3)
    {
        p = fma((xr - xl), static_cast<fp_t>(0.5L), xl);
        fp_t fp = fun(p, coef);
        if (p == xl || p == xr || fp == 0) return p;
        if (fp > 0) xr = p;
        else xl = p;
    }
    return p;
}

/// Метод для нахождения индексов у массивов кубических корней от `betta0, gamma0`, таких что выполнилось условие 4.10
/// - Parameters:
///   - betta0CubeRoots: Кубические корни betta0
///   - gamma0CubeRoots: Кубические корни gamma0
///   - Q: Расчетный коэффициент Q
///   - beta0_I: Индекс корня от `beta0`
///   - gamma0_J: Индекс корня от `gamma0`
template<typename fp_t>
void check4_10Condition(vector<complex<fp_t>> betta0CubeRoots,
    vector<complex<fp_t>> gamma0CubeRoots,
    fp_t Q, int& beta0_I, int& gamma0_J) {
    fp_t min_r = numeric_limits<fp_t>::infinity();
    fp_t min_rr = numeric_limits<fp_t>::infinity();

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            complex<fp_t> helper = betta0CubeRoots[i] * gamma0CubeRoots[j];
            if (abs(helper.imag()) <= min_rr)
            {
                min_rr = abs(helper.imag());
                if (abs(helper.real() - Q) <= min_r)
                {
                    min_r = abs(helper.real() - Q);
                    beta0_I = i;
                    gamma0_J = j;
                    break;
                }
            }
        }
    }
}

/// Метод для нахождения индексов у массивов корней от `betta, gamma delta`, таких что выполнилось условие 4.11
/// - Parameters:
///   - bettaSqrts: Корни от `betta`
///   - gammaSqrts: Корни от `gamma`
///   - deltaSqrts: Корни от `delta`
///   - P: Расчетный коэффициент P
///   - bettaI: Индекс корня от `betta`
///   - gammaJ: Индекс корня от `gamma`
///   - deltaK: Индекс корня от `delta`
template<typename fp_t>
void check4_11Condition(vector<complex<fp_t>> bettaSqrts,
    vector<complex<fp_t>> gammaSqrts,
    vector<complex<fp_t>> deltaSqrts,
    fp_t P, int& bettaI, int& gammaJ,
    int& deltaK)
{
    fp_t min_r1 = numeric_limits<fp_t>::infinity();
    fp_t min_r2 = numeric_limits<fp_t>::infinity();

    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            for (int k = 0; k < 2; ++k)
            {
                complex<fp_t> bettaGammaDeltaSqrts = complex<fp_t>(bettaSqrts[i] * gammaSqrts[j] * deltaSqrts[k]);
                if (abs(bettaGammaDeltaSqrts.imag()) <= min_r2 && abs(bettaGammaDeltaSqrts.real() + P) <= min_r1) {
                    min_r2 = abs(bettaGammaDeltaSqrts.imag());
                    min_r1 = abs(bettaGammaDeltaSqrts.real() + P);
                    bettaI = i;
                    gammaJ = j;
                    deltaK = k;

                    break;
                }
            }
        }
    }
}

// Метод проведения оценки точности результатов работы метода уравнения четвертой степени
template<typename fp_t>
void testQuarticPolynomial(
    int testCount, // Кол-во экспериментов
    long double maxDistance, // Максимальное расстояние между корнями
    int P1, // Кол-во кластеризированных корней
    int P2 // Кол-во кратных корней
)
{
    unsigned P = 4; // Степень исходного полинома
    fp_t low = -1, high = 1; // Интервал на котором заданы корни полинома
    fp_t absMaxError, relMaxError; // Абсолютная и относительная погрешность по итогам пройденного теста
    fp_t absMaxErrorTotal = -1, relMaxErrorTotal = -1; // Итоговая максимальная абсолютная и относительная погрешность по итогам всех тестов
    long double absErrorAvg = 0, relErrorAvg = 0; // Средняя абсолютная и относительная погрешность по итогам всех тестов
    unsigned numberOfFoundRoots; // Количество найденных корней
    unsigned cantFind = 0; // Счетчик количества ситуаций, когда методу не удалось найти корни (numberOfFoundRoots == 0)
    std::vector<fp_t> coefficients(P + 1); // Вектор коэффициентов полинома
    unsigned count = 0; // Счетчик количества ситуаций, когда относительная погрешность больше определенного числа (relMaxError > n)
    int countExcessRoots = 0; // Счетчик лишних корней
    int countLostRoots = 0; // Счетчик потерянных корней
    std::vector<fp_t> errors; // Массив величин наибольшей отн. погрешности среди четырех корней исходного полинома

    for (size_t i = 0; i < testCount; ++i)
    {
        std::vector<fp_t> foundRoots(P);
        std::vector<fp_t> trueRoots(P);
        int excessRoots = 0;
        int lostRoots = 0;

        generate_polynomial<fp_t>(P, 0, P1, P2, static_cast<fp_t>(maxDistance), low, high, trueRoots, coefficients);

        numberOfFoundRoots = ferrari<fp_t>(coefficients[4], coefficients[3], coefficients[2], coefficients[1], coefficients[0], foundRoots);

        if (numberOfFoundRoots > 0)
        {
            compare_roots<fp_t>(numberOfFoundRoots, P, foundRoots, trueRoots, absMaxError, relMaxError, excessRoots, lostRoots, errors);

            absMaxErrorTotal = absMaxError > absMaxErrorTotal ? absMaxError : absMaxErrorTotal;
            absErrorAvg += absMaxError;

            relMaxErrorTotal = relMaxError > relMaxErrorTotal ? relMaxError : relMaxErrorTotal;
            relErrorAvg += relMaxError;

            countExcessRoots += excessRoots;
            countLostRoots += lostRoots;

            count += relMaxError > 0.95 ? 1 : 0;
        }
        else
        {
            countLostRoots += 4;
            cantFind += 1;
        }

    }

    absErrorAvg /= (testCount - cantFind);
    relErrorAvg /= (testCount - cantFind);
    std::sort(errors.begin(), errors.end());
    std::string label("");
    long double quantile = 0;

    if (P1 == 0 && P2 == 0)
    {
        label = std::string("NON-CLUSTER ");
        quantile = 0.01;
    }
    else if (P1 == 0)
    {
        label = std::string("MULTIPLE ");
        quantile = 0.1;
    }
    else if (P2 == 0)
    {
        label = std::string("CLUSTER ");
        quantile = 0.1;
    }

    fp_t quantileMaxErr = errors[errors.size() - errors.size() * quantile];

    if (PRINT)
    {
        std::cout << label << "QUARTIC TEST RESULTS" << std::endl;
        std::cout << "========================================" << std::endl;
        std::cout << "Max distance: " << maxDistance << std::endl;
        std::cout << "Total count of tests: " << testCount << std::endl;
        std::cout << "Couldn't find roots: " << cantFind << " times " << std::endl;
        std::cout << "----------------------------------------" << std::endl;
        std::cout << "Average absolute error: " << absErrorAvg << std::endl;
        std::cout << "Total maximum absolute error: " << absMaxErrorTotal << std::endl;
        std::cout << "Average relative error: " << relErrorAvg << std::endl;
        std::cout << "Total maximum relative error: " << relMaxErrorTotal << std::endl;
        std::cout << "----------------------------------------" << std::endl;
        std::cout << "Total count of lost roots: " << countLostRoots << std::endl;
        std::cout << "Total count of excess roots: " << countExcessRoots << std::endl;
        std::cout << "----------------------------------------" << std::endl;
        std::cout << "relMaxError > 0.95: " << count << " times" << std::endl;
        std::cout << "Quantile max error: " << quantileMaxErr << std::endl;
        std::cout << "========================================" << std::endl << std::endl;
    }
}

template bool isZero(static const float& x);
template bool isZero(static const double& x);
template bool isZero(static const long double& x);

template bool isZero(static const complex<float>&x);
template bool isZero(static const complex<double>&x);
template bool isZero(static const complex<long double>&x);

template bool isinfc(static const complex<float>& x);
template bool isinfc(static const complex<double>& x);
template bool isinfc(static const complex<long double>& x);

template bool isEqual(const float& A, const float& B);
template bool isEqual(const double& A, const double& B);
template bool isEqual(const long double& A, const long double& B);

template float fms(float a, float b, float c, float d);
template double fms(double a, double b, double c, double d);
template long double fms(long double a, long double b, long double c, long double d);

template complex<float> fmac(complex<float> x, complex<float> y, complex<float> z);
template complex<double> fmac(complex<double> x, complex<double> y, complex<double> z);
template complex<long double> fmac(complex<long double> x, complex<long double> y, complex<long double> z);

template complex<float> fmsc(complex<float> a, complex<float> b, complex<float> c, complex<float> d);
template complex<double> fmsc(complex<double> a, complex<double> b, complex<double> c, complex<double> d);
template complex<long double> fmsc(complex<long double> a, complex<long double> b, complex<long double> c, complex<long double> d);

template int sgn(static const float& x);
template int sgn(static const double& x);
template int sgn(static const long double& x);

template int sgn(static const complex<float>& x);
template int sgn(static const complex<double>& x);
template int sgn(static const complex<long double>& x);

template complex<float> epsilonComplex(complex<float> x);
template complex<double> epsilonComplex(complex<double> x);
template complex<long double> epsilonComplex(complex<long double> x);

template bool isComplex(static const complex<float>& x);
template bool isComplex(static const complex<double>& x);
template bool isComplex(static const complex<long double>& x);
 
template float dscrmt(float A, float B, float C);
template double dscrmt(double A, double B, double C);
template long double dscrmt(long double A, long double B, long double C);

template complex<float> dscrmt(complex<float> A, complex<float> B, complex<float> C);
template complex<double> dscrmt(complex<double> A, complex<double> B, complex<double> C);
template complex<long double> dscrmt(complex<long double> A, complex<long double> B, complex<long double> C);

template vector<complex<float>> cubeRoot(complex<float> number);
template vector<complex<double>> cubeRoot(complex<double> number);
template vector<complex<long double>> cubeRoot(complex<long double> number);

template void selectRealRoots(vector<complex<float>> foundRoots, unsigned int& numberOfRealRoots, vector<float>& roots);
template void selectRealRoots(vector<complex<double>> foundRoots, unsigned int& numberOfRealRoots, vector<double>& roots);
template void selectRealRoots(vector<complex<long double>> foundRoots, unsigned int& numberOfRealRoots, vector<long double>& roots);

template unsigned int least_squares(float a, float b, float c, float d, float alpha0, float beta0, vector<float>& coeff);
template unsigned int least_squares(double a, double b, double c, double d, double alpha0, double beta0, vector<double>& coeff);
template unsigned int least_squares(long double a, long double b, long double c, long double d, long double alpha0, long double beta0, vector<long double>& coeff);

template unsigned int backward_optimizer(float a, float b, float c, float d, float alpha0, float beta0, float gamma0, float delta0, unsigned maxIters, vector<float>& coeff);
template unsigned int backward_optimizer(double a, double b, double c, double d, double alpha0, double beta0, double gamma0, double delta0, unsigned maxIters, vector<double>& coeff);
template unsigned int backward_optimizer(long double a, long double b, long double c, long double d, long double alpha0, long double beta0, long double gamma0, long double delta0, unsigned maxIters, vector<long double>& coeff);

template unsigned int FQScoeffs(float n, float a, float b, float c, float d, vector<float>& coeffs);
template unsigned int FQScoeffs(double n, double a, double b, double c, double d, vector<double>& coeffs);
template unsigned int FQScoeffs(long double n, long double a, long double b, long double c, long double d, vector<long double>& coeffs);

template float refine(float left, float right, float* coef, float(*fun)(float, float*));
template double refine(double left, double right, double* coef, double(*fun)(double, double*));
template long double refine(long double left, long double right, long double* coef, long double(*fun)(long double, long double*));

template void check4_10Condition(vector<complex<float>> betta0CubeRoots, vector<complex<float>> gamma0CubeRoots, float Q, int& beta0_I, int& gamma0_J);
template void check4_10Condition(vector<complex<double>> betta0CubeRoots, vector<complex<double>> gamma0CubeRoots, double Q, int& beta0_I, int& gamma0_J);
template void check4_10Condition(vector<complex<long double>> betta0CubeRoots, vector<complex<long double>> gamma0CubeRoots, long double Q, int& beta0_I, int& gamma0_J);

template void check4_11Condition(vector<complex<float>> bettaSqrts, vector<complex<float>> gammaSqrts, vector<complex<float>> deltaSqrts, float P, int& bettaI, int& gammaJ, int& deltaK);
template void check4_11Condition(vector<complex<double>> bettaSqrts, vector<complex<double>> gammaSqrts, vector<complex<double>> deltaSqrts, double P, int& bettaI, int& gammaJ, int& deltaK);
template void check4_11Condition(vector<complex<long double>> bettaSqrts, vector<complex<long double>> gammaSqrts, vector<complex<long double>> deltaSqrts, long double P, int& bettaI, int& gammaJ, int& deltaK);

template void testQuarticPolynomial<float>(int testCount, long double maxDistance, int P1, int P2);
template void testQuarticPolynomial<double>(int testCount, long double maxDistance, int P1, int P2);
template void testQuarticPolynomial<long double>(int testCount, long double maxDistance, int P1, int P2);
