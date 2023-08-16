#ifndef MODULES_H
#define MODULES_H

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

// Функция проверки элемента на ноль (x == 0)
template<typename fp_t>
inline bool isZero(static const fp_t& x);

// Функция проверки комплексного элемента на ноль (x == 0)
template<typename fp_t>
inline bool isZero(static const complex<fp_t>& x);

// Комплексная версия метода isinf (проверка на бесконечность)
template<typename fp_t>
inline bool isinfc(static const complex<fp_t>& x);

// Функция проверки равенства двух переменных (A == B)
template<typename fp_t>
inline bool isEqual(const fp_t& A, const fp_t& B);

// Функция вычисления разности произведений с оптимизацией погрешности (a * b - c * d)
template<typename fp_t>
inline fp_t fms(fp_t a, fp_t b, fp_t c, fp_t d);

// Комплексная версия fma
template<typename fp_t>
inline complex<fp_t> fmac(complex<fp_t> x, complex<fp_t> y, complex<fp_t> z);

// Комплексная версия fms
template<typename fp_t>
inline complex<fp_t> fmsc(complex<fp_t> a, complex<fp_t> b, complex<fp_t> c, complex<fp_t> d);

// Сигнатурная функция
template <typename fp_t>
inline int sgn(static const fp_t& x);

// Сигнатурная функция (для комплексной величины)
template <typename fp_t>
inline int sgn(static const complex<fp_t>& x);

// Функция определяющая является ли число вещественным, путем определения того, насколько мала мнимая часть
template<typename fp_t>
inline complex<fp_t> epsilonComplex(complex<fp_t> x);

// Проверка значения на комплексность 
template<typename fp_t>
inline bool isComplex(static const complex<fp_t>& x);

// Дискриминант 
template<typename fp_t>
fp_t dscrmt(fp_t A, fp_t B, fp_t C);

// Дискриминант (для комплексных величин)
template<typename fp_t>
complex<fp_t> dscrmt(complex<fp_t> A, complex<fp_t> B, complex<fp_t> C);

// Кубический корень от комплексного числа
template<typename fp_t>
vector<complex<fp_t>> cubeRoot(complex<fp_t> number);

// Функция отбора вещественных корней
template<typename fp_t>
inline void selectRealRoots(vector<complex<fp_t>> foundRoots, unsigned int& numberOfRealRoots, vector<fp_t>& roots);

//Быстрый алгоритм наименьших квадратов для вычисления коэффициентов {γ0, δ0}
//для заданного набора{ α0, β0 }
template <typename fp_t>
unsigned int least_squares(fp_t a, fp_t b, fp_t c, fp_t d, fp_t alpha0, fp_t beta0, vector<fp_t>& coeff);

//Итеративный алгоритм уточнения коэффициентов {α, β, γ, δ} из заданного начального
//приближения{ α0, β0, γ0, δ0 }.
template <typename fp_t>
unsigned int backward_optimizer(fp_t a, fp_t b, fp_t c, fp_t d, fp_t alpha0, fp_t beta0, fp_t gamma0, fp_t delta0, unsigned maxIters, vector<fp_t>& coeff);

//Функция вычисления коэффициентов квадратных уравнений
template<typename fp_t>
unsigned int FQScoeffs(fp_t n, fp_t a, fp_t b, fp_t c, fp_t d, vector<fp_t>& coeffs);

// Итеративная функция нахождения корней уравнения методом бисекции
template<typename fp_t>
fp_t refine(fp_t left, fp_t right, fp_t* coef, fp_t(*fun)(fp_t, fp_t*));

/// Метод для нахождения индексов у массивов кубических корней от `betta0, gamma0`, таких что выполнилось условие 4.10
/// - Parameters:
///   - betta0CubeRoots: Кубические корни betta0
///   - gamma0CubeRoots: Кубические корни gamma0
///   - Q: Расчетный коэффициент Q
///   - beta0_I: Индекс корня от `beta0`
///   - gamma0_J: Индекс корня от `gamma0`
template<typename fp_t>
void check4_10Condition(vector<complex<fp_t>> betta0CubeRoots, vector<complex<fp_t>> gamma0CubeRoots, fp_t Q, int& beta0_I, int& gamma0_J);

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
void check4_11Condition(vector<complex<fp_t>> bettaSqrts, vector<complex<fp_t>> gammaSqrts, vector<complex<fp_t>> deltaSqrts, fp_t P, int& bettaI, int& gammaJ, int& deltaK);

// Метод проведения оценки точности результатов работы метода уравнения четвертой степени
template<typename fp_t>
void testQuarticPolynomial(
    int testCount, // Кол-во экспериментов
    long double maxDistance, // Максимальное расстояние между корнями
    int P1, // Кол-во кластеризированных корней
    int P2 // Кол-во кратных корней
);

#endif