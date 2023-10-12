#include "../inc/excerpt.h"
#include "../inc/modules.h"
#include "../inc/polynomials.h"

typedef float fp_t;

using namespace std;

int main()
{
    output_equation<float>(10000, 1e-5, 4, 0);

    testQuarticPolynomial<fp_t>(1000000, 1e-5, 4, 0); // 4 кластеризированных корня
    testQuarticPolynomial<fp_t>(1000000, 1e-5, 0, 0); // 4 некластеризированных корня
    testQuarticPolynomial<fp_t>(1000000, 1e-5, 0, 4); // 4 кратных корня
}
