#include <complex>
#include "matrix.h"

using namespace std::complex_literals;

/* n.b.:
 * these checks are not super rigorous...
 * try to implement better checks in the future */

int main(void)
{
    matrix<double> A(5,4);
    for (int i = 0; i < A.size1; i++) {
        for (int j = 0; j < A.size2; j++) {
            A(i, j) = i + 2.*j + 1;
        }
    }
    A.print("A:");

    matrix<double> B = A.transpose();
    B.print("B: (should be A.transpose)");

    matrix<double> C(B);
    C.print("C: (should be B)");

    matrix<double> D = B;
    D.print("D: (should be B)");

    matrix<std::complex<double>> E(2,2);
    E(0,0) = 1; E(0,1) = 1i; E(1,0) = -1i; E(1,1)= -1;
    E.print("E:");
    E.transpose().print("E.transpose():");

    return 0;
}
