#include <complex>
#include <iostream>
#include "/home/asudler/git/capstone/linalg/matrix/matrix.h"

using namespace std::complex_literals;

/* n.b.:
 * these checks are not super rigorous...
 * try to implement better checks in the future */

int main(void)
{

    matrix<double> A(5,4);
    for(int i = 0; i < A.size1; i++) {
        for(int j = 0; j < A.size2; j++) {
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
    E += 2.*E;
    E.print("E+=2E, so 3E:");

    matrix<std::complex<double>> F(5,4);
    for(int i = 0; i < F.size1; i++) {
        for(int j = 0; j < F.size2; j++) {
            F(i, j) = 1i + 2.*j + 1.;
        }
    }
    F.print("F:");

    matrix<std::complex<double>> G = F+A;
    G.print("G=F+A:");

    (F+A).print("F+A:");

    matrix<std::complex<double>> H(2,2);
    H(0,0) = 1.; H(0,1) = 2.; H(1,0) = 3.; H(1,1) = 4.;
    E += H;
    E.print("E+=H:");
    (E*H).print("E*H:");
    (E*H/3.).print("E*H/3.:");
    (0.5*H-2.*E*E).print("...", ',');

    std::vector<double> v = C.unravel();
    for(double x : v) {
        std::cout << x << ' ';
    }
    std::cout << '\n';

    stack(C.unravel(), C.size1).print("stack of vec above");

    std::cout << "row 2 of the mat above" << '\n';
    std::vector<double> v1 = C.row(2);
    for(double x : v1) {
        std::cout << x << ' ';
    }
    std::cout << '\n';
    std::cout << "column 2 of the mat above" << '\n';
    std::vector<double> v2 = C.column(2);
    for(double x : v2) {
        std::cout << x << ' ';
    }
    std::cout << '\n';

    matrix<std::complex<double>> rho(4,4);
    rho(1,1) = 0.5; rho(2,2) = 0.5;
    matrix<std::complex<double>> ham(4,4);
    ham(0,0) = 1.; ham(2,2) = -ham(0,0);
    ham(0,3) = 0.5 + 0.5i; ham(3,0) = conj(ham(0,3));
    ham(1,3) = 0.2 + 0.5i; ham(3,1) = conj(ham(1,3));
    ham(2,3) = 0.1 + 0.5i; ham(3,2) = conj(ham(2,3));
    rho.print("rho:"); ham.print("ham:");
    auto product = rho*ham; product.print("product");
    product = ham*rho; product.print("reverse order of product", 
            ' ', std::cerr);

    real(product).print("test real(product):");
    imag(product).print("test imag(product):");

    return 0;
}

