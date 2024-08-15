/* n.b.:
 * This library is modeled after Prof. Dmitri Fedorov's matrix class,
 * see http://212.27.24.106:8080/prog/c++/matrix/.
 * Additionally, several of the implementations are adapted from classmate
 * Luca Janken's public git repository, 
 * see https://github.com/LucaJanken/PPNM/tree/main/homework/matrix/.
 * (hopefully the links are working...) */

#ifndef MATRIX_H
#define MATRIX_H

#include <string>
#include <vector>

/* n.b.:
 * maybe organize into public, protected, private members in the future */

template <typename T>
struct matrix
{
    size_t size1, size2; std::vector<T> data;

    /* constructors, destructors, etc. */
    matrix(); // default constructor
    matrix(size_t, size_t); // parameterized constructor
    matrix(const matrix&); // copy constructor
    matrix& operator=(const matrix&); // copy assignment
    matrix(matrix&&) noexcept; // move constructor
    matrix& operator=(matrix&&) noexcept; // move assignment
    ~matrix(); // destructor

    /* getters, setters */
    void set(size_t i, size_t j, T x); // setter
    T get(size_t i, size_t j) const; // getter
    T& operator()(size_t i, size_t j); // access matrix element, e.g. M(i,j)
    const T& operator()(size_t i, size_t j) const; // ^ but read-only
    // T& operator[](size_t i, size_t j); // not possible unless C++23+ :(

    /* member functions */
    matrix copy(); // copy
    matrix transpose(); // transpose (swap rows and columns)
    void print(std::string s="") const; // print matrix, e.g. M.print("M:")

    /* struct member operators */
    matrix& operator+=(const matrix&);
    matrix& operator-=(const matrix&);
    matrix& operator*=(const matrix&);
    matrix& operator*=(T);
    matrix& operator/=(T);
    // matrix& operator^(int); // tbd lol

    matrix operator+(const matrix&, const matrix&);
}; // matrix

/* non-member operators */ 
//template <typename T> matrix<T> operator+(const matrix&, const matrix&);
/*
matrix operator-(const matrix&, const matrix&);
matrix operator*(matrix&, matrix&);
matrix operator*(const matrix&, double x);
matrix operator*(double x, const matrix&);
matrix operator/(const matrix&, double x);
*/
#endif // MATRIX_H. Thanks Dmitri!

