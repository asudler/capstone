/* n.b.:
 * This library is modeled after Prof. Dmitri Fedorov's matrix class
 * see http://212.27.24.106:8080/prog/c++/matrix/
 * Additionally, several of the implementations are adapted from classmate
 * Luca Janken's public git repository 
 * see https://github.com/LucaJanken/PPNM/tree/main/homework/matrix 
 * (hopefully links are working...) */

#ifndef MATRIX_H
#define MATRIX_H

#include <vector>

template <typename T>
struct matrix
{
    std::vector<std::vector<T>> cols;
    size_t size1() const { return cols[0].size(); } // # of rows
    size_t size2() const { return cols.size(); } // # of columns

    matrix(); // default constructor
    matrix(size_t, size_t); // parameterized constructor
    matrix(const matrix&); // copy constructor
    matrix(matrix&&); // move constructor
    ~matrix(); // destructor
/*
    matrix& operator=(const matrix&); // copy assignment
    matrix& operator=(matrix &&); // move assignment
    matrix copy(); 
    matrix transpose(); 

    void set (size_t i, size_t j, T x) { cols[j][i] = x; }
    T get (size_t i, size_t j) const { return cols[j].get(i); }
    T& operator() (size_t i, size_t j) { return cols[j][i]; }
    T& operator[] (size_t i, size_t j) { return cols[j][i]; }
    std::vector<T> get_col(int j);
    void set_col(int j, std::vector<T>& cj);

    matrix& operator+=(const matrix&);
    matrix& operator-=(const matrix&);
    matrix& operator*=(const matrix&);
    matrix& operator*=(double);
    matrix& operator/=(double);
    matrix& operator^(int);

    void print(const char* s="") const; */
}; // matrix
/*
matrix operator+(const matrix&, const matrix&);
matrix operator-(const matrix&, const matrix&);
matrix operator*(const matrix&, const matrix&);
matrix operator*(const matrix&, const double x);
matrix operator*(const double x, const matrix&);
matrix operator*(const matrix&, const std::complex<double> x);
matrix operator*(const std::complex<double> x, const matrix&);
matrix operator/(const matrix&, const double x);
matrix operator/(const matrix&, const std::complex<double> x);

std::vector<double> operator*(matrix&, std::vector<double>&); // !!!!! ????
*/
#endif // MATRIX_H. Thanks Dmitri!
