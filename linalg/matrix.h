/* n.b.:
 * This library is modeled after Prof. Dmitri Fedorov's matrix class,
 * see http://212.27.24.106:8080/prog/c++/matrix/.
 * Additionally, several of the implementations are adapted from classmate
 * Luca Janken's public git repository, 
 * see https://github.com/LucaJanken/PPNM/tree/main/homework/matrix/.
 * (hopefully the links are working...) */

#ifndef MATRIX_H
#define MATRIX_H

#include <stdexcept>
#include <string>
#include <vector>

/* ideas for the future:
 * (1) maybe organize into public, protected, private members
 * (2) implement type casting  */

template <typename T>
struct matrix
{
    size_t size1, size2; std::vector<T> data;

    /* constructors, destructors, etc. */
    matrix(); // default constructor
    matrix(size_t, size_t); // parameterized constructor
    matrix(const matrix&); // copy constructor
    matrix(const std::vector<std::vector<T>>&); // copy vector arr into matrix
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
    void print(std::string s="", char delimiter=' ') const; // print matrix
    
    /* struct member operators */
    matrix& operator+=(const matrix&);
    matrix& operator-=(const matrix&);
    matrix& operator*=(const matrix&);
    matrix& operator*=(T);
    matrix& operator/=(T);
    // matrix& operator^(int); // tbd lol
}; // matrix

/* non-member operators (not sure of the technical term) */ 
template <typename T, typename U>
matrix<decltype(std::declval<T>()+std::declval<U>())>
operator+(const matrix<T>&, const matrix<U>&); // addition
template <typename T, typename U> 
matrix<decltype(std::declval<T>()-std::declval<U>())> 
operator-(const matrix<T>&, const matrix<U>&); // subtraction
template <typename T, typename U> 
matrix<decltype(std::declval<T>()*std::declval<U>())> 
operator*(const matrix<T>&, const matrix<U>&); // matrix*matrix
template <typename T, typename U> 
matrix<decltype(std::declval<T>()*std::declval<U>())> 
operator*(const matrix<T>&, U x); // matrix*number
template <typename T, typename U> 
matrix<decltype(std::declval<T>()*std::declval<U>())> 
operator*(T x, const matrix<U>&); // number*matrix
template <typename T, typename U> 
matrix<decltype(std::declval<T>()/std::declval<U>())> 
operator/(const matrix<T>&, U x); // matrix/number
template <typename T, typename U> 
std::vector<decltype(std::declval<T>()*std::declval<U>())> 
operator*(const matrix<T>&, const std::vector<U>&); // matrix*vector

/* non-member functions */
// template <typename T> // read matrix from input filestream
// matrix<T> read(string::fname, char delimiter=',', int skip_header=0);

#include "matrix.t" // must define public templates in header file

#endif // MATRIX_H. Thanks Dmitri!

