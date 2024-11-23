#include <cassert>
#include <complex>
#include <iostream>
#include "matrix.h"

/* explicit instantiation of matrix templates */
template struct matrix<double>;
template struct matrix<std::complex<double>>;

template <typename T> // default constructor
matrix<T>::matrix() : size1(0), size2(0), data(0) {}

template <typename T> // parameterized constructor
matrix<T>::matrix(size_t n, size_t m) : size1(n), size2(m), data(n*m) {}

template <typename T> // copy constructor
matrix<T>::matrix(const matrix& other) : size1(other.size1), size2(other.size2) 
{
    data.resize(other.size1*other.size2); // allocate memory for new object
    std::copy(other.data.begin(), other.data.end(), data.begin()); // copy
}

template <typename T> // copy vector of vectors into matrix
matrix<T>::matrix(const std::vector<std::vector<T>>& other) : 
    size1(other.size()), size2(other[0].size())
{
    std::vector<T> help = other[0]; // 1D array needed
    for(int i = 1; i < other.size(); i++) // concatenate vector of vectors
    {
        if(other[i].size()!=other[0].size()) throw std::range_error("matrix: "
            "no jagged matrices allowed");
        help.insert(help.end(), other[i].begin(), other[i].end());
    } 
    data.resize(help.size()); // allocate mem
    std::copy(help.begin(), help.end(), data.begin()); // copy
}

template <typename T> // copy assignment
matrix<T>& matrix<T>::operator=(const matrix& other)
{
    if(this == &other) return *this; // self-assignment check

    // resize and copy
    data.resize(other.size1*other.size2);
    std::copy(other.data.begin(), other.data.end(), data.begin());

    // update sizes
    size1 = other.size1;
    size2 = other.size2;

    return *this;
}

template <typename T> // move constructor
matrix<T>::matrix(matrix&& other) noexcept : 
    size1(other.size1), size2(other.size2)
{
    data = std::move(other.data); // move data from "other" to "this"
    other.size1 = 0; // reset "other" to empty state
    other.size2 = 0; 
    other.data.clear();
}

template <typename T> // move assignment
matrix<T>& matrix<T>::operator=(matrix&& other) noexcept
{
    if(this == &other) return *this; // self-assignment check

    // transfer ownership of data
    data = std::move(other.data);
    size1 = other.size1;
    size2 = other.size2;

    // leave moved-from object in unspecified state
    other.size1 = 0;
    other.size2 = 0;
    other.data.clear();

    return *this;
}

template <typename T> // destructor
matrix<T>::~matrix() {} 

template <typename T> // setter
void matrix<T>::set(size_t i, size_t j, T x) { data[i + j*size1] = x; }

template <typename T> // getter
T matrix<T>::get(size_t i, size_t j) const { return data[i + j*size1]; } 

template <typename T> // operator ()
T& matrix<T>::operator()(size_t i, size_t j) { return data[i + j*size1]; }

template <typename T> // operator () again but read only (is it necessary?)
const T& matrix<T>::operator()(size_t i, size_t j) const 
{ 
    return data[i + j*size1]; 
}

template <typename T> // copy member function
matrix<T> matrix<T>::copy() { matrix<T> M(*this); return M; }

template <typename T> // transpose member function
matrix<T> matrix<T>::transpose()
{
    matrix<T> M(size2, size1);
    for(int i = 0; i < size1; i++)
    {
        for(int j = 0; j < size2; j++)
        {
            M.data[j + i*size2] = data[i + j*size1];
        }
    }
    return M;
}

template <typename T> // unravel matrix row-by-row
std::vector<T> matrix<T>::unravel()
{
    std::vector<T> v(size1*size2);
    for(int i = 0; i < size1; i++)
    {
        for(int j = 0; j < size2; j++)
        {
            v[i*size2 + j] = data[i + j*size1];
        }
    }
    return v;
}

template <typename T> // returns row m of the matrix
std::vector<T> matrix<T>::row(int m)
{
    std::vector<T> v(size2);
    for(int i = 0; i < size2; i++) v[i] = data[m + i*size1];
    return v;
}

template <typename T> // returns column n of the matrix
std::vector<T> matrix<T>::column(int n)
{
    std::vector<T> v(size1);
    for(int i = 0; i < size1; i++) v[i] = data[i + n*size1];
    return v;
}

template <typename T> // print member function
void matrix<T>::print(
    std::string s, 
    char delimiter, 
    std::ostream& output_stream
) const 
{
    output_stream << s.c_str() << '\n';
    for(int i = 0; i < size1; i++)
    {
        for(int j = 0; j < size2; j++) 
        {
            output_stream << (*this)(i, j) << delimiter;
        }
        output_stream << '\n';
    }
}

template <typename T> // +=
matrix<T>& matrix<T>::operator+=(const matrix& other)
{
    assert(size1 == other.size1 && size2 == other.size2);
    for(int i = 0; i < size1*size2; i++) data[i] += other.data[i];
    return *this;
}

template <typename T> // -=
matrix<T>& matrix<T>::operator-=(const matrix& other)
{
    assert(size1 == other.size1 && size2 == other.size2);
    for(int i = 0; i < size1*size2; i++) data[i] -= other.data[i];
    return *this;
}

template <typename T> // *= (matrix multiplication on 1D "matrices")
matrix<T>& matrix<T>::operator*=(const matrix& other)
{
    assert(size2 = other.size1);
    matrix M = matrix(size1, other.size2);
    for(int i = 0; i < M.size1; i++)
    {
        for(int j = 0; j < M.size2; j++)
        {
            T sum = 0;
            for(int k = 0; k < size2; k++) {
                sum += data[i*size2 + k]*other.data[k*other.size2 + j];
            }
            M.data[i*M.size2 + j] = sum;
        }
    }
    (*this) = M;
    return *this;
}

template <typename T> // *= (x*matrix)
matrix<T>& matrix<T>::operator*=(T x)
{
    for(int i = 0; i < size1*size2; i++) data[i] *= x;
    return *this;
}

template <typename T> // *= (x*matrix)
matrix<T>& matrix<T>::operator/=(T x)
{
    for(int i = 0; i < size1*size2; i++) data[i] /= x;
    return *this;
}

matrix<double> real(matrix<std::complex<double>> N)
{
    matrix<double> M(N.size1, N.size2);
    for(int i = 0; i < N.size1; i++)
    {
        for(int j = 0; j < N.size2; j++)
        {
            M(i,j) = real(N(i,j));
        }
    }
    return M;
}

matrix<double> imag(matrix<std::complex<double>> N)
{
    matrix<double> M(N.size1, N.size2);
    for(int i = 0; i < N.size1; i++)
    {
        for(int j = 0; j < N.size2; j++)
        {
            M(i,j) = imag(N(i,j));
        }
    }
    return M;
}


