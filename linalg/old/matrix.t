#include <cassert>

template <typename T, typename U> // M3 = M1+M2
matrix<decltype(std::declval<T>()+std::declval<U>())>
operator+(const matrix<T>& M1, const matrix<U>& M2)
{
    assert(M1.size1 == M2.size1 && M1.size2 == M2.size2);
    matrix<decltype(std::declval<T>()+std::declval<U>())>
        M3(M1.size1, M1.size2);
    for(int i = 0; i < M1.size1*M1.size2; i++) M3.data[i] =
        M1.data[i] + M2.data[i];
    return M3;
}

template <typename T, typename U> // M3 = M1-M2
matrix<decltype(std::declval<T>()-std::declval<U>())>
operator-(const matrix<T>& M1, const matrix<U>& M2)
{
    assert(M1.size1 == M2.size1 && M1.size2 == M2.size2);
    matrix<decltype(std::declval<T>()-std::declval<U>())> 
        M3(M1.size1, M1.size2);
    for(int i = 0; i < M1.size1*M1.size2; i++) M3.data[i] =
        M1.data[i] - M2.data[i];
    return M3;
}

template <typename T, typename U> // M3 = M1*M2
matrix<decltype(std::declval<T>()*std::declval<U>())> 
operator*(const matrix<T>& M1, const matrix<U>& M2)
{
    assert(M1.size2 == M2.size1);
    matrix<decltype(std::declval<T>()*std::declval<U>())> 
        M3(M1.size1, M1.size2);
    for(int i = 0; i < M3.size1; i++)
    {
        for(int j = 0; j < M3.size2; j++)
        {
            T sum = 0;
            for(int k = 0; k < M1.size2; k++) {
                sum += M1(i,k)*M2(k,j);
            }
            M3(i,j) = sum;
        }
    }
    return M3;
}

template <typename T, typename U> // M2 = x*M1
matrix<decltype(std::declval<T>()*std::declval<U>())> 
operator*(const matrix<T>& M1, U x)
{
    matrix<decltype(std::declval<T>()*std::declval<U>())> 
        M2(M1.size1, M1.size2);
    for(int i = 0; i < M1.size1*M1.size2; i++) M2.data[i] = M1.data[i]*x;
    return M2;
}

template <typename T, typename U> // M2 = x*M1
matrix<decltype(std::declval<T>()*std::declval<U>())> 
operator*(T x, const matrix<U>& M1)
{
    matrix<decltype(std::declval<T>()*std::declval<U>())> M2 = M1*x;
    return M2;
}

template <typename T, typename U> // M2 = (1/x)*M1 or M1/x
matrix<decltype(std::declval<T>()/std::declval<U>())> 
operator/(const matrix<T>& M1, U x)
{
    matrix<decltype(std::declval<T>()/std::declval<U>())> 
        M2(M1.size1, M1.size2);
    for(int i = 0; i < M1.size1*M1.size2; i++) M2.data[i] = M1.data[i]/x;
    return M2;
}

template <typename T, typename U> // v2 = M*v1
std::vector<decltype(std::declval<T>()*std::declval<U>())> 
operator*(const matrix<T>& M, const std::vector<U>& v1)
{
    assert(M.size2 == v1.size());
    std::vector<decltype(std::declval<T>()*std::declval<U>())> v2(M.size1);
    for(int i = 0; i < M.size1; i++)
    {
        T sum = 0;
        for(int j = 0; j < M.size2; j++)
        {
            sum = M(i,j)*v1[j];
        }
        v2[i] = sum;
    }
    return v2;
}

template <typename T> // M is stacked row-by-row from v
matrix<T> stack(std::vector<T> v, int nrows)
{
    assert(v.size() % nrows == 0); // ensure vector can be stacked in n rows
    matrix<T> M(nrows, v.size()/nrows);
    for(int i = 0; i < M.size1; i++) 
    {
        for(int j = 0; j < M.size2; j++)
        {
            M.data[i + j*M.size1] = v[i*M.size2 + j];
        }
    }
    return M;
}

