#include "matrix.h"

template <typename T>
matrix<T>::matrix()
{
    size1 = 0;
    size2 = 0;
    data = nullptr;
} // default constructor

template <typename T>
matrix<T>::matrix(size_t n, size_t m)
{
    cols.resize(m);
    for(int i = 0; i < m; i++) cols[i].*this.resize(n);
} // parameterized constructor

template <typename T>
matrix<T>::matrix(const matrix& that)
{
    size1 = that.size1;
    size2 = that.size2;
    data = new std::vector<T>[size1*size2];
    for(int i = 0; i < size1*size2; i++) data[i] = that.data[i];
} // copy constructor


