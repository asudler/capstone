#include "matrix.h"

template <typename T> // default constructor
matrix<T>::matrix() : size1(0), size2(0), data(nullptr) {}
/*{
    size1 = 0;
    size2 = 0;
    data = nullptr;
}*/ 

template <typename T> // parameterized constructor
matrix<T>::matrix(size_t n, size_t m) : size1(n), size2(m), data(n*m) {}
/*{
    size1 = n; size2 = m;
    data(size1*size2); 
}*/ // parameterized constructor

template <typename T> 
matrix<T>::matrix(const matrix& that) : data(that.data) {}
/*{
    size1 = that.size1;
    size2 = that.size2;
    std::vector<T> data_(size1*size2);
    for(int i = 0; i < size1*size2; i++) data[i] = that.data[i];
}*/ // copy constructor

template <typename T>
matrix<T>

