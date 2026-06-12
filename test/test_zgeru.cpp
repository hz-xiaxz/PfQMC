#include "../inc/types.h"
#include <iostream>

int main() {
    int n = 2;
    DataType alpha = 1.0;
    DataType x[2] = {1.0, 2.0};
    DataType y[2] = {3.0, 4.0};
    DataType a[4] = {0.0, 0.0, 0.0, 0.0};
    int inc = 1;

    zgeru(&n, &n, &alpha, x, &inc, y, &inc, a, &n);

    std::cout << "zgeru called successfully" << std::endl;
    return 0;
}
