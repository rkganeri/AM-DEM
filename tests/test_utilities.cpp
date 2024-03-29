#include <string>   
#include <limits>
#include <float.h>
#include <cmath>

#include "gtest/gtest.h"
#include "Kokkos_Core.hpp"
#include "fmt/core.h"

#include "utilities.hpp"

// wrapping in the namespace so I can omit prepending to all classes/functions
namespace amdem {

namespace utilities {

// tests the norm2 inline function
TEST(utilities,norm2) {

    // first test the c-array version of the fn
    int size = 4;
    double a[size];
    for (int i=0; i<size; i++) {
        a[i] = (i+1)*2;
    }

    // use default length of 3
    double val = norm2(a);
    double expected_val = sqrt(56); // 2^2 + 4^2 + 6^2
    EXPECT_FLOAT_EQ(val, expected_val);

    // use full length of 4
    val = norm2(a, size);
    expected_val = sqrt(120);
    EXPECT_FLOAT_EQ(val, expected_val);

    // now try passing in components individually (2nd definition of norm2)
    double a0 = 3;
    double a1 = -1.2;
    double a2 = 2.3;
    val = norm2(a0, a1, a2);
    expected_val = sqrt(15.73);
    EXPECT_FLOAT_EQ(val, expected_val);
    
}

// tests the dotProduct inline function
TEST(utilities,dotProduct) {

    int size = 4;
    double a[size];
    double b[size];
    for (int i=0; i<size; i++) {
        a[i] = (i+1)*2;
        b[i] = -(i+1)*1.5;
    }

    // use default length of 3
    double val = dotProduct(a,b);
    double expected_val = -42;
    EXPECT_FLOAT_EQ(val, expected_val);

    val = dotProduct(b,a);
    EXPECT_FLOAT_EQ(val, expected_val);

    // use full length of 4
    val = dotProduct(a, b, size);
    expected_val = -90;
    EXPECT_FLOAT_EQ(val, expected_val);

    // try length of 2 cause why not
    val = dotProduct(b, a, 2);
    expected_val = -15;
    EXPECT_FLOAT_EQ(val, expected_val);
    
}

// tests the min inline functions
TEST(utilities,min_max) {

    int a = 3;
    int b = 2;

    double c = -1.217;
    double d = -0.76;

    EXPECT_EQ(min(a,b), b);
    EXPECT_FLOAT_EQ(min(c,d), c);

    EXPECT_EQ(min(7,6), 6);
    EXPECT_FLOAT_EQ(min(4.708, 9.1), 4.708);

    EXPECT_EQ(max(a,b), a);
    EXPECT_FLOAT_EQ(max(c,d), d);

    EXPECT_EQ(max(-1,0), 0);
    EXPECT_FLOAT_EQ(max(3.76,2.12), 3.76);
    
}

} // namespace utilities

} // namespace amdem




