#include <gtest/gtest.h>
#include "../inc/honeycomb.h"
#include "../inc/pfqmc.h"

TEST(SingleMajoranaTest, GreenFunctionUpdate) {
    mkl_set_num_threads(1);

    int Lx = 4;
    int Ly = 4;
    int LTau = 10;
    int hamiltonianDim = Lx * Ly * 4;
    const MatType identity = MatType::Identity(hamiltonianDim, hamiltonianDim);
    // Enable singleMaj = true
    SpinlessTvHoneycombUtils config(Lx, Ly, 0.1, 0.7, LTau, true);
    rdGenerator rd(114514);
    Honeycomb_tV walker(&config, &rd);
    MatType g;
    MatType A = identity;

    int l = 4 * LTau - 1;
    for (int i = 0; i < l; i++) {
        walker.op_array[i]->left_multiply(A, A);
    }
    walker.op_array[4 * LTau]->right_multiply(A, A);
    walker.op_array[l]->right_multiply(A, A);

    g = (2 * ((identity + A).inverse()));
    
    // Perform updates using singleMajoranaFlip
    DataType sign = 1.0;
    SpinlessVOperator* p = (SpinlessVOperator*) walker.op_array[l];
    
    // Ensure singleMaj is true
    EXPECT_TRUE(p->singleMaj);

    // Try to flip some aux fields
    // We can't easily force an accept, but we can run update and check consistency
    // update() calls singleFlipSingleMajorana because singleMaj is true
    p->update(g);

    // Recompute G from scratch
    // walker.op_array[l] holds a pointer to s, so s is already updated
    
    A = identity;
    for (int i = 0; i < l; i++) {
        walker.op_array[i]->left_multiply(A, A);
    }
    walker.op_array[4 * LTau]->right_multiply(A, A);
    walker.op_array[l]->right_multiply(A, A);
    
    MatType gBrutal = 2 * ((identity + A).inverse());

    double r = (gBrutal - g).squaredNorm();
    EXPECT_NEAR(r, 0.0, 1e-10);
}
