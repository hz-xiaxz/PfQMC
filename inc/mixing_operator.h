#ifndef MIXING_OPERATOR_H
#define MIXING_OPERATOR_H

#include "operator.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <Eigen/Dense>

// Operator that handles swapping between replicas at the boundary
class MixingOperator : public Operator {
public:
    int nReplicas;
    int nDimSingle;
    int nDim;
    std::vector<bool> is_swapped; // true if pair p is swapped
    rdGenerator* rd;

    // Constructor
    MixingOperator(int _nReplicas, int _nDimSingle, rdGenerator* _rd) 
        : nReplicas(_nReplicas), nDimSingle(_nDimSingle), rd(_rd) {
        nDim = nReplicas * nDimSingle;
        is_swapped.resize(nReplicas / 2, false); // Initialize as Identity
    }

    // Destructor
    ~MixingOperator() override {}

    // --- Operator Interface Implementation ---

    // Apply the operator matrix B to A: B = Op * A
    void left_multiply(const MatType &A, MatType &B) override {
        // We need to copy A to B first, but with the swaps applied.
        // Or just write directly to B.
        // B must be same size as A.
        
        for(int p=0; p < nReplicas/2; ++p) {
            int rA = 2*p;
            int rB = 2*p+1;
            int offsetA = rA * nDimSingle;
            int offsetB = rB * nDimSingle;
            
            if(is_swapped[p]) {
                // Row rA in B comes from Row rB in A
                B.block(offsetA, 0, nDimSingle, A.cols()) = A.block(offsetB, 0, nDimSingle, A.cols());
                // Row rB in B comes from Row rA in A
                B.block(offsetB, 0, nDimSingle, A.cols()) = A.block(offsetA, 0, nDimSingle, A.cols());
            } else {
                B.block(offsetA, 0, nDimSingle, A.cols()) = A.block(offsetA, 0, nDimSingle, A.cols());
                B.block(offsetB, 0, nDimSingle, A.cols()) = A.block(offsetB, 0, nDimSingle, A.cols());
            }
        }
    }

    // B = A * Op
    void right_multiply(const MatType &A, MatType &B) override {
        for(int p=0; p < nReplicas/2; ++p) {
            int rA = 2*p;
            int rB = 2*p+1;
            int offsetA = rA * nDimSingle;
            int offsetB = rB * nDimSingle;
            
            if(is_swapped[p]) {
                // Col rA in B comes from Col rB in A
                B.block(0, offsetA, A.rows(), nDimSingle) = A.block(0, offsetB, A.rows(), nDimSingle);
                // Col rB in B comes from Col rA in A
                B.block(0, offsetB, A.rows(), nDimSingle) = A.block(0, offsetA, A.rows(), nDimSingle);
            } else {
                B.block(0, offsetA, A.rows(), nDimSingle) = A.block(0, offsetA, A.rows(), nDimSingle);
                B.block(0, offsetB, A.rows(), nDimSingle) = A.block(0, offsetB, A.rows(), nDimSingle);
            }
        }
    }

    // Since S = S^{-1} = S^T, inverse and adjoint are the same as the operator itself.
    void inv_left_multiply(const MatType &A, MatType &B) override { left_multiply(A, B); }
    void adjoint_left_multiply(const MatType &A, MatType &B) override { left_multiply(A, B); }
    void inv_right_multiply(const MatType &A, MatType &B) override { right_multiply(A, B); }
    void adjoint_inv_right_multiply(const MatType &A, MatType &B) override { right_multiply(A, B); }

    // Propagate: g = Op * g * Op^{-1}
    // Since Op = Op^{-1}, g = S * g * S.
    // This swaps both rows and columns.
    void left_propagate(MatType &g, MatType &gTmp) override {
        // gTmp = S * g
        left_multiply(g, gTmp);
        // g = gTmp * S
        right_multiply(gTmp, g);
    }

    void right_propagate(MatType &g, MatType &gTmp) override {
        // Same as left_propagate because S = S^{-1}
        left_propagate(g, gTmp);
    }

    // Stabilized multiply for UDT
    void stabilizedLeftMultiply(UDT &F) override {
        if (F.isBlock) {
            // Handle block UDT
            // We assume the block structure aligns with the mixing operator's structure.
            // Specifically, for 4 replicas and nBlocks=2, each block contains a pair of replicas (0,1) and (2,3).
            // The mixing operator swaps replicas within each pair.
            
            for (size_t b = 0; b < F.blocks.size(); b++) {
                // Assuming 1 pair per block for now (consistent with nBlocks=2 for 4 replicas)
                int p = b; 
                if (p < is_swapped.size() && is_swapped[p]) {
                    MatType& U = F.blocks[b].U;
                    int halfSize = U.rows() / 2;
                    // Swap top and bottom halves of U
                    for (int i = 0; i < halfSize; i++) {
                        U.row(i).swap(U.row(i + halfSize));
                    }
                }
            }
        } else {
            // F = S * F
            // We can just swap rows in F.U
            MatType U_new = F.U; // Copy
            left_multiply(F.U, U_new);
            F.U = U_new;
        }
    }

    // --- Update Logic ---

    // Returns the sign change of the weight
    DataType update(MatType &g) override {
        DataType signRatio = 1.0;
        
        for (int p = 0; p < nReplicas / 2; ++p) {
            // Propose flip
            // The update kernel K for switching between I and S (or S and I) is:
            // K = S_{sub} (I_{sub} - G_{sub}) + G_{sub}
            // where sub refers to the 2N x 2N block for the pair.
            
            int rA = 2 * p;
            int rB = 2 * p + 1;
            int offsetA = rA * nDimSingle;
            int offsetB = rB * nDimSingle;
            int N = nDimSingle;

            // Extract G blocks
            MatType G_AA = g.block(offsetA, offsetA, N, N);
            MatType G_AB = g.block(offsetA, offsetB, N, N);
            MatType G_BA = g.block(offsetB, offsetA, N, N);
            MatType G_BB = g.block(offsetB, offsetB, N, N);
            
            MatType I = MatType::Identity(N, N);
            
            // Construct K matrix (2N x 2N)
            // K = [ G_AA - G_BA,       I - G_BB + G_AB ]
            //     [ I - G_AA + G_BA,   G_BB - G_AB     ]
            MatType K(2*N, 2*N);
            K.topLeftCorner(N, N) = G_AA - G_BA;
            K.topRightCorner(N, N) = I - G_BB + G_AB;
            K.bottomLeftCorner(N, N) = I - G_AA + G_BA;
            K.bottomRightCorner(N, N) = G_BB - G_AB;
            
            // Calculate determinant (ratio)
            DataType r = K.determinant();
            
            // Metropolis check
            double prob = std::abs(r);
            double rand_val = rd->rdUniform01();
            
            if (rand_val < prob) {
                // Accept
                is_swapped[p] = !is_swapped[p];
                signRatio *= (r / prob); // Phase change
                
                // Update Green's function
                // G' = G - L K^{-1} R
                
                MatType K_inv = K.inverse();
                
                // Construct L (N_total x 2N)
                // L = (I - G) restricted to columns rA, rB.
                // L = -G_cols except diagonal elements +1.
                MatType L(nDim, 2*N);
                L.leftCols(N) = -g.block(0, offsetA, nDim, N);
                L.rightCols(N) = -g.block(0, offsetB, nDim, N);
                
                // Add Identity to diagonal parts of L
                for(int i=0; i<N; ++i) {
                    L(offsetA + i, i) += 1.0;
                    L(offsetB + i, N + i) += 1.0;
                }
                
                // Construct R (2N x N_total)
                // R = [ -G_A + G_B ]
                //     [  G_A - G_B ]
                MatType R(2*N, nDim);
                R.topRows(N) = -g.block(offsetA, 0, N, nDim) + g.block(offsetB, 0, N, nDim);
                R.bottomRows(N) = g.block(offsetA, 0, N, nDim) - g.block(offsetB, 0, N, nDim);
                
                // Update G
                g -= L * K_inv * R;
            }
        }
        
        return signRatio;
    }

    void getGreensMat(MatType &g0) override {
        // Not used in this context
    }
    
    DataType getSignOfWeight() override {
        return 1.0;
    }
};

// Trivial Swap Operator: Always swaps pairs (0,1), (2,3), etc.
class TrivialSwapOperator : public MixingOperator {
public:
    TrivialSwapOperator(int _nReplicas, int _nDimSingle, rdGenerator* _rd)
        : MixingOperator(_nReplicas, _nDimSingle, _rd) {
        // Set all to false (Identity) to verify no mixing
        std::fill(is_swapped.begin(), is_swapped.end(), false);
    }

    DataType update(MatType &g) override {
        // Do nothing, keep swapped
        return 1.0;
    }
};

#endif // MIXING_OPERATOR_H
