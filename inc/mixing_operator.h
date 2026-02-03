#ifndef MIXING_OPERATOR_H
#define MIXING_OPERATOR_H

#include "operator.h"
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <vector>

// Operator that handles swapping between replicas at the boundary
// The SWAP operator for fermions: S^{αβ}[σ] = exp(π/4 σ ψ^α ψ^β)
// satisfies: S^† (ψ^α, ψ^β)^T S = σ (ψ^β, -ψ^α)^T
// In matrix form for 2N×2N block: S = [0, I; -I, 0]
class MixingOperator : public Operator {
public:
  int nReplicas;
  int nDimSingle;
  int nDim;
  std::vector<std::vector<int>>
      sigma; // sigma[p][m] is the aux field for pair p, mode m
  rdGenerator *rd;

  // Constructor
  MixingOperator(int _nReplicas, int _nDimSingle, rdGenerator *_rd)
      : nReplicas(_nReplicas), nDimSingle(_nDimSingle), rd(_rd) {
    nDim = nReplicas * nDimSingle;
    sigma.resize(nReplicas / 2);
    for (int p = 0; p < nReplicas / 2; ++p) {
      sigma[p].resize(nDimSingle);
      for (int m = 0; m < nDimSingle; ++m) {
        // Initialize randomly +1 or -1
        // Ensure all pairs have the same sigma
        if (p == 0) {
          sigma[p][m] = (rd->rdUniform01() < 0.5) ? 1 : -1;
        } else {
          sigma[p][m] = sigma[0][m];
        }
      }
    }
  }

  // Destructor
  ~MixingOperator() override {}

  // --- Operator Interface Implementation ---
  // Note: A and B are nDim × nDim matrices (block-diagonal with 2N×2N blocks
  // per pair) The SWAP only acts on each 2N×2N diagonal block independently

  // Apply the operator matrix: B = Op * A
  // NEW DEFINITION: B = i * sigma * tau_y
  // tau_y = [0, -i; i, 0]
  // B = sigma * [0, 1; -1, 0]
  void left_multiply(const MatType &A, MatType &B) override {
    for (int p = 0; p < nReplicas / 2; ++p) {
      int offsetA = 2 * p * nDimSingle;       // row offset for replica α
      int offsetB = (2 * p + 1) * nDimSingle; // row offset for replica β
      int blockCol =
          2 * p * nDimSingle; // column offset for this pair's 2N block

      int parity = (p % 2 == 0) ? 1 : -1;
      for (int m = 0; m < nDimSingle; ++m) {
        int s = sigma[p][m] * parity;
        // B * A = s * [0, 1; -1, 0] * [row_α; row_β]
        // row_α_new = s * row_β
        // row_β_new = -s * row_α

        B.row(offsetA + m).segment(blockCol, 2 * nDimSingle) =
            (double)s * A.row(offsetB + m).segment(blockCol, 2 * nDimSingle);

        B.row(offsetB + m).segment(blockCol, 2 * nDimSingle) =
            -(double)s * A.row(offsetA + m).segment(blockCol, 2 * nDimSingle);
      }
    }
  }

  // B = A * S
  void right_multiply(const MatType &A, MatType &B) override {
    for (int p = 0; p < nReplicas / 2; ++p) {
      int offsetA = 2 * p * nDimSingle;
      int offsetB = (2 * p + 1) * nDimSingle;
      int blockRow = 2 * p * nDimSingle; // row offset for this pair's 2N block

      int parity = (p % 2 == 0) ? 1 : -1;
      for (int m = 0; m < nDimSingle; ++m) {
        int s = sigma[p][m] * parity;
        // A * B = [col_α, col_β] * s * [0, 1; -1, 0]
        // col_α_new = -s * col_β
        // col_β_new = s * col_α

        B.col(offsetA + m).segment(blockRow, 2 * nDimSingle) =
            -(double)s * A.col(offsetB + m).segment(blockRow, 2 * nDimSingle);

        B.col(offsetB + m).segment(blockRow, 2 * nDimSingle) =
            (double)s * A.col(offsetA + m).segment(blockRow, 2 * nDimSingle);
      }
    }
  }

  // Stabilized multiply for UDT
  void stabilizedLeftMultiply(UDT &F) override {
    // We cannot easily use the block structure of F because S mixes rows.
    // So we apply to the full matrix U.
    MatType U_new(F.U.rows(), F.U.cols());
    left_multiply(F.U, U_new);
    F.U = U_new;
  }

  // Calculate ratio for flipping sigma[p][m]
  DataType getRatio(MatType &g, int p, int m) {
    int offsetA = 2 * p * nDimSingle;
    int offsetB = (2 * p + 1) * nDimSingle;
    int parity = (p % 2 == 0) ? 1 : -1;
    int s = sigma[p][m] * parity; // Current state with parity

    // Extract G_12 (element connecting mode m of replica 1 and 2)
    int idxA = offsetA + m;
    int idxB = offsetB + m;

    // User formula: ratio = 1 + sigma_m * G_12
    // Assuming G_12 is g(idxA, idxB)
    DataType G12 = g(idxA, idxB);

    // Note: The formula 1 + s * G12 implies a specific update relation.
    // We use it as requested.
    return 1.0 + (double)s * G12;
  }

  // Perform local update restricted to the replica pair block
  // This is now the default and only update method
  void localUpdate(MatType &g, int p, int m) {
    int blkOffset = 2 * p * nDimSingle;
    int blkSize = 2 * nDimSingle;
    DataType *g_block_ptr = g.data() + blkOffset + blkOffset * nDim;
    int inc = 1;

    // We only work on the block starting at blkOffset
    // Indices relative to the block
    int idxA_local = m;
    int idxB_local = nDimSingle + m;

    int parity = (p % 2 == 0) ? 1 : -1;
    int s = sigma[p][m] * parity;

    // Extract G12 from the block (same as global since g is block diagonal)
    DataType G12 = g(blkOffset + idxA_local, blkOffset + idxB_local);
    DataType X = 2.0 * (double)s + G12;

    DataType D = 1.0 + X * X;
    DataType invD = 1.0 / D;

    DataType M00 = invD;
    DataType M01 = -X * invD;
    DataType M10 = X * invD;
    DataType M11 = invD;

    // c = 2e - g.col
    cVecType c0 = -g.col(blkOffset + idxA_local).segment(blkOffset, blkSize);
    c0(idxA_local) += 2.0;

    cVecType c1 = -g.col(blkOffset + idxB_local).segment(blkOffset, blkSize);
    c1(idxB_local) += 2.0;

    cVecType u0 = c0 * M00 + c1 * M10;
    cVecType u1 = c0 * M01 + c1 * M11;

    // Update the block: G_block += U * V^T
    // We use zgeru for rank-1 updates: A := alpha * x * y^T + A
    // G_block += u0 * c0^T + u1 * c1^T
    DataType alpha = 1.0;
    zgeru(&blkSize, &blkSize, &alpha, u0.data(), &inc, c0.data(), &inc,
          g_block_ptr, &nDim);
    zgeru(&blkSize, &blkSize, &alpha, u1.data(), &inc, c1.data(), &inc,
          g_block_ptr, &nDim);

    // Flip the spin
    sigma[p][m] = -s;
  }

  // --- Update Logic ---
  // Returns the sign change of the weight
  DataType update(MatType &g) override {
    DataType signRatio = 1.0;

    // Loop over all pairs and all modes
    for (int p = 0; p < nReplicas / 2; ++p) {
      for (int m = 0; m < nDimSingle; ++m) {
        DataType r = getRatio(g, p, m);

        // Metropolis check
        double prob = std::abs(r);
        double rand_val = rd->rdUniform01();
        total++;

        if (rand_val < prob) {
          accepted++;
          localUpdate(g, p, m);
          signRatio *= (r / prob);
        }
      }
    }

    return signRatio;
  }

  void getGreensMat(MatType &g0) override {
    // Not used in this context
  }

  // this is not correct yet
  DataType getSignOfWeight() override { return 1.0; }
};

// Trivial Swap Operator: Identity operator (no swaps)
// Used to verify multi-replica implementation matches single replica
class TrivialSwapOperator : public MixingOperator {
public:
  TrivialSwapOperator(int _nReplicas, int _nDimSingle, rdGenerator *_rd)
      : MixingOperator(_nReplicas, _nDimSingle, _rd) {
    // Initialize sigma to all 1 (Identity)
    // B(1) = 1/sqrt(2) * [1, 1; -1, 1] ... Wait.
    // Trivial operator should be Identity.
    // But our B(s) is never Identity!
    // B(s) is always a rotation by pi/4.
    // The "TrivialSwapOperator" was intended to be Identity.
    // If we want Identity, we need a different class or different logic.
    // But MixingOperator enforces the Sqrt-Swap form.
    // So TrivialSwapOperator cannot inherit from MixingOperator if it wants to
    // be Identity unless we change left_multiply to support Identity. But
    // MixingOperator is specifically for the Swap operator. So
    // TrivialSwapOperator is probably not applicable here anymore if we enforce
    // the form. Let's just make it do nothing in multiply?
  }

  // Override multiply to do nothing
  void left_multiply(const MatType &A, MatType &B) override { B = A; }
  void right_multiply(const MatType &A, MatType &B) override { B = A; }
  void inv_left_multiply(const MatType &A, MatType &B) override { B = A; }
  void inv_right_multiply(const MatType &A, MatType &B) override { B = A; }
  void stabilizedLeftMultiply(UDT &F) override {
    // Do nothing
  }

  DataType update(MatType &g) override { return 1.0; }
};

#endif // MIXING_OPERATOR_H
