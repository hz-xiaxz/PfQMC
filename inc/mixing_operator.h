#ifndef MIXING_OPERATOR_H
#define MIXING_OPERATOR_H

#include "operator.h"
#include "types.h"

/**
 * MixingOperator: Implements inter-replica mixing at τ=0 for 4-replica PFQMC.
 *
 * The mixing operator has a 4*nDimSingle × 4*nDimSingle block structure:
 *
 *   M_mix = | M_12    0   |
 *           |   0   M_34  |
 *
 * Where:
 * - M_12 (2*nDimSingle × 2*nDimSingle): mixes replicas 0 and 1
 * - M_34 (2*nDimSingle × 2*nDimSingle): mixes replicas 2 and 3
 *
 * This operator is applied only at τ=0 during sweeps to couple replica pairs.
 */
class MixingOperator : public Operator {
public:
    const int nDimSingle;   // single-replica dimension
    const int nDimTotal;    // = 4 * nDimSingle

    // M_12 block (2*nDimSingle × 2*nDimSingle) - mixes replicas 0 and 1
    MatType M_12;
    MatType M_12_inv;

    // M_34 block (2*nDimSingle × 2*nDimSingle) - mixes replicas 2 and 3
    MatType M_34;
    MatType M_34_inv;

    // Full mixing matrix (4*nDimSingle × 4*nDimSingle)
    MatType mat;
    MatType mat_inv;

    // Green's function for this operator: g0 = 2(1+M)^{-1} - 1
    MatType g0;
    MatType g0_inv;
    DataType signOfWeight;
    DataType signPf_g0_inv;

    /**
     * Constructor for MixingOperator.
     *
     * @param _nDimSingle Single-replica Hilbert space dimension
     * @param _M12 Mixing matrix for replica pair (0,1), size 2*nDimSingle × 2*nDimSingle
     * @param _M34 Mixing matrix for replica pair (2,3), size 2*nDimSingle × 2*nDimSingle
     */
    MixingOperator(int _nDimSingle, const MatType& _M12, const MatType& _M34)
        : nDimSingle(_nDimSingle), nDimTotal(4 * _nDimSingle)
    {
        assert(_M12.rows() == 2 * nDimSingle && _M12.cols() == 2 * nDimSingle);
        assert(_M34.rows() == 2 * nDimSingle && _M34.cols() == 2 * nDimSingle);

        M_12 = _M12;
        M_34 = _M34;
        M_12_inv = _M12.inverse();
        M_34_inv = _M34.inverse();

        // Build full block-diagonal matrix
        mat = MatType::Zero(nDimTotal, nDimTotal);
        mat.block(0, 0, 2 * nDimSingle, 2 * nDimSingle) = M_12;
        mat.block(2 * nDimSingle, 2 * nDimSingle, 2 * nDimSingle, 2 * nDimSingle) = M_34;

        mat_inv = MatType::Zero(nDimTotal, nDimTotal);
        mat_inv.block(0, 0, 2 * nDimSingle, 2 * nDimSingle) = M_12_inv;
        mat_inv.block(2 * nDimSingle, 2 * nDimSingle, 2 * nDimSingle, 2 * nDimSingle) = M_34_inv;

        // Compute Green's function: g0 = 2(1+M)^{-1} - 1
        MatType identity = MatType::Identity(nDimTotal, nDimTotal);
        g0 = 2.0 * (identity + mat).inverse() - identity;
        g0_inv = g0.inverse();

        // Compute sign of Pfaffian
        MatType tmp = g0_inv;
        signPf_g0_inv = signOfPfaf(tmp);
        signOfWeight = 1.0;  // Mixing operator has unit weight
    }

    /**
     * Factory method: Create identity mixing operator (no inter-replica coupling).
     */
    static MixingOperator* createIdentity(int nDimSingle) {
        MatType I_pair = MatType::Identity(2 * nDimSingle, 2 * nDimSingle);
        return new MixingOperator(nDimSingle, I_pair, I_pair);
    }

    // ========== Operator interface implementation ==========

    void left_multiply(const MatType& A, MatType& B) override {
        B = mat * A;
    }

    void inv_left_multiply(const MatType& A, MatType& B) override {
        B = mat_inv * A;
    }

    void adjoint_left_multiply(const MatType& A, MatType& B) override {
        B = mat.adjoint() * A;
    }

    void right_multiply(const MatType& A, MatType& B) override {
        B = A * mat;
    }

    void inv_right_multiply(const MatType& A, MatType& B) override {
        B = A * mat_inv;
    }

    void adjoint_inv_right_multiply(const MatType& A, MatType& B) override {
        B = A * mat_inv.adjoint();
    }

    /**
     * Left propagate: g → M * g * M^{-1}
     */
    void left_propagate(MatType& g, MatType& tmp) override {
        tmp = mat * g;
        g = tmp * mat_inv;
    }

    /**
     * Right propagate: g → M^{-1} * g * M
     */
    void right_propagate(MatType& g, MatType& tmp) override {
        tmp = mat_inv * g;
        g = tmp * mat;
    }

    /**
     * No auxiliary field updates for mixing operator - returns 1.0.
     */
    DataType update(MatType& g) override {
        return 1.0;
    }

    DataType getSignOfWeight() override {
        return signOfWeight;
    }

    void getGreensMat(MatType& g) override {
        g = g0;
    }

    void getGreensMatInv(MatType& g) override {
        g = g0_inv;
    }

    void stabilizedLeftMultiply(UDT& F) override {
        F = mat * F;
    }

    ~MixingOperator() override = default;
};

#endif // MIXING_OPERATOR_H
