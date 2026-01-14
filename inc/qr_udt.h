#ifndef QR_UDT_H
#define QR_UDT_H

#include "types.h"
#include <vector>

class UDT
{
public:
    int nDim;
    MatType U;
    dVecType D;
    MatType T;

    // Block support
    bool isBlock = false;
    std::vector<UDT> blocks;

    UDT() = default;
    UDT(int nDim)
    {
        this->nDim = nDim;
        U = MatType::Identity(nDim, nDim);
        D = dVecType::Ones(nDim);
        T = MatType::Identity(nDim, nDim);
    }

    UDT(const MatType &_U, const dVecType &_D, const MatType &_T)
    {
        this->nDim = _U.rows();
        U = _U;
        D = _D;
        T = _T;
    }

    UDT &operator=(const UDT &other)
    {
        nDim = other.nDim;
        U = other.U;
        D = other.D;
        T = other.T;
        isBlock = other.isBlock;
        blocks = other.blocks;
        return *this;
    }

    UDT(const UDT &other)
    {
        *this = other;
    }

    UDT &operator=(UDT &&other)
    {
        if (this != &other)
        {
            nDim = other.nDim;
            U = std::move(other.U);
            D = std::move(other.D);
            T = std::move(other.T);
            isBlock = other.isBlock;
            blocks = std::move(other.blocks);
        }
        return (*this);
    }

    UDT(UDT &&other)
    {
        *this = std::move(other);
    }

    // use qr to get UDT decomposition
    // nBlocks > 1 triggers block decomposition
    explicit UDT(MatType &A, int nBlocks = 1)
    {
        nDim = A.rows();
        
        if (nBlocks > 1 && nDim % nBlocks == 0) {
            isBlock = true;
            int blkSize = nDim / nBlocks;
            blocks.resize(nBlocks);
            
            // Decompose each block
            for (int i = 0; i < nBlocks; i++) {
                MatType blk = A.block(i * blkSize, i * blkSize, blkSize, blkSize);
                blocks[i] = UDT(blk, 1); 
            }
            return;
        }

        T = MatType::Zero(nDim, nDim);
        D = dVecType(nDim);
        int jpvt[nDim];
        DataType tau[nDim];
        for (int j = 0; j < nDim; j++)
        {
            jpvt[j] = 0;
        }

        LAPACKE_zgeqp3(LAPACK_COL_MAJOR, nDim, nDim, A.data(), nDim, jpvt, tau);

        double alpha;
        for (int i = 0; i < nDim; i++)
        {
            D(i) = std::abs(A(i, i).real()); // A(i, i)'s are all real
            alpha = 1.0 / D(i);
            for (int j = i; j < nDim; j++)
            {
                A(i, j) = A(i, j) * alpha;
            }
        }

        int j;
        for (int i = 0; i < nDim; i++)
        {
            j = jpvt[i] - 1;
            T(i, j) = 1;
        }

        cblas_ztrmm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, nDim, nDim, &one, A.data(), nDim, T.data(), nDim);

        LAPACKE_zungqr(LAPACK_COL_MAJOR, nDim, nDim, nDim, A.data(), nDim, tau);

        U = A;
    }

    // g = 2 * (1 + UDT)^{-1}
    inline void onePlusInv(MatType &g) const
    {
        if (isBlock) {
            g = MatType::Zero(nDim, nDim);
            int blkSize = nDim / blocks.size();
            for (size_t i = 0; i < blocks.size(); i++) {
                MatType gBlk;
                blocks[i].onePlusInv(gBlk);
                g.block(i * blkSize, i * blkSize, blkSize, blkSize) = gBlk;
            }
            return;
        }

        MatType Xinv = T.inverse();
        dVecType Dpinv(nDim);
        dVecType Dm(nDim);
        for (int i = 0; i < nDim; i++)
        {
            Dpinv(i) = 1.0 / std::max(D(i), 1.0);
            Dm(i) = std::min(D(i), 1.0);
        }
        MatType tmp1 = Xinv * Dpinv.asDiagonal();
        MatType tmp2 = (tmp1 + U * Dm.asDiagonal()).inverse();
        g = 2.0 * tmp1 * tmp2;
    }
};

inline UDT operator*(const UDT &udtL, const UDT &udtR)
{
    if (udtL.isBlock && udtR.isBlock) {
        UDT res;
        res.nDim = udtL.nDim;
        res.isBlock = true;
        res.blocks.resize(udtL.blocks.size());
        for (size_t i = 0; i < udtL.blocks.size(); i++) {
            res.blocks[i] = udtL.blocks[i] * udtR.blocks[i];
        }
        return res;
    }

    MatType mat = udtL.T * udtR.U;
    mat = udtL.D.asDiagonal() * mat;
    mat = mat * udtR.D.asDiagonal();
    UDT tmp(mat);
    tmp.U = udtL.U * tmp.U;
    tmp.T = tmp.T * udtR.T;
    return tmp;
}

// TO DO: lazy evaluation, check udta = B * udta
inline UDT operator*(const MatType &B, const UDT &udtR)
{
    if (udtR.isBlock) {
        UDT res;
        res.nDim = udtR.nDim;
        res.isBlock = true;
        res.blocks.resize(udtR.blocks.size());
        int blkSize = udtR.nDim / udtR.blocks.size();
        
        for (size_t i = 0; i < udtR.blocks.size(); i++) {
            MatType Bblk = B.block(i * blkSize, i * blkSize, blkSize, blkSize);
            res.blocks[i] = Bblk * udtR.blocks[i];
        }
        return res;
    }

    MatType mat = B * udtR.U;
    mat = mat * udtR.D.asDiagonal();
    UDT tmp(mat);
    tmp.T = tmp.T * udtR.T;
    return tmp;
}

// return (1+udtR@(udtL).adjoint)^{-1}
inline MatType onePlusInv(UDT &udtL, UDT &udtR)
{
    if (udtL.isBlock && udtR.isBlock) {
        int nDim = udtL.nDim;
        MatType res = MatType::Zero(nDim, nDim);
        int blkSize = nDim / udtL.blocks.size();
        for (size_t i = 0; i < udtL.blocks.size(); i++) {
            MatType blk = onePlusInv(udtL.blocks[i], udtR.blocks[i]);
            res.block(i * blkSize, i * blkSize, blkSize, blkSize) = blk;
        }
        return res;
    }

    int n = udtR.U.cols();
    MatType tem1 = udtR.U.adjoint() * udtL.U;
    MatType tem2 = udtR.T * udtL.T.adjoint();
    auto DrPinv = dVecType(n);
    auto DrM = dVecType(n);
    auto DlPinv = dVecType(n);
    auto DlM = dVecType(n);
    for (int i = 0; i < n; i++)
    {
        if (std::abs(udtR.D(i)) > 1.0)
        {
            DrPinv(i) = 1.0 / udtR.D(i);
            DrM(i) = 1.0;
        }
        else
        {
            DrPinv(i) = 1.0;
            DrM(i) = udtR.D(i);
        }
        if (std::abs(udtL.D(i)) > 1.0)
        {
            DlPinv(i) = 1.0 / udtL.D(i);
            DlM(i) = 1.0;
        }
        else
        {
            DlPinv(i) = 1.0;
            DlM(i) = udtL.D(i);
        }
    }

    tem2 = DrPinv.asDiagonal() * tem1 * DlPinv.asDiagonal() + DrM.asDiagonal() * tem2 * DlM.asDiagonal();
    tem1 = DlPinv.asDiagonal() * tem2.inverse() * DrPinv.asDiagonal();

    return 2. * udtL.U * tem1 * (udtR.U.adjoint());
}
#endif