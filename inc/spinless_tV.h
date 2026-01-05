#ifndef Spinless_tV_H
#define Spinless_tV_H

#include "operator.h"
#include "types.h"

class SpinlessTvUtils {
   public:
    int Lx, Ly;
    // int nsites;
    int nDim;
    // std::vector<int> nBond;
    double dt;
    double V;
    int l;
    double lambdaV, chlV, shlV, thlV, etaM;

    const int hsScheme;  

    const bool singleMaj;

    SpinlessTvUtils(int _Lx, int _Ly, double _dt, double _V, int _l, int _nDim,
                    bool _singleMaj = false, int _hsScheme = 0)
        : singleMaj(_singleMaj), hsScheme(_hsScheme) {
        // model configuration
        Lx = _Lx;
        Ly = _Ly;
        dt = _dt;
        V = _V;
        l = _l;  // imaginary time slices
        nDim = _nDim;

        lambdaV = acosh(exp(0.5 * V * dt));
        chlV = cosh(lambdaV);
        shlV = sinh(lambdaV);
        thlV = tanh(lambdaV);
        etaM = chlV * chlV;
    }

    inline int unitCellCoord2Idx(int ix, int iy) const { return ix * Ly + iy; }

    // virtual inline void KineticGenerator(MatType &H, DataType t) const {};

    // virtual inline DataType energyFromGreensFunc(const MatType &g) {};
    virtual inline void aux2MajoranaIdx(int idxAux, int imaj, int bType,
                                        int &idx1, int &idx2) const = 0;

    // Generate the Greens function for single slice
    inline void InteractionTanhGenerator(MatType &H, const iVecType &s,
                                         const int bondType,
                                         bool inv = false) const {
        DataType tmp = (1.0i) * tanh(0.5 * lambdaV);
        if (inv) {
            tmp = (-1.0) / tmp;
        }

        if (hsScheme == 0) {
            int idx1, idx2;
            if (!singleMaj) {
                for (int i = 0; i < s.size(); i++) {
                    for (int k = 0; k < 2; k++) {
                        aux2MajoranaIdx(i, k, bondType, idx1, idx2);
                        H(idx1, idx2) += -tmp * double(s(i));
                        H(idx2, idx1) += +tmp * double(s(i));
                    }
                }
            } else {
                for (int i = 0; i < s.size(); i++) {
                    aux2MajoranaIdx(i, 0, bondType, idx1, idx2);
                    H(idx1, idx2) += -tmp * double(s(i));
                    H(idx2, idx1) += +tmp * double(s(i));
                }
            }
        } else if (hsScheme == 1) {
            int idxi1, idxi2, idxj1, idxj2;
            assert(!singleMaj);
            for (int i = 0; i < s.size(); i++) {
                aux2MajoranaIdx(i, 0, bondType, idxi1, idxj1);
                aux2MajoranaIdx(i, 1, bondType, idxi2, idxj2);
                H(idxi1, idxi2) += -tmp * double(s(i));
                H(idxi2, idxi1) += +tmp * double(s(i));
                H(idxj1, idxj2) += +tmp * double(s(i));
                H(idxj2, idxj1) += -tmp * double(s(i));
            }
        }
    }

    // Directly generate B by directly writing each 2*2 block
    // B should be initialized as Identity
    inline void InteractionBGenerator(MatType &B, const iVecType &s,
                                      const int bondType,
                                      bool inv = false) const {
        DataType ch = chlV;
        DataType ish = (1.0i) * shlV;
        if (inv) ish = -ish;

        if (hsScheme == 0) {
            // B = exp(λ/2 * s * (i γi1 γj1 + i γi2 γj2))
            int idx1, idx2;
            if (!singleMaj) {
                for (int i = 0; i < s.size(); i++) {
                    for (int k = 0; k < 2; k++) {
                        aux2MajoranaIdx(i, k, bondType, idx1, idx2);
                        B(idx1, idx1) = ch;
                        B(idx2, idx2) = ch;
                        B(idx1, idx2) = +ish * double(s(i));
                        B(idx2, idx1) = -ish * double(s(i));
                    }
                }
            } else {
                for (int i = 0; i < s.size(); i++) {
                    aux2MajoranaIdx(i, 0, bondType, idx1, idx2);
                    B(idx1, idx1) = ch;
                    B(idx2, idx2) = ch;
                    B(idx1, idx2) = +ish * double(s(i));
                    B(idx2, idx1) = -ish * double(s(i));
                }
            }
        } else if (hsScheme == 1) {
            // B = exp(λ/2 * s * (i γi1 γi2 - i γj1 γj2))
            int idxi1, idxi2, idxj1, idxj2;
            assert(!singleMaj);
            for (int i = 0; i < s.size(); i++) {
                aux2MajoranaIdx(i, 0, bondType, idxi1, idxj1);
                aux2MajoranaIdx(i, 1, bondType, idxi2, idxj2);
                B(idxi1, idxi1) = ch;
                B(idxi2, idxi2) = ch;
                B(idxi1, idxi2) = +ish * double(s(i));
                B(idxi2, idxi1) = -ish * double(s(i));

                B(idxj1, idxj1) = ch;
                B(idxj2, idxj2) = ch;
                B(idxj1, idxj2) = -ish * double(s(i));
                B(idxj2, idxj1) = +ish * double(s(i));
            }
        }
    }
};

class SpinlessVOperator : public Operator {
   protected:
    const SpinlessTvUtils *config;
    const double etaM;
    // delayed update for spinless t-V requires additional diagonalization
    // const int delay_max = 32;
    // cannot be larger than 64 (unless you change the threads value in
    // kernel_linalg.cpp)
   public:
    // const int nUnitcell;
    const int bondType;

    // ===== Replica support =====
    // ===== Replica support =====
    const int nReplicas;    // Number of replicas (default 1 for backward compatibility)
    const int nDimSingle;   // Single-replica dimension
    const int nDim;         // Total dimension = nReplicas * nDimSingle

    const int hsScheme; // 0: hopping channel, 1: density channel

    // For multi-replica: vector of aux field pointers (one per replica)
    std::vector<iVecType*> s;

    // Per-replica B matrices (each nDimSingle × nDimSingle)
    std::vector<MatType> B_replica;

    // Full block-diagonal B matrix (nDim × nDim)
    MatType B;
    MatType B_inv;
    rdGenerator *rd;
    bool assumeBlockDiagonal = false; // Optimization flag

    const bool singleMaj;

    // Multi-replica constructor
    // _s: vector of aux field pointers, one per replica
    SpinlessVOperator(const SpinlessTvUtils *_config,
                      std::vector<iVecType*> _s,
                      int _bondType,
                      rdGenerator *_rd,
                      int _nReplicas = 1)
        : config(_config),
          etaM(_config->etaM),
          bondType(_bondType),
          nReplicas(_nReplicas),
          nDimSingle(_config->nDim),
          nDim(_nReplicas * _config->nDim),
          singleMaj(_config->singleMaj),
          hsScheme(_config->hsScheme) {
        assert(static_cast<int>(_s.size()) == _nReplicas);
        s = _s;
        rd = _rd;

        // Initialize per-replica B matrices
        B_replica.resize(nReplicas);
        for (int r = 0; r < nReplicas; r++) {
            B_replica[r] = MatType::Identity(nDimSingle, nDimSingle);
            config->InteractionBGenerator(B_replica[r], *s[r], bondType, false);
        }

        // Build full block-diagonal B
        rebuildFullB();
    }

    // Backward-compatible single-replica constructor
    // _s: aux fields, Z_2 variable, length = nUnitcell
    SpinlessVOperator(const SpinlessTvUtils *_config, iVecType *_s,
                      int _bondType, rdGenerator *_rd)
        : SpinlessVOperator(_config, std::vector<iVecType*>{_s}, _bondType, _rd, 1) {}

    ~SpinlessVOperator() {
        for (auto* ptr : s) {
            delete ptr;
        }
    }

    // Rebuild the full block-diagonal B matrix from per-replica matrices
    void rebuildFullB() {
        if (nReplicas == 1) {
            B = B_replica[0];
            return;
        }

        B = MatType::Zero(nDim, nDim);
        for (int r = 0; r < nReplicas; r++) {
            B.block(r * nDimSingle, r * nDimSingle, nDimSingle, nDimSingle) = B_replica[r];
        }
    }

    void reCalcInv() {
        if (nReplicas == 1) {
            B_inv = MatType::Identity(nDimSingle, nDimSingle);
            config->InteractionBGenerator(B_inv, *s[0], bondType, true);
            return;
        }

        B_inv = MatType::Zero(nDim, nDim);
        for (int r = 0; r < nReplicas; r++) {
            MatType B_inv_r = MatType::Identity(nDimSingle, nDimSingle);
            config->InteractionBGenerator(B_inv_r, *s[r], bondType, true);
            B_inv.block(r * nDimSingle, r * nDimSingle, nDimSingle, nDimSingle) = B_inv_r;
        }
    }

    void setAssumeBlockDiagonal(bool val) {
        assumeBlockDiagonal = val;
    }

    // virtual inline void aux2MajoranaIdx(int idxAux, int imaj, int& idx1, int&
    // idx2) {};

    void singleFlip(MatType &g, int replica, int idxAux, double rand, bool &flag,
                    DataType &signCur) {
        // Compute offset for this replica in the global matrix
        int offset = replica * nDimSingle;
        
        DataType r;
        // auto m = mConfig->idxCell2Coord(idxCell);
        int auxCur = (*s[replica])(idxAux);
        int idx1, idx2, idx3, idx4;
        DataType tmp[2];
        const int inc = 1;
        DataType alpha;

        config->aux2MajoranaIdx(idxAux, 0, bondType, idx1, idx2); // 1j, 1k
        config->aux2MajoranaIdx(idxAux, 1, bondType, idx3, idx4); // 2j, 2k
        
        // Convert to global indices
        int idx1_g = offset + idx1;
        int idx2_g = offset + idx2;
        int idx3_g = offset + idx3;
        int idx4_g = offset + idx4;

        if (hsScheme == 0) {
            tmp[0] =
                (1.0 - ((1.0i) * (config->thlV) * double(auxCur) * g(idx1_g, idx2_g)));
            tmp[1] =
                (1.0 - ((1.0i) * (config->thlV) * double(auxCur) * g(idx3_g, idx4_g)));
            r = tmp[0] * tmp[1];
            // std::cout << "r1" << r << " ";
            r += (config->thlV * config->thlV) * ((g(idx1_g, idx3_g) * g(idx2_g, idx4_g)) -
                                                (g(idx2_g, idx3_g) * g(idx1_g, idx4_g)));
            // std::cout << " r2" << r << "\n";
            r *= etaM;
        } else if (hsScheme == 1) {
            tmp[0] =
                (1.0 - ((1.0i) * (config->thlV) * double(auxCur) * g(idx1_g, idx3_g)));
            tmp[1] =
                (1.0 + ((1.0i) * (config->thlV) * double(auxCur) * g(idx2_g, idx4_g)));
            r = tmp[0] * tmp[1];
            // std::cout << "r1" << r << " ";
            r -= (config->thlV * config->thlV) * ((g(idx1_g, idx2_g) * g(idx3_g, idx4_g)) -
                                                (g(idx3_g, idx2_g) * g(idx1_g, idx4_g)));
            // std::cout << " r2" << r << "\n";
            r *= etaM;
        }
        // for (int imaj = 0; imaj < 2; imaj ++) {
        //     config->aux2MajoranaIdx(idxAux, imaj, bondType, idx1, idx2);
        //     // idx1 = mConfig->majoranaCoord2Idx(m.ix, m.iy, 0, imaj);
        //     // idx2 = mConfig->neighborSiteIdx(m.ix, m.iy, imaj, bondType);
        //     // tmp = [1 + i \sigma_{12} \tanh(\lambda / 2) G_{12}]
        //     tmp[imaj] = ( 1.0 - ( (1.0i) * (config->thlV) * double(auxCur) *
        //     g(idx1, idx2) ) ); r *= tmp[imaj];
        // }

        flag = rand < std::abs(r);
        // std::cout << rand << "=rand " << "r = " << r << "\n";

        if (flag) {
            // DataType t = (r / std::abs(r));
            // if (std::abs(t.imag()) > 1e-5) {
            //     std::cout << "r= " << r << " sign(r)= " << t << "\n";
            // }
            signCur *= (r / std::abs(r));

            if (hsScheme == 0) {
                for (int imaj = 0; imaj < 2; imaj++) {
                    config->aux2MajoranaIdx(idxAux, imaj, bondType, idx1, idx2);
                    int idx1_g = offset + idx1;
                    int idx2_g = offset + idx2;
                    
                    if (imaj == 1) {
                        tmp[1] = (1.0 - ((1.0i) * (config->thlV) * double(auxCur) *
                                        g(idx1_g, idx2_g)));
                    }
                    // update aux field and B matrix
                    (*s[replica])(idxAux) = -auxCur;
                    B_replica[replica](idx1, idx2) = -B_replica[replica](idx1, idx2);
                    B_replica[replica](idx2, idx1) = -B_replica[replica](idx2, idx1);

                    // update Green's function
                    if (assumeBlockDiagonal) {
                        // Optimized update for block diagonal case
                        // Only update the block corresponding to this replica
                        // Indices are already local to the block if we subtract offset?
                        // No, idx1_g is global. We need to map it to local or use block operations.
                        // Actually, since g is stored as one large matrix, we just need to limit the 
                        // rank-1 update to the columns/rows of this replica.
                        
                        // x1 and x2 are columns. We only need the segment [offset, offset+nDimSingle]
                        // But g.col() returns the whole column.
                        // We can manually implement the update loop over the relevant range.
                        
                        int startIdx = offset;
                        int endIdx = offset + nDimSingle;
                        int localDim = nDimSingle;
                        
                        // We need x1 and x2 segments.
                        // x1 = -g.col(idx1_g) + 2*e_{idx1_g}
                        // x2 = -g.col(idx2_g) + 2*e_{idx2_g}
                        
                        // Optimized update for block diagonal case
                        // Only update the block corresponding to this replica
                        
                        // We need x1 and x2 segments.
                        // x1 = -g.col(idx1_g) + 2*e_{idx1_g}
                        // x2 = -g.col(idx2_g) + 2*e_{idx2_g}
                        
                        // Extract relevant segments of columns to local vectors
                        cVecType x1_seg = -g.col(idx1_g).segment(offset, nDimSingle);
                        cVecType x2_seg = -g.col(idx2_g).segment(offset, nDimSingle);
                        
                        x1_seg(idx1) += 2.0; // idx1 is local index
                        x2_seg(idx2) += 2.0; // idx2 is local index
                        
                        alpha = (+1.0i) * double(auxCur) * (config->thlV) / tmp[imaj];
                        
                        // Rank-1 update on the block: G_block += alpha * x1 * x2^T
                        // We use zgeru. The block starts at (offset, offset).
                        // The leading dimension (stride between columns) is nDim (total dimension).
                        DataType* g_block_ptr = g.data() + offset + offset * nDim;
                        
                        zgeru(&nDimSingle, &nDimSingle, reinterpret_cast<MKL_Complex16*>(&alpha), 
                              reinterpret_cast<MKL_Complex16*>(x1_seg.data()), &inc, 
                              reinterpret_cast<MKL_Complex16*>(x2_seg.data()), &inc,
                              reinterpret_cast<MKL_Complex16*>(g_block_ptr), &nDim);
                        
                        alpha = -alpha;
                        // Rank-1 update: G_block += -alpha * x2 * x1^T
                        zgeru(&nDimSingle, &nDimSingle, reinterpret_cast<MKL_Complex16*>(&alpha), 
                              reinterpret_cast<MKL_Complex16*>(x2_seg.data()), &inc, 
                              reinterpret_cast<MKL_Complex16*>(x1_seg.data()), &inc,
                              reinterpret_cast<MKL_Complex16*>(g_block_ptr), &nDim);
                            
                    } else {
                        // Full dense update (original)
                        cVecType x1 = -g.col(idx1_g);
                        cVecType x2 = -g.col(idx2_g);
                        x1(idx1_g) += 2;
                        x2(idx2_g) += 2;
                        alpha = (+1.0i) * double(auxCur) * (config->thlV) / tmp[imaj];
                        zgeru(&nDim, &nDim, &alpha, x1.data(), &inc, x2.data(), &inc,
                            g.data(), &nDim);
                        alpha = -alpha;
                        zgeru(&nDim, &nDim, &alpha, x2.data(), &inc, x1.data(), &inc,
                            g.data(), &nDim);
                    }
                }
            } else if (hsScheme == 1) {
                int idxj1, idxk1, idxj2, idxk2;
                config->aux2MajoranaIdx(idxAux, 0, bondType, idxj1, idxk1);
                config->aux2MajoranaIdx(idxAux, 1, bondType, idxj2, idxk2);
                for (int iaux = 0; iaux < 2; iaux++) {
                    if (iaux == 0) {
                        idx1 = idxj1;
                        idx2 = idxj2;
                    } else {
                        idx1 = idxk1;
                        idx2 = idxk2;
                        tmp[1] = (1.0 + ((1.0i) * (config->thlV) * double(auxCur) *
                                        g(idx1_g, idx2_g)));
                    }
                    int idx1_g = offset + idx1;
                    int idx2_g = offset + idx2;

                    // update aux field and B matrix
                    (*s[replica])(idxAux) = -auxCur;
                    B_replica[replica](idx1, idx2) = -B_replica[replica](idx1, idx2);
                    B_replica[replica](idx2, idx1) = -B_replica[replica](idx2, idx1);

                    // update Green's function
                    cVecType x1 = -g.col(idx1_g);
                    cVecType x2 = -g.col(idx2_g);
                    x1(idx1_g) += 2;
                    x2(idx2_g) += 2;
                    alpha = (+1.0i) * double(auxCur) * (config->thlV) / tmp[iaux];
                    if (iaux == 1) alpha = -alpha;
                    zgeru(&nDim, &nDim, &alpha, x1.data(), &inc, x2.data(), &inc,
                        g.data(), &nDim);
                    alpha = -alpha;
                    zgeru(&nDim, &nDim, &alpha, x2.data(), &inc, x1.data(), &inc,
                        g.data(), &nDim);
                }
            }
        }
    };

    void singleFlipSingleMajorana(MatType &g, int replica, int idxAux, double rand,
                                  bool &flag, DataType &signCur) {
        int offset = replica * nDimSingle;
        DataType r;
        // auto m = mConfig->idxCell2Coord(idxCell);
        int auxCur = (*s[replica])(idxAux);
        int idx1, idx2;
        DataType tmp;
        const int inc = 1;
        DataType alpha;

        config->aux2MajoranaIdx(idxAux, 0, bondType, idx1, idx2);
        int idx1_g = offset + idx1;
        int idx2_g = offset + idx2;
        
        tmp =
            (1.0 - ((1.0i) * (config->thlV) * double(auxCur) * g(idx1_g, idx2_g)));
        r = tmp * tmp * etaM;

        flag = rand < std::abs(r);
        // std::cout << rand << "=rand " << "r = " << r << "\n";

        if (flag) {
            // std::cout << "tmp = " << tmp << " r = " << r << "\n";
            signCur *= (tmp / std::abs(tmp));
            config->aux2MajoranaIdx(idxAux, 0, bondType, idx1, idx2);
            // update aux field and B matrix
            (*s[replica])(idxAux) = -auxCur;
            B_replica[replica](idx1, idx2) = -B_replica[replica](idx1, idx2);
            B_replica[replica](idx2, idx1) = -B_replica[replica](idx2, idx1);

            // update Green's function
            cVecType x1 = -g.col(idx1_g);
            cVecType x2 = -g.col(idx2_g);
            x1(idx1_g) += 2;
            x2(idx2_g) += 2;
            alpha = (+1.0i) * double(auxCur) * (config->thlV) / tmp;
            zgeru(&nDim, &nDim, &alpha, x1.data(), &inc, x2.data(), &inc,
                  g.data(), &nDim);
            alpha = -alpha;
            zgeru(&nDim, &nDim, &alpha, x2.data(), &inc, x1.data(), &inc,
                  g.data(), &nDim);
        }
    };

    DataType update(MatType &g) override {
        double rand;
        bool flag;
        DataType signCur = 1.0;
        
        // Loop over all replicas
        for (int r = 0; r < nReplicas; r++) {
            if (singleMaj) {
                for (int i = 0; i < s[r]->size(); i++) {
                    // random real between (0, 1)
                    rand = rd->rdUniform01();
                    singleFlipSingleMajorana(g, r, i, rand, flag, signCur);
                }
            } else {
                for (int i = 0; i < s[r]->size(); i++) {
                    // random real between (0, 1)
                    rand = rd->rdUniform01();
                    singleFlip(g, r, i, rand, flag, signCur);
                }
            }
        }
        return signCur;
    }

    void right_multiply(const MatType &AIn, MatType &AOut) override {
        AOut = AIn * B;
    }

    void left_multiply(const MatType &AIn, MatType &Aout) override {
        Aout = B * AIn;
    }

    void left_propagate(MatType &g, MatType &gTmp) override {
        reCalcInv();
        gTmp = B * g;
        g = gTmp * B_inv;
    }
    void right_propagate(MatType &g, MatType &gTmp) override {
        reCalcInv();
        gTmp = B_inv * g;
        g = gTmp * B;
    }

    // Return the first replica's aux field for backward compatibility with single-replica code/tests
    iVecType *getAuxField() override { return s[0]; }

    int getType() override { return bondType; }

    void getGreensMat(MatType &g) override {
        if (nReplicas == 1) {
            g = MatType::Zero(nDim, nDim);
            config->InteractionTanhGenerator(g, *s[0], bondType, false);
            return;
        }

        g = MatType::Zero(nDim, nDim);
        for (int r = 0; r < nReplicas; r++) {
            MatType g_r = MatType::Zero(nDimSingle, nDimSingle);
            config->InteractionTanhGenerator(g_r, *s[r], bondType, false);
            g.block(r * nDimSingle, r * nDimSingle, nDimSingle, nDimSingle) = g_r;
        }
    }

    void getGreensMatInv(MatType &g) override {
        if (nReplicas == 1) {
            g = MatType::Zero(nDim, nDim);
            config->InteractionTanhGenerator(g, *s[0], bondType, true);
            return;
        }

        g = MatType::Zero(nDim, nDim);
        for (int r = 0; r < nReplicas; r++) {
            MatType g_r = MatType::Zero(nDimSingle, nDimSingle);
            config->InteractionTanhGenerator(g_r, *s[r], bondType, true);
            g.block(r * nDimSingle, r * nDimSingle, nDimSingle, nDimSingle) = g_r;
        }
    }

    // inline DataType getSignPfGInv() override;

    void stabilizedLeftMultiply(UDT &F) override {
        // F.bMultUpdate(B);
        F = B * F;
    }
};

class Spinless_tV {
   public:
    std::vector<Operator *> op_array;
    int nDim;
    int nDimSingle; // NEW: single replica dimension
    ~Spinless_tV() {
        for (int i = 0; i < op_array.size(); i++) {
            delete op_array[i];
        }
    }
};

#endif