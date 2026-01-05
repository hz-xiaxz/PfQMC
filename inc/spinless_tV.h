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

    DataType getRatio(MatType &g, int replica, int idxAux) {
        int offset = replica * nDimSingle;
        int auxCur = (*s[replica])(idxAux);
        int idx1, idx2, idx3, idx4;
        DataType tmp[2];
        DataType r;

        config->aux2MajoranaIdx(idxAux, 0, bondType, idx1, idx2);
        config->aux2MajoranaIdx(idxAux, 1, bondType, idx3, idx4);

        int idx1_g = offset + idx1;
        int idx2_g = offset + idx2;
        int idx3_g = offset + idx3;
        int idx4_g = offset + idx4;

        if (hsScheme == 0) {
            tmp[0] = (1.0 - ((1.0i) * (config->thlV) * double(auxCur) * g(idx1_g, idx2_g)));
            tmp[1] = (1.0 - ((1.0i) * (config->thlV) * double(auxCur) * g(idx3_g, idx4_g)));
            r = tmp[0] * tmp[1];
            r += (config->thlV * config->thlV) * ((g(idx1_g, idx3_g) * g(idx2_g, idx4_g)) -
                                                (g(idx2_g, idx3_g) * g(idx1_g, idx4_g)));
            r *= etaM;
        } else { // hsScheme == 1
            tmp[0] = (1.0 - ((1.0i) * (config->thlV) * double(auxCur) * g(idx1_g, idx3_g)));
            tmp[1] = (1.0 + ((1.0i) * (config->thlV) * double(auxCur) * g(idx2_g, idx4_g)));
            r = tmp[0] * tmp[1];
            r -= (config->thlV * config->thlV) * ((g(idx1_g, idx2_g) * g(idx3_g, idx4_g)) -
                                                (g(idx3_g, idx2_g) * g(idx1_g, idx4_g)));
            r *= etaM;
        }
        return r;
    }

    void updateReplica(MatType &g, int replica, int idxAux) {
        int offset = replica * nDimSingle;
        int auxCur = (*s[replica])(idxAux);
        const int inc = 1;
        
        if (hsScheme == 0) {
            for (int imaj = 0; imaj < 2; imaj++) {
                int idx1, idx2;
                config->aux2MajoranaIdx(idxAux, imaj, bondType, idx1, idx2);
                int idx1_g = offset + idx1;
                int idx2_g = offset + idx2;
                
                DataType tmp_update = (1.0 - ((1.0i) * (config->thlV) * double(auxCur) * g(idx1_g, idx2_g)));

                (*s[replica])(idxAux) = -auxCur;
                B_replica[replica](idx1, idx2) = -B_replica[replica](idx1, idx2);
                B_replica[replica](idx2, idx1) = -B_replica[replica](idx2, idx1);

                cVecType x1 = -g.col(idx1_g);
                cVecType x2 = -g.col(idx2_g);
                x1(idx1_g) += 2;
                x2(idx2_g) += 2;
                DataType alpha = (+1.0i) * double(auxCur) * (config->thlV) / tmp_update;
                
                zgeru(&nDim, &nDim, reinterpret_cast<MKL_Complex16*>(&alpha), 
                    reinterpret_cast<MKL_Complex16*>(x1.data()), &inc, 
                    reinterpret_cast<MKL_Complex16*>(x2.data()), &inc, 
                    reinterpret_cast<MKL_Complex16*>(g.data()), &nDim);
                alpha = -alpha;
                zgeru(&nDim, &nDim, reinterpret_cast<MKL_Complex16*>(&alpha), 
                    reinterpret_cast<MKL_Complex16*>(x2.data()), &inc, 
                    reinterpret_cast<MKL_Complex16*>(x1.data()), &inc, 
                    reinterpret_cast<MKL_Complex16*>(g.data()), &nDim);
            }
        } else { // hsScheme == 1
            int idxj1, idxk1, idxj2, idxk2;
            config->aux2MajoranaIdx(idxAux, 0, bondType, idxj1, idxk1);
            config->aux2MajoranaIdx(idxAux, 1, bondType, idxj2, idxk2);
            for (int iaux = 0; iaux < 2; iaux++) {
                int idx1, idx2;
                DataType tmp_update;
                 if (iaux == 0) {
                    idx1 = idxj1; idx2 = idxj2;
                    tmp_update = (1.0 - ((1.0i) * (config->thlV) * double(auxCur) * g(offset + idx1, offset + idx2)));
                } else {
                    idx1 = idxk1; idx2 = idxk2;
                    tmp_update = (1.0 + ((1.0i) * (config->thlV) * double(auxCur) * g(offset + idx1, offset + idx2)));
                }
                int idx1_g = offset + idx1;
                int idx2_g = offset + idx2;

                (*s[replica])(idxAux) = -auxCur;
                B_replica[replica](idx1, idx2) = -B_replica[replica](idx1, idx2);
                B_replica[replica](idx2, idx1) = -B_replica[replica](idx2, idx1);

                cVecType x1 = -g.col(idx1_g);
                cVecType x2 = -g.col(idx2_g);
                x1(idx1_g) += 2;
                x2(idx2_g) += 2;
                DataType alpha = (+1.0i) * double(auxCur) * (config->thlV) / tmp_update;
                if (iaux == 1) alpha = -alpha;
                
                zgeru(&nDim, &nDim, reinterpret_cast<MKL_Complex16*>(&alpha), 
                    reinterpret_cast<MKL_Complex16*>(x1.data()), &inc, 
                    reinterpret_cast<MKL_Complex16*>(x2.data()), &inc, 
                    reinterpret_cast<MKL_Complex16*>(g.data()), &nDim);
                alpha = -alpha;
                zgeru(&nDim, &nDim, reinterpret_cast<MKL_Complex16*>(&alpha), 
                    reinterpret_cast<MKL_Complex16*>(x2.data()), &inc, 
                    reinterpret_cast<MKL_Complex16*>(x1.data()), &inc, 
                    reinterpret_cast<MKL_Complex16*>(g.data()), &nDim);
            }
        }
    }

    void updateBlockPair(MatType &g, int rA, int rB, int idxAux) {
        // Optimized update for block diagonal case (pair block)
        // We assume rA and rB form a pair (0,1) or (2,3) etc.
        // The block offset is determined by the pair index.
        int p = rA / 2; // Assuming rA is even and rB is rA+1
        int blkOffset = p * 2 * nDimSingle;
        int blkSize = 2 * nDimSingle;
        DataType* g_block_ptr = g.data() + blkOffset + blkOffset * nDim;
        const int inc = 1;

        for (int replica : {rA, rB}) {
            int offset = replica * nDimSingle;
            int auxCur = (*s[replica])(idxAux);

            if (hsScheme == 0) {
                for (int imaj = 0; imaj < 2; imaj++) {
                    int idx1, idx2;
                    config->aux2MajoranaIdx(idxAux, imaj, bondType, idx1, idx2);
                    int idx1_g = offset + idx1;
                    int idx2_g = offset + idx2;
                    
                    DataType tmp_update = (1.0 - ((1.0i) * (config->thlV) * double(auxCur) * g(idx1_g, idx2_g)));

                    (*s[replica])(idxAux) = -auxCur;
                    B_replica[replica](idx1, idx2) = -B_replica[replica](idx1, idx2);
                    B_replica[replica](idx2, idx1) = -B_replica[replica](idx2, idx1);

                    DataType alpha = (+1.0i) * double(auxCur) * (config->thlV) / tmp_update;

                    // Extract relevant segments of columns to local vectors
                    cVecType x1_seg = -g.col(idx1_g).segment(blkOffset, blkSize);
                    cVecType x2_seg = -g.col(idx2_g).segment(blkOffset, blkSize);
                    
                    int idx1_local = idx1_g - blkOffset;
                    int idx2_local = idx2_g - blkOffset;
                    
                    x1_seg(idx1_local) += 2.0; 
                    x2_seg(idx2_local) += 2.0; 
                    
                    zgeru(&blkSize, &blkSize, reinterpret_cast<MKL_Complex16*>(&alpha), 
                        reinterpret_cast<MKL_Complex16*>(x1_seg.data()), &inc, 
                        reinterpret_cast<MKL_Complex16*>(x2_seg.data()), &inc,
                        reinterpret_cast<MKL_Complex16*>(g_block_ptr), &nDim);
                    
                    alpha = -alpha;
                    zgeru(&blkSize, &blkSize, reinterpret_cast<MKL_Complex16*>(&alpha), 
                        reinterpret_cast<MKL_Complex16*>(x2_seg.data()), &inc, 
                        reinterpret_cast<MKL_Complex16*>(x1_seg.data()), &inc,
                        reinterpret_cast<MKL_Complex16*>(g_block_ptr), &nDim);
                }
            } else { // hsScheme == 1
                int idxj1, idxk1, idxj2, idxk2;
                config->aux2MajoranaIdx(idxAux, 0, bondType, idxj1, idxk1);
                config->aux2MajoranaIdx(idxAux, 1, bondType, idxj2, idxk2);
                for (int iaux = 0; iaux < 2; iaux++) {
                    int idx1, idx2;
                    DataType tmp_update;
                    if (iaux == 0) {
                        idx1 = idxj1; idx2 = idxj2;
                        tmp_update = (1.0 - ((1.0i) * (config->thlV) * double(auxCur) * g(offset + idx1, offset + idx2)));
                    } else {
                        idx1 = idxk1; idx2 = idxk2;
                        tmp_update = (1.0 + ((1.0i) * (config->thlV) * double(auxCur) * g(offset + idx1, offset + idx2)));
                    }
                    int idx1_g = offset + idx1;
                    int idx2_g = offset + idx2;

                    (*s[replica])(idxAux) = -auxCur;
                    B_replica[replica](idx1, idx2) = -B_replica[replica](idx1, idx2);
                    B_replica[replica](idx2, idx1) = -B_replica[replica](idx2, idx1);

                    DataType alpha = (+1.0i) * double(auxCur) * (config->thlV) / tmp_update;
                    if (iaux == 1) alpha = -alpha;

                    cVecType x1_seg = -g.col(idx1_g).segment(blkOffset, blkSize);
                    cVecType x2_seg = -g.col(idx2_g).segment(blkOffset, blkSize);
                    
                    int idx1_local = idx1_g - blkOffset;
                    int idx2_local = idx2_g - blkOffset;
                    
                    x1_seg(idx1_local) += 2.0; 
                    x2_seg(idx2_local) += 2.0; 
                    
                    zgeru(&blkSize, &blkSize, reinterpret_cast<MKL_Complex16*>(&alpha), 
                        reinterpret_cast<MKL_Complex16*>(x1_seg.data()), &inc, 
                        reinterpret_cast<MKL_Complex16*>(x2_seg.data()), &inc,
                        reinterpret_cast<MKL_Complex16*>(g_block_ptr), &nDim);
                    
                    alpha = -alpha;
                    zgeru(&blkSize, &blkSize, reinterpret_cast<MKL_Complex16*>(&alpha), 
                        reinterpret_cast<MKL_Complex16*>(x2_seg.data()), &inc, 
                        reinterpret_cast<MKL_Complex16*>(x1_seg.data()), &inc,
                        reinterpret_cast<MKL_Complex16*>(g_block_ptr), &nDim);
                }
            }
        }
    }

    void singleFlip(MatType &g, int idxAux, bool &flag, DataType &signCur) {
        // Loop over replica pairs
        for (int p = 0; p < nReplicas / 2; p++) {
            int rA = 2 * p;
            int rB = 2 * p + 1;
            
            // Generate ONE random number for the pair
            double rand = rd->rdUniform01();
            
            DataType r_A = getRatio(g, rA, idxAux);
            DataType r_B = getRatio(g, rB, idxAux);

            // --- Combined Acceptance ---
            DataType r_total = r_A * r_B;
            bool accept = rand < std::abs(r_total);
            
            if (accept) {
                flag = true;
                signCur *= (r_total / std::abs(r_total));
                
                if (assumeBlockDiagonal) {
                    updateBlockPair(g, rA, rB, idxAux);
                } else {
                    updateReplica(g, rA, idxAux);
                    updateReplica(g, rB, idxAux);
                }
            }
        }
    };

    DataType getRatioSingleMajorana(MatType &g, int replica, int idxAux) {
        int offset = replica * nDimSingle;
        int auxCur = (*s[replica])(idxAux);
        int idx1, idx2;
        DataType tmp;
        
        config->aux2MajoranaIdx(idxAux, 0, bondType, idx1, idx2);
        int idx1_g = offset + idx1;
        int idx2_g = offset + idx2;
        
        tmp = (1.0 - ((1.0i) * (config->thlV) * double(auxCur) * g(idx1_g, idx2_g)));
        return tmp * tmp * etaM;
    }

    void updateReplicaSingleMajorana(MatType &g, int replica, int idxAux) {
        int offset = replica * nDimSingle;
        int auxCur = (*s[replica])(idxAux);
        int idx1, idx2;
        DataType tmp;
        const int inc = 1;
        DataType alpha;

        config->aux2MajoranaIdx(idxAux, 0, bondType, idx1, idx2);
        int idx1_g = offset + idx1;
        int idx2_g = offset + idx2;
        
        tmp = (1.0 - ((1.0i) * (config->thlV) * double(auxCur) * g(idx1_g, idx2_g)));

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
        
        zgeru(&nDim, &nDim, reinterpret_cast<MKL_Complex16*>(&alpha), 
              reinterpret_cast<MKL_Complex16*>(x1.data()), &inc, 
              reinterpret_cast<MKL_Complex16*>(x2.data()), &inc,
              reinterpret_cast<MKL_Complex16*>(g.data()), &nDim);
        alpha = -alpha;
        zgeru(&nDim, &nDim, reinterpret_cast<MKL_Complex16*>(&alpha), 
              reinterpret_cast<MKL_Complex16*>(x2.data()), &inc, 
              reinterpret_cast<MKL_Complex16*>(x1.data()), &inc,
              reinterpret_cast<MKL_Complex16*>(g.data()), &nDim);
    }

    void updateBlockPairSingleMajorana(MatType &g, int rA, int rB, int idxAux) {
        int p = rA / 2;
        int blkOffset = p * 2 * nDimSingle;
        int blkSize = 2 * nDimSingle;
        DataType* g_block_ptr = g.data() + blkOffset + blkOffset * nDim;
        const int inc = 1;

        for (int replica : {rA, rB}) {
            int offset = replica * nDimSingle;
            int auxCur = (*s[replica])(idxAux);
            int idx1, idx2;
            DataType tmp;
            DataType alpha;

            config->aux2MajoranaIdx(idxAux, 0, bondType, idx1, idx2);
            int idx1_g = offset + idx1;
            int idx2_g = offset + idx2;
            
            tmp = (1.0 - ((1.0i) * (config->thlV) * double(auxCur) * g(idx1_g, idx2_g)));

            // update aux field and B matrix
            (*s[replica])(idxAux) = -auxCur;
            B_replica[replica](idx1, idx2) = -B_replica[replica](idx1, idx2);
            B_replica[replica](idx2, idx1) = -B_replica[replica](idx2, idx1);

            // update Green's function (block optimized)
            alpha = (+1.0i) * double(auxCur) * (config->thlV) / tmp;

            cVecType x1_seg = -g.col(idx1_g).segment(blkOffset, blkSize);
            cVecType x2_seg = -g.col(idx2_g).segment(blkOffset, blkSize);
            
            int idx1_local = idx1_g - blkOffset;
            int idx2_local = idx2_g - blkOffset;
            
            x1_seg(idx1_local) += 2.0; 
            x2_seg(idx2_local) += 2.0; 
            
            zgeru(&blkSize, &blkSize, reinterpret_cast<MKL_Complex16*>(&alpha), 
                reinterpret_cast<MKL_Complex16*>(x1_seg.data()), &inc, 
                reinterpret_cast<MKL_Complex16*>(x2_seg.data()), &inc,
                reinterpret_cast<MKL_Complex16*>(g_block_ptr), &nDim);
            
            alpha = -alpha;
            zgeru(&blkSize, &blkSize, reinterpret_cast<MKL_Complex16*>(&alpha), 
                reinterpret_cast<MKL_Complex16*>(x2_seg.data()), &inc, 
                reinterpret_cast<MKL_Complex16*>(x1_seg.data()), &inc,
                reinterpret_cast<MKL_Complex16*>(g_block_ptr), &nDim);
        }
    }

    void singleFlipSingleMajorana(MatType &g, int idxAux, bool &flag, DataType &signCur) {
        for (int p = 0; p < nReplicas / 2; p++) {
            int rA = 2 * p;
            int rB = 2 * p + 1;
            
            double rand = rd->rdUniform01();
            
            DataType r_A = getRatioSingleMajorana(g, rA, idxAux);
            DataType r_B = getRatioSingleMajorana(g, rB, idxAux);

            DataType r_total = r_A * r_B;
            bool accept = rand < std::abs(r_total);
            
            if (accept) {
                flag = true;
                signCur *= (r_total / std::abs(r_total));
                
                if (assumeBlockDiagonal) {
                    updateBlockPairSingleMajorana(g, rA, rB, idxAux);
                } else {
                    updateReplicaSingleMajorana(g, rA, idxAux);
                    updateReplicaSingleMajorana(g, rB, idxAux);
                }
            }
        }
    };

    DataType update(MatType &g) override {
        bool flag = false;
        DataType signCur = 1.0;
        
        // Loop over all replicas
        // Loop over all aux indices
        for (int i = 0; i < s[0]->size(); i++) {
            if (singleMaj) {
                singleFlipSingleMajorana(g, i, flag, signCur);
            } else {
                singleFlip(g, i, flag, signCur);
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