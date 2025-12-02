#ifndef CHAIN1D_TV_H
#define CHAIN1D_TV_H

#include "types.h"
#include "operator.h"
#include "skewMatUtils.h"
#include "spinless_tV.h"

// A simple 1D t-V chain with:
// - Lx sites
// - 2 Majorana modes per site (spinless fermions)
// - Nearest-neighbor hopping (t)
// - Nearest-neighbor interaction (V)
// - Boundary conditions: 0 = PBC, 1 = OBC
class SpinlessTvChain1dUtils: public SpinlessTvUtils
{
public:
    int boundaryType;  // 0 = PBC, 1 = OBC
    int nBonds;        // Number of bonds

    SpinlessTvChain1dUtils(int _Lx, double _dt, double _V, int _l,
                           int _boundaryType = 0, int _hsScheme = 0)
        : SpinlessTvUtils(_Lx, 1, _dt, _V, _l, _Lx * 2, false, _hsScheme),
          boundaryType(_boundaryType) {

        // Number of bonds
        nBonds = (boundaryType == 0) ? Lx : (Lx - 1);
    }

    // Map site index and Majorana species to flat index
    // idx = site * 2 + imaj
    inline int siteToIdx(int site, int imaj) const {
        return site * 2 + imaj;
    }

    // Map auxiliary field index to Majorana indices
    // For 1D chain, auxiliary field lives on bonds
    inline void aux2MajoranaIdx(int idAux, int imaj, int bType, int& idx1, int& idx2) const override {
        // idAux is the bond index
        // Bond idAux connects site idAux to site (idAux+1)
        idx1 = siteToIdx(idAux, imaj);
        idx2 = siteToIdx((idAux + 1) % Lx, imaj);
    }

    // Generate kinetic Hamiltonian matrix
    // H_kin = -t * sum_i (c_i^dag c_{i+1} + h.c.)
    // In Majorana: = i*t * sum_i (gamma1_i * gamma2_{i+1} - gamma2_i * gamma1_{i+1})
    inline void KineticGenerator(MatType &H, DataType t) const {
        H.setZero();
        DataType tmp = (1.0i) * t;

        int nSitesToSum = (boundaryType == 0) ? Lx : (Lx - 1);

        for (int i = 0; i < nSitesToSum; i++) {
            int j = (i + 1) % Lx;

            // gamma1_i * gamma2_j term
            int idx1 = siteToIdx(i, 0);
            int idx2 = siteToIdx(j, 1);
            H(idx1, idx2) = tmp;
            H(idx2, idx1) = -tmp;

            // -gamma2_i * gamma1_j term
            idx1 = siteToIdx(i, 1);
            idx2 = siteToIdx(j, 0);
            H(idx1, idx2) = -tmp;
            H(idx2, idx1) = tmp;
        }
    }

    // Compute energy from Green's function
    // E = <H_kin> + <H_int>
    inline DataType energyFromGreensFunc(const MatType &g) const {
        DataType r = 0.0;

        // Kinetic energy: <H_kin> = sum_<ij> t_ij * i/2 * G_ji
        DataType tmp = (0.5i);
        int nSitesToSum = (boundaryType == 0) ? Lx : (Lx - 1);

        for (int i = 0; i < nSitesToSum; i++) {
            int j = (i + 1) % Lx;

            // gamma1_i * gamma2_j term
            int idx1 = siteToIdx(i, 0);
            int idx2 = siteToIdx(j, 1);
            r += tmp * g(idx2, idx1);

            // -gamma2_i * gamma1_j term
            idx1 = siteToIdx(i, 1);
            idx2 = siteToIdx(j, 0);
            r -= tmp * g(idx2, idx1);
        }

        // Interaction energy: <H_int> = V/4 * sum_<ij> (G_i1j1 * G_i2j2 + G_i1j2 * G_j1i2 - G_i1i2 * G_j1j2)
        tmp = 0.25 * V;
        for (int bond = 0; bond < nBonds; bond++) {
            int i = bond;
            int j = (bond + 1) % Lx;

            int idxi1 = siteToIdx(i, 0);
            int idxi2 = siteToIdx(i, 1);
            int idxj1 = siteToIdx(j, 0);
            int idxj2 = siteToIdx(j, 1);

            r += tmp * g(idxi1, idxj1) * g(idxi2, idxj2);
            r += tmp * g(idxi1, idxj2) * g(idxj1, idxi2);
            r -= tmp * g(idxi1, idxi2) * g(idxj1, idxj2);
        }

        return r;
    }

    // Particle number: N = sum_i n_i = sum_i (1 - i*gamma1_i*gamma2_i)/2
    inline DataType particleNumber(const MatType &g) const {
        DataType r = 0.0;
        for (int i = 0; i < Lx; i++) {
            int idx1 = siteToIdx(i, 0);
            int idx2 = siteToIdx(i, 1);
            // n_i = (1 - i*gamma1*gamma2)/2
            // <n_i> = 1/2 - i/2 * <gamma1*gamma2> = 1/2 + 1/2 * G_21
            r += 0.5 * (1.0 + g(idx2, idx1));
        }
        return r;
    }

    // Density-density structure factor S_rho(q) = 1/L * sum_{i,j} e^{iq(i-j)} <n_i n_j>
    inline DataType densityStructureFactor(const MatType &g, double q) const {
        DataType result = 0.0;

        for (int i = 0; i < Lx; i++) {
            for (int j = 0; j < Lx; j++) {
                DataType phase = std::exp(1.0i * q * static_cast<double>(i - j));

                int idx_i1 = siteToIdx(i, 0);
                int idx_i2 = siteToIdx(i, 1);
                int idx_j1 = siteToIdx(j, 0);
                int idx_j2 = siteToIdx(j, 1);

                // <n_i n_j> = 1/4 * (1 + G_i21 + G_j21 + G_i21*G_j21 - G_i12*G_j21 - G_i21*G_j12 + G_i11*G_j22 + G_i12*G_j12)
                DataType ni_nj = 0.25 * (
                    1.0 + g(idx_i2, idx_i1) + g(idx_j2, idx_j1)
                    + g(idx_i2, idx_i1) * g(idx_j2, idx_j1)
                    - g(idx_i1, idx_i2) * g(idx_j2, idx_j1)
                    - g(idx_i2, idx_i1) * g(idx_j1, idx_j2)
                    + g(idx_i1, idx_j1) * g(idx_i2, idx_j2)
                    + g(idx_i1, idx_i2) * g(idx_j1, idx_j2)
                );

                result += phase * ni_nj;
            }
        }

        return result / static_cast<double>(Lx);
    }
};


class Chain1d_tV: public Spinless_tV
{
public:
    const SpinlessTvChain1dUtils* modelConfig;
    double dt;
    int l;
    int nSites, nBonds;
    rdGenerator* rd;

    Chain1d_tV(SpinlessTvChain1dUtils* _config, rdGenerator* _rd) {
        modelConfig = _config;
        dt = modelConfig->dt;
        l = modelConfig->l;
        nSites = modelConfig->Lx;
        nBonds = modelConfig->nBonds;
        rd = _rd;

        nDim = nSites * 2;  // 2 Majorana modes per site

        // Generate kinetic operators
        MatType Ht(nDim, nDim);
        modelConfig->KineticGenerator(Ht, 1.0);
        MatType expK = expm(Ht, -dt);
        MatType expKhalf = expm(Ht, -dt / 2.0);

        // Build operator array: expK/2, (expK, V_interaction) x l, expK/2
        op_array = std::vector<Operator*>(2 * l + 1);
        iVecType* s;

        // First half kinetic step
        op_array[0] = new DenseOperator(expKhalf, 1.0);

        // Alternating kinetic and interaction operators
        for (int i = 0; i < l; i++) {
            // Full kinetic step
            if (i > 0) {
                op_array[2 * i] = new DenseOperator(expK, 1.0);
            }

            // Interaction operator with auxiliary fields on all bonds
            s = new iVecType(nBonds);
            for (int k = 0; k < nBonds; k++) {
                (*s)(k) = rd->rdZ2();
            }
            op_array[2 * i + 1] = new SpinlessVOperator(modelConfig, s, 0, rd);
        }

        // Final half kinetic step
        op_array[2 * l] = new DenseOperator(expKhalf, 1.0);
    }
};

#endif