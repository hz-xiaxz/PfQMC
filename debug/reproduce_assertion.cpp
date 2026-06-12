#include <iostream>
#include <vector>
#include <complex>
#include "inc/honeycomb.h"
#include "inc/spinless_tV.h"
#include "inc/types.h"
#include "inc/mixing_operator.h"
#include "inc/pfqmc.h"

int main() {
    int Lx = 2;
    int Ly = 2;
    int LTau = 10;
    double dt = 0.1;
    double V = 0.0;
    int nReplicas = 4;

    SpinlessTvHoneycombUtils config(Lx, Ly, dt, V, LTau);
    rdGenerator rd(42);

    int nSites = Lx * Ly * 2;
    int nDimSingle = nSites * 2; // 16
    int nDim = nReplicas * nDimSingle; // 64

    std::cout << "nDimSingle: " << nDimSingle << std::endl;
    std::cout << "nDim: " << nDim << std::endl;

    // Create 4-replica operators
    std::vector<Operator*> op_array(4 * LTau + 1);
    
    // Helper to create 4-replica DenseOperator
    auto createDenseOp = [&](const MatType& m) {
        MatType m_full = MatType::Zero(4*nDimSingle, 4*nDimSingle);
        for(int r=0; r<4; r++) m_full.block(r*nDimSingle, r*nDimSingle, nDimSingle, nDimSingle) = m;
        return new DenseOperator(m_full, 1.0);
    };

    MatType Ht(nDimSingle, nDimSingle);
    config.KineticGenerator(Ht, 1.0);
    MatType expK = expm(Ht, -dt);
    MatType expKhalf = expm(Ht, -dt / 2.0);

    for (int i=0; i<LTau; i++) {
        if (i == 0) {
            op_array[0] = createDenseOp(expKhalf);
        } else {
            op_array[4*i] = createDenseOp(expK);
        }

        for (int j=0; j<3; j++) {
            std::vector<iVecType*> s_vec(4);
            for(int r=0; r<4; r++) {
                s_vec[r] = new iVecType(Lx*Ly); // nUnitcell
                for (int k=0; k<Lx*Ly; k++) (*s_vec[r])(k) = rd.rdZ2();
            }
            op_array[4*i + j + 1] = new SpinlessVOperator(&config, s_vec, j, &rd, 4);
        }
    }
    op_array[4*LTau] = createDenseOp(expKhalf);

    // Create a temporary walker struct to pass to PfQMC
    Spinless_tV walker;
    walker.op_array = op_array;
    walker.nDim = 4 * nDimSingle; // Total dimension
    walker.nDimSingle = nDimSingle; // Single replica dimension

    MixingOperator* mixingOp = new TrivialSwapOperator(nReplicas, nDimSingle, &rd);
    
    std::cout << "Creating PfQMC..." << std::endl;
    PfQMC pfqmc(&walker, 10, nReplicas, mixingOp);
    std::cout << "PfQMC created." << std::endl;

    std::cout << "Running rightSweep..." << std::endl;
    pfqmc.rightSweep();
    std::cout << "rightSweep finished." << std::endl;

    // Check g matrix structure
    std::cout << "Checking g matrix structure..." << std::endl;
    MatType& g = pfqmc.g;
    int N = nDimSingle;
    
    for(int r1=0; r1<nReplicas; r1++) {
        for(int r2=0; r2<nReplicas; r2++) {
            double norm = g.block(r1*N, r2*N, N, N).norm();
            std::cout << "Block (" << r1 << "," << r2 << ") norm: " << norm << std::endl;
        }
    }

    return 0;
}
