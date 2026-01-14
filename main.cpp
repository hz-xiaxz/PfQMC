#include "mpi.h"
#include <omp.h>
#include <string.h>
#include <fstream>

#include "inc/honeycomb.h"
#include "inc/pfqmc.h"
#include "inc/square.h"
#include "inc/singleMajoranaHoneycomb.h"
#include "inc/kitaevChain.h"
#include "inc/mixing_operator.h"

void fixSign(DataType& sign, DataType signRaw, double threshold) {
    // unreliable signRaw
    if (std::abs( 1.0 - std::abs(signRaw.real()) ) > threshold)
    {
        // reliable sign
        if (std::abs( 1.0 - std::abs( sign.real() ) ) < threshold) {
            if (sign.real() > 0.0) {
                sign = 1.0;
            } else {
                sign = -1.0;
            }
        } else {
            std::cout << "=== error in sign === " << " sign = " << sign << " ≠ " << signRaw << "==== \n"; 
        }
    } else {
        // reliable signRaw
        if (signRaw.real() > 0.0) {
            sign = 1.0;
        } else {
            sign = -1.0;
        }
    }
}


int main_honeycomb() {
    double start_time = omp_get_wtime();
    mkl_set_num_threads(8);

    int Lx = 4;
    int Ly = 4;
    int LTau = 100;
    double dt = 0.1;
    double V = 0.7;
    int stabilizationTime = 10;
    int thermalLength = 200;
    int evaluationLength = 1000;

    // int nDim = Lx * Ly * 4;
    SpinlessTvHoneycombUtils config(Lx, Ly, dt, V, LTau);
    rdGenerator rd(42);
    Honeycomb_tV walker(&config, &rd);
    PfQMC pfqmc(&walker, stabilizationTime);
    for (int i = 0; i < thermalLength; i++) {
        pfqmc.rightSweep();
        pfqmc.leftSweep();
        // std::cout << i << std::endl;
    }

    DataType energy = 0.0;
    DataType sign = 1.0;
    for (int i = 0; i < evaluationLength; i++) {
        pfqmc.rightSweep();
        pfqmc.leftSweep();
        // check if the sign problem free condition is met
        if (std::abs(1.0 - pfqmc.sign) > 1e-2) {
            std::cout << "\n=== error in sign at round = " << i << " sign = " << pfqmc.sign << " ≠ " << 1.0 << "==== \n";
        }
        // sign = pfqmc.getSignRaw();
        energy += sign * config.energyFromGreensFunc(pfqmc.g);

        if ((i % 10) == 0)
            std::cout << "iter = " << i << " energy=" << energy / double(i + 1)
                      << "\n";
    }
    // std::cout << pfqmc.g << "\n";

    double time = omp_get_wtime() - start_time;
    std::cout << "total time " << time << std::endl;
    return 0;
}

int main_honeycomb_param(int Lx, int Ly, int LTau, double dt, double V, int nthreads, int nseed, int evaluationLength, char* filename) {
    double start_time = omp_get_wtime();
    mkl_set_num_threads(nthreads);

    std::fstream fout(filename, std::fstream::out);

    int stabilizationTime = 10;
    int thermalLength = 200;

    fout << "=== Standard Honeycomb Model ===\n";
    fout << "Lx = " << Lx << " Ly = " << Ly << " LTau = " << LTau << " dt = " << dt << " V = " << V << " seed = " << nseed << " nthreads = " << nthreads << " evaluationLength = " << evaluationLength << std::endl;

    SpinlessTvHoneycombUtils config(Lx, Ly, dt, V, LTau);
    rdGenerator rd(nseed);
    Honeycomb_tV walker(&config, &rd);
    PfQMC pfqmc(&walker, stabilizationTime);

    for (int i = 0; i < thermalLength; i++) {
        fout << i << " " << std::flush;
        pfqmc.rightSweep();
        pfqmc.leftSweep();
    }
    fout << std::endl;

    DataType energy = 0.0;
    DataType sign, signRaw;
    DataType signTot = 0.0;
    DataType energyTot = 0.0;

    for (int i = 0; i < evaluationLength; i++) {
        pfqmc.rightSweep();
        pfqmc.leftSweep();
        sign = pfqmc.sign;

        if (i % 20 == 0) {
            signRaw = pfqmc.getSignRaw();
            double threshold = 1e-2;
            if (std::abs(sign - signRaw) > threshold) {
                fixSign(sign, signRaw, threshold);
                pfqmc.sign = sign;
            }
        }

        energy = config.energyFromGreensFunc(pfqmc.g);
        energyTot += sign * energy;
        signTot += sign;

        fout << "iter = " << i << " energy = " << energy << " sign = " << sign << std::endl;

        if ((i % 10) == 0)
            std::cout << "iter = " << i << " energy=" << energyTot / signTot << " sign=" << signTot / double(i + 1) << "\n";
    }

    std::cout << "Final: AveEnergy = " << energyTot / signTot << " AveSign = " << signTot / double(evaluationLength) << std::endl;

    double time = omp_get_wtime() - start_time;
    fout << "=== End of Standard Honeycomb Model ===\n";
    fout << "total time " << time << std::endl;
    return 0;
}

int main_honeycombSingleMajorana(int Lx, int Ly, int LTau, double dt, double V, int nthreads, int nseed, int evaluationLength, char* filename) {
    double start_time = omp_get_wtime();
    mkl_set_num_threads(nthreads);

    std::fstream fout(filename, std::fstream::out);

    int stabilizationTime = 10;
    int thermalLength = 200;
    // int evaluationLength = 1000;
    fout << "=== Honeycomb Single Majorana Model ===\n";
    fout << "Lx = " << Lx << " Ly = " << Ly << " LTau = " << LTau << " dt = " << dt << " V = " << V << " seed = " << nseed << " nthreads = " << nthreads << " evaluationLength = " << evaluationLength << std::endl;
    SpinlessTvHoneycombSingleMajoranaUtils config(Lx, Ly, dt, V, LTau);
    rdGenerator rd(nseed);
    HoneycombSingleMajorana_tV walker(&config, &rd);
    PfQMC pfqmc(&walker, stabilizationTime);
    for (int i = 0; i < thermalLength; i++) {
	    fout << i << " " << std::flush;
        pfqmc.rightSweep();
        pfqmc.leftSweep();
    }
    fout << std::endl;

    // DataType energy = 0.0;
    DataType structureFactorCDW = 0.0;
    DataType structureFactorCDWq = 0.0;
    DataType sign, signRaw;

    DataType srSignTot = 0.0;
    DataType CDWTot = 0.0;
    DataType CDWqTot = 0.0;

    for (int i = 0; i < evaluationLength; i++) {
        pfqmc.rightSweep();
        pfqmc.leftSweep();
        sign = pfqmc.sign;

        if (i % 20 == 0) {
            signRaw = pfqmc.getSignRaw();
            double threshold = 1e-2;
            if (std::abs(sign - signRaw) > threshold) {
                fixSign(sign, signRaw, threshold);
                pfqmc.sign = sign;
            }
	        // pfqmc.sign = signRaw;
        }

        // srSignTot += sign;
        // energy = config.energyFromGreensFunc(pfqmc.g);
        structureFactorCDW = config.structureFactorCDW(pfqmc.g);
        structureFactorCDWq = config.structureFactorCDWoffset(pfqmc.g);

        fout << "i = " << i << " CDW = " << structureFactorCDW << " MRsign = " << sign << " CDWq = " << structureFactorCDWq << std::endl;

        srSignTot += sign;
        CDWTot += structureFactorCDW;
        CDWqTot += structureFactorCDWq;

        if (i == evaluationLength - 1) {
            std::cout << "Average MRsign = " << srSignTot / double(evaluationLength) << " AveCDW = " << CDWTot / double(evaluationLength) << " AveCDWq = " << CDWqTot / double(evaluationLength) << std::endl;
        }

        
    }
    // std::cout << pfqmc.g << "\n";

    double time = omp_get_wtime() - start_time;
    fout << "=== End of Honeycomb Single Majorana Model ===\n";
    fout << "total time " << time << std::endl;
    return 0;
}

int main_honeycomb_4replica(int Lx, int Ly, int LTau, double dt, double V, int nthreads, int nseed, int evaluationLength, char* filename, bool useSwap = false) {
    double start_time = omp_get_wtime();
    mkl_set_num_threads(nthreads);

    std::fstream fout(filename, std::fstream::out);

    int stabilizationTime = 10;
    int thermalLength = 200;
    int nReplicas = 4;
    
    fout << "=== Honeycomb 4-Replica Model ===\n";
    fout << "Lx = " << Lx << " Ly = " << Ly << " LTau = " << LTau << " dt = " << dt << " V = " << V << " seed = " << nseed << " nthreads = " << nthreads << " evaluationLength = " << evaluationLength << std::endl;
    
    SpinlessTvHoneycombUtils config(Lx, Ly, dt, V, LTau);
    rdGenerator rd(nseed);
    
    // Create mixing operator (identity for now, or specific mixing)
    // For Renyi entropy S2, we need SWAP operator between replicas
    // But for now let's just use Identity to verify stability
    int nDimSingle = config.nDim;
    MixingOperator* mixingOp;
    if (useSwap) {
        mixingOp = new TrivialSwapOperator(nReplicas, nDimSingle, &rd);
    } else {
        mixingOp = new MixingOperator(nReplicas, nDimSingle, &rd);
    }
    
    // Create walker with 4 replicas
    // We need to modify Honeycomb_tV to support replicas or create a new one
    // But Honeycomb_tV uses SpinlessVOperator which now supports replicas.
    // However, Honeycomb_tV constructor initializes op_array with single replica operators.
    // We need a way to initialize it with 4-replica operators.
    
    // For now, let's assume we can use a modified walker or just manually build op_array here?
    // Or better, modify Honeycomb_tV to accept nReplicas.
    
    // Let's modify Honeycomb_tV in inc/honeycomb.h first to support nReplicas.
    // But since I cannot modify header in this step, I will use a local construction here if possible
    // or just note that Honeycomb_tV needs update.
    
    // Actually, Honeycomb_tV is simple enough to inline here or just use if updated.
    // Let's assume Honeycomb_tV is NOT updated yet.
    // I will construct the walker manually here for 4 replicas.
    
    int nSites = Lx * Ly * 2;
    int nDim = nSites * 2; // single replica dimension
    
    MatType Ht(nDim, nDim);
    config.KineticGenerator(Ht, 1.0);
    MatType expK = expm(Ht, -dt);
    MatType expKhalf = expm(Ht, -dt / 2.0);
    
    // Create 4-replica operators
    std::vector<Operator*> op_array(4 * LTau + 1);
    iVecType* s;
    
    // Helper to create 4-replica DenseOperator
    auto createDenseOp = [&](const MatType& m) {
        MatType m_full = MatType::Zero(4*nDim, 4*nDim);
        for(int r=0; r<4; r++) m_full.block(r*nDim, r*nDim, nDim, nDim) = m;
        return new DenseOperator(m_full, 1.0);
    };

    for (int i=0; i<LTau; i++) {
        if (i == 0) {
            op_array[0] = createDenseOp(expKhalf);
        } else {
            op_array[4*i] = createDenseOp(expK);
        }

        for (int j=0; j<3; j++) {
            std::vector<iVecType*> s_vec(4);
            
            // Generate for replica 0
            s_vec[0] = new iVecType(Lx*Ly);
            for (int k=0; k<Lx*Ly; k++) (*s_vec[0])(k) = rd.rdZ2();
            
            // Copy to other replicas to ensure identical initialization
            for(int r=1; r<4; r++) {
                s_vec[r] = new iVecType(Lx*Ly);
                *s_vec[r] = *s_vec[0];
            }
            
            op_array[4*i + j + 1] = new SpinlessVOperator(&config, s_vec, j, &rd, 4);
        }
    }
    op_array[4*LTau] = createDenseOp(expKhalf);
    
    // Create a temporary walker struct to pass to PfQMC
    Spinless_tV walker;
    walker.op_array = op_array;
    walker.nDim = 4 * nDim; // Total dimension
    walker.nDimSingle = nDim; // Single replica dimension
    
    PfQMC pfqmc(&walker, stabilizationTime, nReplicas, mixingOp);
    
    for (int i = 0; i < thermalLength; i++) {
        fout << i << " " << std::flush;
        pfqmc.rightSweep();
        pfqmc.leftSweep();
    }
    fout << std::endl;

    DataType energy = 0.0;
    DataType sign = 1.0;
    
    for (int i = 0; i < evaluationLength; i++) {
        pfqmc.rightSweep();
        pfqmc.leftSweep();

        // Simple energy measurement (averaged over replicas)
        // We need to extract blocks from pfqmc.g
        DataType currentEnergy = 0.0;
        for(int r=0; r<nReplicas; r++) {
            currentEnergy += config.energyFromGreensFuncReplica(pfqmc.g, r, nDim);
        }
        // currentEnergy /= nReplicas; // User requested total energy

        energy += currentEnergy;

        if ((i % 10) == 0)
            std::cout << "iter = " << i << " energy=" << energy / double(i + 1) << "\n";

        fout << "iter = " << i << " energy = " << currentEnergy << std::endl;
    }

    std::cout << "Final: AveEnergy = " << energy / double(evaluationLength) << std::endl;

    double time = omp_get_wtime() - start_time;
    fout << "total time " << time << std::endl;
    return 0;
}

int main_square() {
    double start_time = omp_get_wtime();
    mkl_set_num_threads(8);

    int Lx = 4;
    int Ly = 4;
    int LTau = 200;
    double dt = 0.1;
    double V = 0.7;
    double delta = 0.5;
    int stabilizationTime = 10;
    int thermalLength = 200;
    int evaluationLength = 1000;

    // int nDim = Lx * Ly * 4;
    SpinlessTvSquareUtils config(Lx, Ly, dt, V, LTau, delta);
    rdGenerator rd(42);
    Square_tV walker(&config, &rd);
    PfQMC pfqmc(&walker, stabilizationTime);
    // std::cout << "here!" << std::endl;
    for (int i = 0; i < thermalLength; i++) {
        pfqmc.rightSweep();
        pfqmc.leftSweep();
    }

    DataType energy = 0.0;
    DataType sign = 1.0;
    DataType signFast;
    DataType signTot = 0.0;
    DataType rawSignTot = 0.0;
    DataType rawSign = 1.0;
    // pfqmc.sign = 1.0;
    for (int i = 0; i < evaluationLength; i++) {
        pfqmc.rightSweep();
        pfqmc.leftSweep();
        sign = pfqmc.sign;

        // rawSign = pfqmc.getSignRaw();
        rawSignTot += rawSign;
        // double deviation = std::abs(sign - rawSign);
        // if (deviation > 1e-2) {
        //     pfqmc.sign = rawSign; // update sign
        //     std::cout << "\n=== error in sign at round = " << i << " raw sign = " << rawSign << " ≠ " << sign << "==== \n"; 
        // }

        energy += sign * config.energyFromGreensFunc(pfqmc.g);
        signTot += sign;
        // std::cout << "sign = " << sign << " raw sign = " << rawSign << "\n";
        if ((i % 10) == 0)
            std::cout << "iter = " << i << " sign ave = " << signTot / double(i + 1) << " raw sign ave= " << rawSignTot / double(i + 1) 
                      << " energy=" << energy / signTot << " curSign = "  <<  sign << "\n";
    }
    // std::cout << pfqmc.g << "\n";

    double time = omp_get_wtime() - start_time;
    std::cout << "total time " << time << std::endl;
    return 0;
}

int main_chain(int Lx, int LTau, double dt, double V, double delta,int nthreads, int nseed, int evaluationLength, char* filename, double mu=0.0, int hsScheme=0, int boundary=1) {
    double start_time = omp_get_wtime();
    mkl_set_num_threads(nthreads);

    std::fstream fout(filename, std::fstream::out);

    int stabilizationTime = 10;
    int thermalLength = 1000;
    // int boundary = 0;
    int aveLength = 100;
    // int evaluationLength = 1000;

    // int nDim = Lx * Ly * 4;
    SpinlessTvChainUtils config(Lx, dt, V, LTau, boundary, delta, mu, hsScheme);
    rdGenerator rd(nseed);
    Chain_tV walker(&config, &rd);
    PfQMC pfqmc(&walker, stabilizationTime);
    // std::cout << "here!\n";
    fout << "=== p-wave Chain Model with chemical potential mu = " << mu << " ===\n";
    fout << "Lx = " << Lx << " LTau = " << LTau << " dt = " << dt << " V = " << V << " seed = " << nseed << " nthreads = " << nthreads
    << " delta = " << delta << " boundary = " << boundary << " evaluationLength = " << evaluationLength << std::endl;
    // std::cout << "here!" << std::endl;
    for (int i = 0; i < thermalLength; i++) {
        fout << i << " " << std::flush;
        // std::cout << "sign = " << pfqmc.sign << "\n";
        pfqmc.rightSweep();
        pfqmc.leftSweep();
    }
    fout << std::endl;

    // DataType energy = 0.0;
    DataType structureFactorCDW, structureFactorCDWq;
    DataType sign, signRaw;
    DataType obsZ2, obsEnergy, obsEdgeCorrelator, obsEdgeCorrelatorZ2;

    DataType SignTot = 0.0;
    DataType z2PlusSignTot = 0.0;
    DataType z2MinusSignTot = 0.0;
    DataType structureFactorCDWTot = 0.0;

    DataType obsEnergyTot = 0.0;
    DataType obsSignTotTrue = 0.0;
    DataType obsEdgeCorrelatorTot = 0.0;
    DataType obsEdgeCorrelatorZ2Tot = 0.0;
    DataType obsZ2Tot = 0.0;
    DataType obsNelectronsTot = 0.0;
    // DataType obsEdgeCorrelatorZ2PlusTot = 0.0;
    // DataType obsEdgeCorrelatorZ2MinusTot = 0.0;
    DataType obsStructureFactorCDWTot = 0.0;
    // DataType obsStructureFactorCDWM4Tot = 0.0;
    DataType obsStructureFactorCDWqTot = 0.0;

    // fourier transform of the CDW structure factor
    // CDW should corresponds to obsCDWFT[\pi]
    // CDWq should corresponds to obsCDWFT[\pi + 2\pi / Lx]
    cVecType obsCDWFTTot = cVecType::Zero(Lx);
    cVecType obsCDWFT;

    // pfqmc.sign = pfqmc.getSignRaw(); // initialize sign
    for (int i = 0; i < evaluationLength; i++) {
        pfqmc.rightSweep();
        pfqmc.leftSweep();
        sign = pfqmc.sign;

        if (i % 20 == 0) {
            signRaw = pfqmc.getSignRaw();
            double threshold = 1e-2;
            if (std::abs(sign - signRaw) > threshold) {
                fixSign(sign, signRaw, threshold);
                pfqmc.sign = sign;
            }
	        // pfqmc.sign = signRaw;
        }
        
        obsEnergy = config.energyFromGreensFunc(pfqmc.g);
        obsEnergyTot += sign * obsEnergy;
        obsSignTotTrue += sign;

        obsEdgeCorrelator = config.EdgeCorrelator(pfqmc.g);
        obsZ2 = config.Z2FermionParity(pfqmc.g);
        obsEdgeCorrelatorZ2 = config.Z2FermionParityEdgeCorrelator(pfqmc.g);
        structureFactorCDW = config.StructureFactorCDW(pfqmc.g);
        // structureFactorCDWM4 = config.StructureFactorCDWM4(pfqmc.g);
        structureFactorCDWq = config.StructureFactorCDWOffset(pfqmc.g);
        config.StructureFactorCDWFFT(pfqmc.g, obsCDWFT);

        // ------------- output each iteration -----------
        fout << "iter = " << i << " sign = " << sign << " z2 = " << obsZ2 << " edgeCorrelator = " << obsEdgeCorrelator << " edgeZ2Correlator = " << obsEdgeCorrelatorZ2 << " CDW = " << structureFactorCDW << " CDWq = " << structureFactorCDWq << std::endl;

        // ------------- end of output each iteration ----


        // fout << "iter = " << i << " sign = " << sign << " z2 = " << obsZ2 << " edgeCorrelator = " << obsEdgeCorrelator << " edgeZ2Correlator = " << obsEdgeCorrelatorZ2 
        // << " CDW = " << structureFactorCDW 
        // << " CDWM4 = " << structureFactorCDWM4 << " CDWq = " << structureFactorCDWq << std::endl;

        // SignTot += sign;
        // obsZ2Tot += sign * obsZ2;
        // // z2PlusSignTot += sign * (1.0 + obsZ2) / 2.0;
        // // z2MinusSignTot += sign * (1.0 - obsZ2) / 2.0;
        // obsEdgeCorrelatorTot += sign * obsEdgeCorrelator;
        // obsEdgeCorrelatorZ2Tot += sign * obsEdgeCorrelatorZ2;
        // // obsEdgeCorrelatorZ2PlusTot += sign * (obsEdgeCorrelator + obsEdgeCorrelatorZ2) / 2.0;
        // // obsEdgeCorrelatorZ2MinusTot += sign * (obsEdgeCorrelator - obsEdgeCorrelatorZ2) / 2.0;
        // obsStructureFactorCDWTot += sign * structureFactorCDW;
        // // obsStructureFactorCDWM4Tot += sign * structureFactorCDWM4;
        // obsStructureFactorCDWqTot += sign * structureFactorCDWq;
        // obsCDWFTTot += sign * obsCDWFT;

        // std::cout << "iter = " << i << " cdw = " << structureFactorCDW << " cdw also = " << obsCDWFT[4] << std::endl;

        // if ((i+1) % aveLength == 0) {
        //     double aveLengthD = double(aveLength);
        //     fout << "iter = " << i << " sign = " << SignTot / aveLengthD << " z2 = " << obsZ2Tot / aveLengthD << " edgeCorrelator = " << obsEdgeCorrelatorTot / aveLengthD << " edgeZ2Correlator = " << obsEdgeCorrelatorZ2Tot / aveLengthD << " CDW = " << obsStructureFactorCDWTot / aveLengthD << " CDWq = " << obsStructureFactorCDWqTot / aveLengthD << std::endl;

        //     // fout << "obsCDWFT = ";
        //     for (int j = 0; j < Lx; j++) {
        //         fout << obsCDWFTTot[j].real() / aveLengthD << " ";
        //     }
        //     fout << std::endl;
        //     SignTot = 0.0;
        //     obsZ2Tot = 0.0;
        //     z2PlusSignTot = 0.0;
        //     z2MinusSignTot = 0.0;
        //     obsEdgeCorrelatorTot = 0.0;
        //     obsEdgeCorrelatorZ2Tot = 0.0;
        //     obsStructureFactorCDWTot = 0.0;
        //     // obsStructureFactorCDWM4Tot = 0.0;
        //     obsStructureFactorCDWqTot = 0.0;
        //     obsCDWFTTot = cVecType::Zero(Lx);
        // }

        if (i == evaluationLength - 1) {
        // if ( i % 20 == 0) {
            std::cout << filename << " finished" << std::endl;
            std::cout << "AveEnergy = " << obsEnergyTot / obsSignTotTrue << " sign = " << obsSignTotTrue / double(evaluationLength) << "\n"; 
            // << "AveSign = " << SignTot / double(i + 1) 
            // << " AveZ2 = " << obsZ2Tot / SignTot 
            // << " AveZ2PlusSign = " << z2PlusSignTot / double(i + 1) 
            // << " AveZ2MinusSign = " << z2MinusSignTot / double(i + 1) 
            // << " AveEdgeCorrelator = " << obsEdgeCorrelatorTot / SignTot 
            // << " AveEdgeCorrelatorZ2Plus = " << obsEdgeCorrelatorZ2PlusTot / z2PlusSignTot 
            // << " AveEdgeCorrelatorZ2Minus = " << obsEdgeCorrelatorZ2MinusTot / z2MinusSignTot
            // << " AveStructureFactorCDW = " << obsStructureFactorCDWTot / SignTot
            // << " AveStructureFactorCDWM4 = " << obsStructureFactorCDWM4Tot / SignTot
            // << " AveStructureFactorCDWq = " << obsStructureFactorCDWqTot / SignTot
            //           << "\n";
        }
    }
    // std::cout << pfqmc.g << "\n";

    double time = omp_get_wtime() - start_time;
    fout << "=== End of p-wave Chain Model ===\n";
    fout << "total time " << time << std::endl;
    return 0;
}


int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " --square | --honeycomb | --MRhoneycomb [PARAMS]\n";
        return 0;
    }

    MPI_Init(&argc,&argv);
    int ok = 0;

    if (std::strcmp(argv[1], "--square") == 0) {
        return main_square();
    } else if (std::strcmp(argv[1], "--honeycomb") == 0) {
        return main_honeycomb();
    } else if (std::strcmp(argv[1], "--honeycomb2") == 0) {
        // Parameterized standard honeycomb model
        // Usage: --honeycomb2 Lx Ly LTau dt V nthreads nseed evaluationLength filepath
        int Lx, Ly, LTau, nthreads, nseed, evaluationLength;
        char* filepath;
        double dt, V;
        char filename[100];

        Lx = std::stoi(argv[2]);
        Ly = std::stoi(argv[3]);
        LTau = std::stoi(argv[4]);
        dt = std::stod(argv[5]);
        V = std::stod(argv[6]);
        nthreads = std::stoi(argv[7]);
        nseed = std::stoi(argv[8]);
        evaluationLength = std::stoi(argv[9]);
        filepath = argv[10];

        int numprocs, myid;
        MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);

        srand(nseed);
        int nseeds[numprocs];
        for (int i = 0; i < numprocs; i++) {
            nseeds[i] = rand();
        }

        sprintf(filename, "%shoneycomb2-%d-%d-%d-%.2lf-%.2lf-%d-%d.out",
            filepath, Lx, Ly, LTau, dt, V, myid, nseeds[myid]);
        std::cout << "filename = " << filename << std::endl;
        ok = main_honeycomb_param(Lx, Ly, LTau, dt, V, nthreads, nseeds[myid], evaluationLength, filename);

    } else if (std::strcmp(argv[1], "--MRhoneycomb") == 0) {
        int Lx, Ly, LTau, nthreads, nseed, evaluationLength;
        char* filepath;
        double dt, V;
        char filename[100];

        Lx = std::stoi(argv[2]);
        Ly = std::stoi(argv[3]);
        LTau = std::stoi(argv[4]);
        dt = std::stod(argv[5]);
        V = std::stod(argv[6]);
        nthreads = std::stoi(argv[7]);
        nseed = std::stoi(argv[8]);
        evaluationLength = std::stoi(argv[9]);
        filepath = argv[10];

        int numprocs, myid;
        MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);

        srand(nseed);
        int nseeds[numprocs];
        for (int i = 0; i < numprocs; i++) {
            nseeds[i] = rand();
        }
        
        // filename = "honeycomb-Lx-Ly-LTau-dt-V-myid-seed.out"
        sprintf(filename, "%shoneycomb-%d-%d-%d-%.2lf-%.2lf-%d-%d.out", 
            filepath, Lx, Ly, LTau, dt, V, myid, nseeds[myid]);
        std::cout << "filename = " << filename << std::endl;
        ok = main_honeycombSingleMajorana(Lx, Ly, LTau, dt, V, nthreads, nseeds[myid], evaluationLength, filename);

        // fstream fin()
    } else if (std::strcmp(argv[1], "--MRhoneycomb4") == 0) {
        int Lx, Ly, LTau, nthreads, nseed, evaluationLength;
        char* filepath;
        double dt, V;
        char filename[100];

        Lx = std::stoi(argv[2]);
        Ly = std::stoi(argv[3]);
        LTau = std::stoi(argv[4]);
        dt = std::stod(argv[5]);
        V = std::stod(argv[6]);
        nthreads = std::stoi(argv[7]);
        nseed = std::stoi(argv[8]);
        evaluationLength = std::stoi(argv[9]);
        filepath = argv[10];

        int numprocs, myid;
        MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);

        srand(nseed);
        int nseeds[numprocs];
        for (int i = 0; i < numprocs; i++) {
            nseeds[i] = rand();
        }
        
        sprintf(filename, "%shoneycomb4-%d-%d-%d-%.2lf-%.2lf-%d-%d.out", 
            filepath, Lx, Ly, LTau, dt, V, myid, nseeds[myid]);
        std::cout << "filename = " << filename << std::endl;
        ok = main_honeycomb_4replica(Lx, Ly, LTau, dt, V, nthreads, nseeds[myid], evaluationLength, filename);

    } else if (std::strcmp(argv[1], "--MRhoneycomb4-swap") == 0) {
        int Lx, Ly, LTau, nthreads, nseed, evaluationLength;
        char* filepath;
        double dt, V;
        char filename[100];

        Lx = std::stoi(argv[2]);
        Ly = std::stoi(argv[3]);
        LTau = std::stoi(argv[4]);
        dt = std::stod(argv[5]);
        V = std::stod(argv[6]);
        nthreads = std::stoi(argv[7]);
        nseed = std::stoi(argv[8]);
        evaluationLength = std::stoi(argv[9]);
        filepath = argv[10];

        int numprocs, myid;
        MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);

        srand(nseed);
        int nseeds[numprocs];
        for (int i = 0; i < numprocs; i++) {
            nseeds[i] = rand();
        }
        
        sprintf(filename, "%shoneycomb4swap-%d-%d-%d-%.2lf-%.2lf-%d-%d.out", 
            filepath, Lx, Ly, LTau, dt, V, myid, nseeds[myid]);
        std::cout << "filename = " << filename << std::endl;
        // Call main_honeycomb_4replica with a flag or a new function?
        // Let's overload main_honeycomb_4replica to accept a boolean for trivial swap
        ok = main_honeycomb_4replica(Lx, Ly, LTau, dt, V, nthreads, nseeds[myid], evaluationLength, filename, true);

    } else if (std::strcmp(argv[1], "--chain") == 0) {
        int Lx, LTau, nthreads, nseed, evaluationLength, hsScheme, boundary;
        double dt, V, delta, mu;
        char* filepath;
        char filename[100];

        Lx = std::stoi(argv[2]);
        LTau = std::stoi(argv[3]);
        dt = std::stod(argv[4]);
        V = std::stod(argv[5]);
        delta = std::stod(argv[6]);
        nthreads = std::stoi(argv[7]);
        nseed = std::stoi(argv[8]);
        evaluationLength = std::stoi(argv[9]);
        filepath = argv[10];
        // if (argc <= 11) { 
        //     mu = 0.0;
        // } else {
        //     mu = std::stod(argv[11]);
        // }
        mu = std::stod(argv[11]);
        hsScheme = std::stoi(argv[12]);
        boundary = std::stoi(argv[13]);

        int numprocs, myid;
        MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);

        srand(nseed);
        int nseeds[numprocs];
        for (int i = 0; i < numprocs; i++) {
            nseeds[i] = rand();
        }
        
        // filename = "chain-Lx-Ltau-V-delta-myid-seed.out"
        sprintf(filename, "%schain-%d-%d-%.2lf-%.2lf-%.2lf-%d-%d-%.2lf-BC=%d.out", 
            filepath, Lx, LTau, dt, V, delta, myid, nseeds[myid], mu, boundary);
        std::cout << "filename = " << filename << std::endl;
        // if (argc == 11) { 
        //     ok = main_chain(Lx, LTau, dt, V, delta, nthreads, nseeds[myid], evaluationLength, filename);
        // } else {
        mu = std::stod(argv[11]);
        ok = main_chain(Lx, LTau, dt, V, delta, nthreads, nseeds[myid], evaluationLength, filename, mu, hsScheme, boundary);
        // }
    } else {
        std::cout << argv[1] << ": invalid arguments\n";
        ok = 0;
    }

    MPI_Finalize();
    return ok;
}
