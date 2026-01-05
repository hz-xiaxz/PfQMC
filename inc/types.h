#ifndef TYPES_H
#define TYPES_H

#ifndef EIGEN_USE_MKL_ALL
#define EIGEN_USE_MKL_ALL
#endif

#include<float.h>
#include<complex>
#include<random>
#include<iostream>

#define MKL_Complex16 std::complex<double>

#include "mkl_types.h"
#include "mkl.h"
#include "mkl_vsl.h"
#include <Eigen/Dense>

typedef std::complex<double> DataType;
using namespace std::complex_literals;
const DataType zero{0.0, 0.0};
const DataType one{1.0, 0.0};
// const DataType im{0.0, 1.0};
typedef Eigen::MatrixXcd MatType;
typedef Eigen::VectorXcd cVecType;
typedef Eigen::VectorXd dVecType;
typedef Eigen::VectorXi iVecType;
#define ind(i, j, N) (((j)*N) + i)
const double thresholdDBL = 1.0e-300;

// ============================================================
// Replica Constants and Helpers for 4-Replica PFQMC
// ============================================================

// Number of replicas for replica trick calculations
constexpr int N_REPLICAS = 4;

// Number of replica pairs (replicas 0-1 and 2-3 are paired)
constexpr int N_REPLICA_PAIRS = 2;

// Convert replica-local index to global index in the enlarged Green's function
// replica: replica index (0 to N_REPLICAS-1)
// localIdx: index within single replica (0 to nDimSingle-1)
// nDimSingle: dimension of single-replica Hilbert space
// Returns: global index in the 4*nDimSingle dimensional space
inline int replicaIdx(int replica, int localIdx, int nDimSingle) {
    return replica * nDimSingle + localIdx;
}

// Get the pair index for a given replica (0 or 1)
// Replicas 0,1 belong to pair 0; Replicas 2,3 belong to pair 1
inline int replicaPairIdx(int replica) {
    return replica / 2;
}

// Get the local index within a pair (0 or 1)
// Replica 0,2 -> 0; Replica 1,3 -> 1
inline int replicaLocalInPair(int replica) {
    return replica % 2;
}

// Get the offset for a replica pair in the global matrix
// Pair 0: offset 0; Pair 1: offset 2*nDimSingle
inline int pairOffset(int pair, int nDimSingle) {
    return pair * 2 * nDimSingle;
}

// ============================================================

inline DataType logDet(const MatType &H) {
    DataType ld = 0.0;
    Eigen::PartialPivLU<MatType> lu(H);
    auto& LU = lu.matrixLU();
    DataType c = lu.permutationP().determinant(); // -1 or 1
    // std::cout << "c=" << c << "\n";
    for (unsigned i = 0; i < LU.rows(); ++i) {
      const DataType& lii = LU(i,i);
    //   std::cout << lii << " q\n";
      ld += log(lii);
    }
    ld += log(c);
    return ld;
}

class rdGenerator {

private:
    std::uniform_int_distribution<int> rdDistZ2;
    std::uniform_real_distribution<double> rdDistUniform01;
    std::normal_distribution<double> rdDistNormal;
    std::mt19937 rdEng;
public:
    rdGenerator(int seed=114514) {
        // random generator, Z_2 auxillary field
        rdDistZ2 = std::uniform_int_distribution<int>(0, 1);
        rdDistUniform01 = std::uniform_real_distribution<double>(0.0, 1.0);
        rdDistNormal = std::normal_distribution<double>(0.0, 1.0);
        rdEng.seed(seed);
    }

    inline int rdZ2() {
        return 2*rdDistZ2(rdEng) - 1;
    }

    inline double rdUniform01() {
        // double r = rdDistUniform01(rdEng);
        // return r;
        return rdDistUniform01(rdEng);
    }

    inline double rdNormal() {
        return rdDistNormal(rdEng);
    }

};

#endif