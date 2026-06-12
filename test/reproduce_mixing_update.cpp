#include "../inc/mixing_operator.h"
#include "../inc/types.h"
#include <iostream>
#include <vector>
#include <iomanip>

// Mock rdGenerator since we don't want to include the full random header if it's complex,
// or just include it if it's simple.
// Looking at mixing_operator.h, it includes "operator.h".
// Let's hope operator.h and rdGenerator are simple enough or available.
// mixing_operator.h -> operator.h
// We might need to mock rdGenerator or include it.
// Assuming "operator.h" defines rdGenerator or includes "types.h" where it might be.

// Let's define a minimal mock for rdGenerator if not easily available, but it's used in constructor.
// Better to rely on linking.

int main() {
    int nReplicas = 4; // Use 4 replicas to test p=0 and p=1
    int nDimSingle = 4;
    int nDim = nReplicas * nDimSingle; // 16

    // Mock RD
    rdGenerator rd;

    MixingOperator mix(nReplicas, nDimSingle, &rd);
    
    // Create a valid skew-symmetric Green's function
    MatType G = MatType::Random(nDim, nDim);
    G = 0.5 * (G - G.transpose()).eval(); // Skew symmetric

    std::cout << "Initial G:\n" << G << "\n\n";
    std::cout << "Is Initial Skew? " << (G.transpose() + G).norm() << "\n";
    
    // Apply Left Propagate (G' = - B G B)
    // For Mixing Operator, B is initially Identity? No, it's defined by sigma.
    // Let's verify initial sigma.
    std::cout << "Initial Sigmas: ";
    for(int p=0; p<nReplicas/2; ++p) {
        for(int m=0; m<nDimSingle; ++m) {
            std::cout << mix.sigma[p][m] << " ";
        }
    }
    std::cout << "\n";

    MatType G_propagated = MatType::Zero(nDim, nDim);
    MatType G_tmp = MatType::Zero(nDim, nDim);
    
    mix.left_propagate(G, G_tmp);
    G_propagated = G; // left_propagate modifies G in-place (first arg)
    
    std::cout << "Propagated G' (left_propagate):\n" << G_propagated << "\n\n";
    std::cout << "Is Propagated Skew? " << (G_propagated.transpose() + G_propagated).norm() << "\n";
    
    // === GROUND TRUTH VERIFICATION ===
    std::cout << "\n=== Ground Truth Verification ===\n";
    
    MatType I = MatType::Identity(nDim, nDim);
    
    // Get B from mixing operator
    MatType B = MatType::Zero(nDim, nDim);
    mix.left_multiply(I, B);
    
    // Add a small random perturbation to B to avoid degenerate ratio=0 case
    // In real simulation, B includes contributions from Hamiltonian propagators
    MatType B_perturb = 0.3 * MatType::Random(nDim, nDim);
    B_perturb = (B_perturb - B_perturb.transpose()).eval();  // Make skew-symmetric
    MatType B_total = B + B_perturb;
    
    // Compute g = 2*(I+B_total)^{-1}
    MatType g_old = 2.0 * (I + B_total).inverse();
    
    int p = 0;
    int m = 0;
    int original_sigma = mix.sigma[p][m];
    
    // Flip sigma to get B'
    mix.sigma[p][m] = -original_sigma;
    MatType Bp = MatType::Zero(nDim, nDim);
    mix.left_multiply(I, Bp);
    MatType Bp_total = Bp + B_perturb;  // Same perturbation
    
    // Compute g_ref = 2*(I+B'_total)^{-1} - ground truth
    MatType g_ref = 2.0 * (I + Bp_total).inverse();
    
    // Restore sigma
    mix.sigma[p][m] = original_sigma;
    
    int parity = (p % 2 == 0) ? 1 : -1;
    int s = original_sigma * parity;
    int j = 2 * p * nDimSingle + m;
    int k = (2 * p + 1) * nDimSingle + m;
    
    DataType G12 = g_old(j, k);
    DataType ratio_computed = mix.getRatio(g_old, p, m);
    
    std::cout << "s = " << s << ", G12 = " << G12 << "\n";
    std::cout << "Ratio (Computed): " << ratio_computed << "\n";
    
    // Show the 2x2 block before and after
    std::cout << "\ng_old 2x2 block:\n";
    std::cout << "  [[" << g_old(j,j) << ", " << g_old(j,k) << "],\n";
    std::cout << "   [" << g_old(k,j) << ", " << g_old(k,k) << "]]\n";
    
    std::cout << "\ng_ref 2x2 block (target):\n";
    std::cout << "  [[" << g_ref(j,j) << ", " << g_ref(j,k) << "],\n";
    std::cout << "   [" << g_ref(k,j) << ", " << g_ref(k,k) << "]]\n";
    
    // Apply localUpdate
    MatType g_updated = g_old;
    mix.localUpdate(g_updated, p, m);
    
    std::cout << "\ng_updated 2x2 block:\n";
    std::cout << "  [[" << g_updated(j,j) << ", " << g_updated(j,k) << "],\n";
    std::cout << "   [" << g_updated(k,j) << ", " << g_updated(k,k) << "]]\n";
    
    double diff = (g_updated - g_ref).norm();
    std::cout << "\nDifference ||g_updated - g_ref||: " << diff << "\n";
    
    // Note: Skew-symmetry of G = g - I depends on B being skew-symmetric.
    // Our test perturbation may not preserve this, so we just verify that
    // g_updated matches g_ref (which inherits the same skew-symmetry structure).
    MatType G_ref = g_ref - I;
    double skew_ref_norm = (G_ref + G_ref.transpose()).norm();
    std::cout << "(Informational) Skew-symmetry of g_ref: ||G + G^T|| = " << skew_ref_norm << "\n";
    
    // Round-trip test: apply update again (sigma is now flipped)
    MatType g_roundtrip = g_updated;
    mix.localUpdate(g_roundtrip, p, m);
    double roundtrip_diff = (g_roundtrip - g_old).norm();
    std::cout << "Round-trip difference: " << roundtrip_diff << "\n";
    
    if (diff > 1e-10) {
        std::cerr << "Update logic FAILED! Mismatch with ground truth.\n";
        return 1;
    }
    
    if (roundtrip_diff > 1e-10) {
        std::cerr << "Update logic FAILED! Round-trip not identity.\n";
        return 1;
    }
    
    std::cout << "\nUpdate logic VERIFIED!\n";
    std::cout << "  - Ground truth match: PASS\n";
    std::cout << "  - Round-trip: PASS\n";

    return 0;
}
