#include "../inc/spinless_tV.h"
#include "../inc/types.h"
#include <iostream>
#include <vector>

// Mock Utils to avoid full dependency
class MockUtils : public SpinlessTvUtils {
public:
    MockUtils(int _Lx, int _Ly, double _dt, double _V, int _l, int _nDim)
        : SpinlessTvUtils(_Lx, _Ly, _dt, _V, _l, _nDim) {}
    
    void aux2MajoranaIdx(int idxAux, int imaj, int bType, int &idx1, int &idx2) const override {
        // Simple mapping for testing
        // idxAux is bond index. 
        // Let's say bond i connects 2*i and 2*i+1
        if (imaj == 0) {
            idx1 = 4 * idxAux;
            idx2 = 4 * idxAux + 1;
        } else {
            idx1 = 4 * idxAux + 2;
            idx2 = 4 * idxAux + 3;
        }
    }
};

int main() {
    int Lx = 2, Ly = 2;
    int nDim = 8; // Small dimension
    double dt = 0.1;
    double V = 1.0;
    int l = 10;
    
    MockUtils utils(Lx, Ly, dt, V, l, nDim);
    rdGenerator rd;
    
    // Create aux field
    iVecType* s = new iVecType(2); // 2 bonds
    (*s)(0) = 1;
    (*s)(1) = -1;
    
    SpinlessVOperator op(&utils, s, 0, &rd);
    
    // Get initial G
    MatType G_initial(nDim, nDim);
    op.getGreensMat(G_initial);

    std::cout << "Diagonal of G_initial:" << std::endl;
    for(int i=0; i<nDim; ++i) {
        std::cout << G_initial(i,i) << " ";
    }
    std::cout << std::endl;

    // Check consistency of G and B
    // G = 2(I+B)^-1 - I
    MatType B_mat = op.B; // Access B (public in SpinlessVOperator)
    MatType I = MatType::Identity(nDim, nDim);
    MatType G_from_B = 2.0 * (I + B_mat).inverse() - I;
    
    double consistency_diff = (G_initial - G_from_B).norm();
    std::cout << "Consistency diff (G vs 2(I+B)^-1 - I): " << consistency_diff << std::endl;


    
    // Perform single flip on bond 0
    DataType sign = 1.0;
    MatType G_updated = G_initial;
    
    // We force a flip to check the update logic
    // We need to access the update method. 
    // singleFlip does random check. We want deterministic update.
    // We can call updateSingleReplica directly if we make it public or friend.
    // Or we can just copy the logic here.
    
    // Let's copy the logic from updateSingleReplica for testing
    // to verify the formula with factor of 2.
    
    int replica = 0;
    int idxAux = 0;
    int offset = 0;
    int auxCur = (*s)(idxAux);
    int nDim_local = nDim;
    int inc = 1;
    
    // Logic from spinless_tV.h
    for (int imaj = 0; imaj < 2; imaj++) {
        int idx1, idx2;
        utils.aux2MajoranaIdx(idxAux, imaj, 0, idx1, idx2);
        int idx1_g = offset + idx1;
        int idx2_g = offset + idx2;
        
        DataType tmp_update = (1.0 - ((1.0i) * (utils.thlV) * double(auxCur) * G_updated(idx1_g, idx2_g)));
        
        // Update aux and B (simulated)
        // (*s)(idxAux) = -auxCur; // Don't update s yet, need it for next imaj loop? 
        // Wait, the loop updates s inside? 
        // In spinless_tV.h:
        // (*s[replica])(idxAux) = -auxCur;
        // This flips s. So for imaj=1, s is already flipped?
        // No, auxCur is local variable. But *s is modified.
        // If *s is modified, then for imaj=1, does it use new s?
        // auxCur is captured at start.
        // But inside loop: (*s[replica])(idxAux) = -auxCur;
        // So s is flipped in first iteration.
        // Second iteration: config->aux2MajoranaIdx uses idxAux.
        // Does aux2MajoranaIdx depend on s? No.
        // So it's fine.
        
        DataType alpha = (+1.0i) * double(auxCur) * (utils.thlV) / tmp_update;
        
        cVecType x1 = -G_updated.col(idx1_g);
        cVecType x2 = -G_updated.col(idx2_g);
        
        x1(idx1_g) += 1.0; // CHANGED TO 1.0 TO TEST
        x2(idx2_g) += 1.0;
        
        zgeru(&nDim_local, &nDim_local, &alpha, x1.data(), &inc, x2.data(), &inc, G_updated.data(), &nDim_local);
        
        alpha = -alpha;
        zgeru(&nDim_local, &nDim_local, &alpha, x2.data(), &inc, x1.data(), &inc, G_updated.data(), &nDim_local);
    }
    
    // Now flip s manually to calculate reference G
    (*s)(idxAux) = -auxCur;
    
    // Rebuild B and calculate G from scratch
    op.rebuildFullB(); // This uses the *s which is now flipped
    MatType G_ref(nDim, nDim);
    op.getGreensMat(G_ref);
    
    // Compare
    double diff = (G_updated - G_ref).norm();
    std::cout << "Difference: " << diff << std::endl;
    
    if (diff < 1e-10) {
        std::cout << "Test PASSED: Update formula matches brute force." << std::endl;
    } else {
        std::cout << "Test FAILED: Update formula incorrect." << std::endl;
    }
    
    return 0;
}
