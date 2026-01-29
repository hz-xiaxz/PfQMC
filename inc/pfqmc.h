#ifndef PFQMC_H
#define PFQMC_H

#include "spinless_tV.h"
#include "qr_udt.h"

class PfQMC
{
public:
    int stb;
    int nReplicas;        // NEW: number of replicas (default 1)
    int nDimSingle;       // NEW: single-replica dimension
    int nDim;             // CHANGED: = nReplicas * nDimSingle

    MatType g;            // nDim × nDim Green's function

    std::vector<Operator *> op_array;
    Operator* mixingOp;   // NEW: mixing operator at τ=0 (nullptr if nReplicas==1)

    int op_length;
    std::vector<bool> need_stabilization;
    int checkpoints;
    std::vector<UDT> udtL;
    std::vector<UDT> udtR;

    DataType sign;
    
    long long acceptedUpdates = 0;
    long long totalUpdates = 0;

    PfQMC(Spinless_tV *walker, int _stb = 10, int _nReplicas = 1, Operator* _mixingOp = nullptr);

    void rightInit()
    {
        MatType tmp = MatType::Identity(nDim, nDim);
        MatType Aseg = MatType::Identity(nDim, nDim);
        int curSeg = 0;
        
        // Apply mixing operator first if it exists
        UDT mixUDT;
        bool hasMix = (mixingOp != nullptr);
        if (hasMix) {
            mixingOp->left_multiply(Aseg, tmp);
            std::swap(Aseg, tmp);
            
            // Immediate UDT
            int nBlocks = (nReplicas > 1) ? 2 : 1;
            mixUDT = UDT(Aseg, nBlocks);
            Aseg = MatType::Identity(nDim, nDim);
        }

        for (int l = 0; l < op_length; l++)
        {
            op_array[l]->left_multiply(Aseg, tmp);
            std::swap(Aseg, tmp);
            // %op_length is important, cannot be remove. or else the last segment will not be calculated
            if (need_stabilization[(l + 1) % op_length])
            {
                if (curSeg == 0)
                {
                    int nBlocks = 1;
                    if (nReplicas > 1) {
                        nBlocks = (mixingOp != nullptr) ? 2 : nReplicas;
                    }
                    udtR[curSeg] = UDT(Aseg, nBlocks); // TODO: performance check
                    if (hasMix) udtR[curSeg] = udtR[curSeg] * mixUDT;
                }
                else
                {
                    udtR[curSeg] = Aseg * udtR[curSeg - 1];
                }
                Aseg = MatType::Identity(nDim, nDim);
                curSeg++;
            }
        }
        udtR[checkpoints - 1].onePlusInv(g);
    }

    void leftInit()
    {
        MatType tmp = MatType::Identity(nDim, nDim);
        MatType Aseg = MatType::Identity(nDim, nDim);
        int curSeg = checkpoints - 1;
        for (int l = op_length - 1; l > -1; l--)
        {
            op_array[l]->right_multiply(Aseg, tmp);
            std::swap(Aseg, tmp);
            
            // Apply mixing operator at l=0 (end of time evolution in this direction)
            if (l == 0 && mixingOp != nullptr) {
                mixingOp->right_multiply(Aseg, tmp);
                std::swap(Aseg, tmp);
            }

            if (need_stabilization[l])
            {
                Aseg.adjointInPlace();
                if (curSeg == (checkpoints - 1))
                {
                    int nBlocks = 1;
                    if (nReplicas > 1) {
                        nBlocks = (mixingOp != nullptr) ? 2 : nReplicas;
                    }
                    udtL[curSeg] = UDT(Aseg, nBlocks); // TODO: performance check
                }
                else
                {
                    udtL[curSeg] = Aseg * udtL[curSeg + 1];
                }
                Aseg = MatType::Identity(nDim, nDim);
                curSeg--;
            }
        }
        udtL[0].onePlusInv(g);
        g.adjointInPlace();
    }

    // after each sweep
    // the greens function g
    // is automatically updated
    void rightSweep();
    void leftSweep();

    // get sign by computing the pfaffian of
    // a 4N * 4N matrix
    DataType getSignRaw();

    // should provide same result as getSignRaw
    // but by computing the pfaffian of a 2N * 2N matrix
    // TODO: this method currently has fundamental flaws
    // therefore should not be used
    // DataType getSign();
    
    // ~PfQMC()
    // {
        // for (int i = 0; i < udtR.size(); i++)
        // {
        //     delete udtR[i];
        // }
        // delete Al;
    // }
};

#endif