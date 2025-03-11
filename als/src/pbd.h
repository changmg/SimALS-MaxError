/**
 * @file pbd.h
 * @author Chang Meng (chang.meng@epfl.ch)
 * @brief Partial Boolean difference (PBD) computation
 * @date 2025-02-21
 * 
 */


#pragma once


#include "my_abc.h"
// #include "lac.h"
// #include "cut.h"

namespace mecals_v1 {
class PBDMan {
private:
    NetMan net;
    std::vector<AbcObjVect> oldFos;
    std::vector<AbcObjVect> appFos;
    AbcObjVect oneCuts;
    std::vector<AbcObjVect> cutNtks;
    IntVect topoIds;
    DblVect flows;
public:
    explicit PBDMan() = default;
    ~PBDMan() = default;
    PBDMan(const PBDMan &) = delete;
    PBDMan(PBDMan &&) = delete;
    PBDMan & operator = (const PBDMan &) = delete;
    PBDMan & operator = (PBDMan &&)  = delete;

    void BuildMit(const NetMan& accNet, NetMan& appNet, NetMan& mitNet);
    void BuildPBD(double exactPBDPerc);
    int Synth(int fSASIMI);
    // int CheckWithSweepSASIMI(std::vector<std::pair<std::string, std::string>> & lacNames);
    // std::vector<std::pair<std::string, std::string>> CheckWithSATSASIMI();
    NetMan PostProc();
    AbcObjVect GetTFO(AbcObj* pObj) const;
    void GetTFORec(AbcObj* pObj, AbcObjVect& nodes) const;
    AbcObj* GetLocPBD(AbcObj* pV, AbcObj* pU);
};


// NetMan BuildMit(NetMan& accNet, NetMan& appNet, NetMan& mitNet);
} // namespace mecals_v1