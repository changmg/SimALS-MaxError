/**
 * @file simulator.h
 * @author Chang Meng (chang.meng@epfl.ch)
 * @brief Circuit simulator
 * @date 2025-01-25
 * 
 */
#pragma once


#include "header.h"
#include "my_abc.h"


/**
 * @brief Distribution type
 */
enum class DISTR_TYPE {
    UNIF, ENUM
};
template <> 
struct fmt::formatter<DISTR_TYPE> {
    constexpr auto parse(fmt::format_parse_context& ctx) {
        return ctx.begin();
    }
    template <typename FormatContext>
    auto format(const DISTR_TYPE& value, FormatContext& ctx) {
        const std::string strs[2] = {"UNIFORM", "ENUMERATION"};
        return fmt::format_to(ctx.out(), "{}", strs[static_cast<int>(value)]);
    }
};


/**
 * @brief Get the value encoded by the binary pattern dat[iPatt] (fast version)
 * @brief Assume the value is an unsigned integer
 * 
 * @param dat          the data
 * @param iPatt        the pattern index
 * @return long long   the value
 */
static inline ll GetValueFast(const std::vector<BitVect>& dat, int iPatt) {
    int msb = static_cast<int>(dat.size()) - 1;
    assert(msb < 63);
    ll ret = 0;
    for (ll k = msb; k >= 0; --k) {
        ret <<= 1;
        ret |= dat[k][iPatt];
    }
    return ret;
}


/**
 * @brief Get the value encoded by the binary pattern dat[iPatt]
 * @brief Assume the value is an unsigned integer
 * 
 * @param dat          the data
 * @param iPatt        the pattern index
 * @param ret          the value
 * @retval ret         the value
 * @return void
 */
static inline void GetValue(const std::vector<BitVect>& dat, int iPatt, BigInt& ret) {
    int msb = static_cast<int>(dat.size()) - 1;
    assert(msb < 500);
    ret = 0;
    for (ll k = msb; k >= 0; --k) {
        ret <<= 1;
        ret |= dat[k][iPatt];
    }
}


/**
 * @brief Simulator for a circuit network
 */
class Simulator: public NetMan {
private:
    unsigned seed;                  // random seed for generating input patterns
    int nFrame;                     // number of simulation frames (patterns)
    DISTR_TYPE distrType;           // input distribution type
    std::vector<BitVect> dat;       // dat[pObj->Id], simulation patterns for the node pObj
    std::vector<BitVect> tempDat;   // tempDat[pObj->Id], temporary simulation patterns for the node pObj

public:
    explicit Simulator(const NetMan& net_man, unsigned _seed, int n_frame, DISTR_TYPE distr_type = DISTR_TYPE::UNIF);
    ~Simulator() = default;
    Simulator(const Simulator&) = delete;
    Simulator(Simulator&&) = delete;
    Simulator& operator = (const Simulator&) = delete;
    Simulator& operator = (Simulator&&) = delete;

    void InitConstNodes();
    void GenInpUnif();
    void GenInpUnifFast();
    void GenInpEnum();
    void GenInpFromOthSmlt(const Simulator& othSmlt);
    void GenInpFromBitVects(const std::vector<BitVect>& piPatts);
    void ReplInp(int iPatt, const IntVect& piVals);
    void AppendInp(const IntVect& piVals);
    void UpdNodeAndPoPatts();
    void SimSop(const IntVect& faninIds, const std::string& sop, BitVect& res);
    void UpdAigNode(AbcObj* pObj);
    void UpdSopNode(AbcObj* pObj);
    void UpdGateNode(AbcObj* pObj);
    void UpdStrashNode(AbcObj* pObj);
    void UpdSop(AbcObj* pObj, char* pSop);
    BigInt GetInput(int iPatt, int lsb, int msb) const;
    void GetInpVect(int iPatt, IntVect& piVals) const;
    void PrintInpStream(int iPatt) const;
    void GetOutput(int iPatt, BigInt& ret) const;
    ll GetOutputFast(int iPatt) const;
    ll GetTempOutputFast(int iPatt) const;
    void PrintOutpStream(int iPatt) const;
    double GetSignalProb(int objId) const;
    void PrintSignalProb() const;
    double GetErrRate(const Simulator& oth_smlt, bool isCheck = false) const;
    double GetMeanErrDist(const Simulator& oth_smlt, bool isCheck = false) const;
    ll GetMaxErrDistFast(const Simulator& oth_smlt, bool isCheck = false) const;
    void GetMaxErrDist(const Simulator& oth_smlt, bool isCheck, BigInt& maxErrLowBound) const;
    void CalcBoolDiff(const AbcObjVect& topoNodes, int idx, std::vector<BitVect>& bdPosWrtNode);
    void CalcLocBoolDiff(AbcObj* pObj, std::list <AbcObj*>& disjCut, std::vector <AbcObj*>& cutNtk, std::vector <BitVect>& bdCut2Node);
    void UpdSopNodeForBoolDiff(AbcObj* pObj);
    void UpdGateNodeForBoolDiff(AbcObj* pObj);
    void UpdSopForBoolDiff(AbcObj* pObj, char* pSop);
    void PrintDat() const;

    inline void GenInpPatts() {if (distrType == DISTR_TYPE::UNIF) GenInpUnifFast(); else if (distrType == DISTR_TYPE::ENUM) GenInpEnum(); else assert(0);}
    inline int GetFrameNumb() const {return nFrame;}
    inline void LogicSim() {GenInpPatts(); UpdNodeAndPoPatts();}
    inline const BitVect& GetDat(int id) const {return dat[id];}
    inline void SetPiConst(int ithPi, int const1) {if (const1) dat[GetPiId(ithPi)].set(); else dat[GetPiId(ithPi)].reset();}
};