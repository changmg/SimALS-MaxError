/**
 * @file error.h
 * @author Chang Meng
 * @brief Error analysis for approximate circuits
 * @date 2024-12-29
 * 
 */

#pragma once


#include "my_abc.h"
#include "simulator.h"
#include "lac.h"
#include "sat_wrapper.hpp"


/**
 * @brief Error metric type
 */
enum class METR_TYPE{
    ER, MED, MSE, MHD, MAXED, MAXHD
};
template <> 
struct fmt::formatter<METR_TYPE> {
    constexpr auto parse(fmt::format_parse_context& ctx) {
        return ctx.begin();
    }
    template <typename FormatContext>
    auto format(const METR_TYPE& value, FormatContext& ctx) {
        const std::string strs[6] = {"ER", "MED", "MSE", "MHD", "MAXED", "MAXHD"};
        return fmt::format_to(ctx.out(), "{}", strs[static_cast<int>(value)]);
    }
};
static inline METR_TYPE Str2MetrType(const std::string& str) {
    if (str == "ER")
        return METR_TYPE::ER;
    else if (str == "MED")
        return METR_TYPE::MED;
    else if (str == "MSE")
        return METR_TYPE::MSE;
    else if (str == "MHD")
        return METR_TYPE::MHD;
    else if (str == "MAXED")
        return METR_TYPE::MAXED;
    else if (str == "MAXHD")
        return METR_TYPE::MAXHD;
    else {
        fmt::print(stderr, "Error: unsupported metric type {}\n", str);
        assert(0);
    }
}


/**
 * @brief Error measurement manager
 */
class ErrMan {
private:
    const NetMan& net0;                        // accurate network
    const NetMan& net1;                        // approximate network
    std::shared_ptr<Simulator> pSmlt0;         // simulator for accurate network
    std::shared_ptr<Simulator> pSmlt1;         // simulator for approximate network
    std::shared_ptr<NetMan> pErrMit;           // error miter
    std::shared_ptr<CMSat::SATSolver> pSolver; // SAT solver
    IntVect cnfVarIdOfIthPi;                   // CNF variable ID of the i-th PI of the error miter

public:
    ErrMan(const NetMan& netMan0, const NetMan& netMan1);
    ErrMan(const NetMan& netMan0, const NetMan& netMan1, const NetMan& devNet);
    ~ErrMan() = default;
    ErrMan(const ErrMan&) = delete;
    ErrMan(ErrMan&&) = delete;
    ErrMan& operator = (const ErrMan&) = delete;
    ErrMan& operator = (ErrMan&&) = delete;

    void LogicSim(unsigned seed, int nFrame, DISTR_TYPE distrType);
    BigInt ComputeMaxErr(METR_TYPE metrType);
    std::shared_ptr<NetMan> BuildErrMit(const NetMan& accNet, const NetMan& appNet, const NetMan& devNet);
    std::shared_ptr<CMSat::SATSolver> BuildSatSolver_Naive(NetMan& net);
    std::shared_ptr<CMSat::SATSolver> BuildSatSolver_Abc(NetMan& net, IntVect& cnfVarIdOfIthPi);
    BigInt SolveSatsForMaxErrBinSearch(NetMan& net, CMSat::SATSolver& solver, int refEdWidth);
    void AddUnitClauseOfPi(int iPi, bool fVarCompl);

    inline CMSat::lbool SolveSat(bool printTime = false) {assert(pSolver != nullptr); return SolveSat_(*pSolver, printTime);}
    inline CMSat::lbool SolveSat(IntVect& counterExample, bool printTime = false) {assert(pSolver != nullptr); return SolveSatAndGetCountEx(*pSolver, std::vector<CMSat::Lit>{}, cnfVarIdOfIthPi, counterExample, printTime);}
    inline CMSat::lbool SolveSat(const std::vector<CMSat::Lit>& assumpts, IntVect& counterExample, bool printTime = false) {assert(pSolver != nullptr); return SolveSatAndGetCountEx(*pSolver, assumpts, cnfVarIdOfIthPi, counterExample, printTime);}
    inline const NetMan& GetErrMit() const {assert(pErrMit != nullptr); return *pErrMit;}
    inline double GetErrRate(unsigned seed, int nFrame, DISTR_TYPE distrType) {LogicSim(seed, nFrame, distrType); return pSmlt0->GetErrRate(*pSmlt1);}
    inline double GetMeanErrDist(unsigned seed, int nFrame, bool isSign, DISTR_TYPE distrType) {LogicSim(seed, nFrame, distrType); return pSmlt0->GetMeanErrDist(*pSmlt1, isSign);}
    inline ll GetMaxErrDistUsingEnum() {assert(net0.GetPiNum() < 20); LogicSim(0, 1ll << net0.GetPiNum(), DISTR_TYPE::ENUM); return pSmlt0->GetMaxErrDistFast(*pSmlt1);}
    inline void GetMaxErrDistLowBound(unsigned seed, int nFrame, BigInt& maxErrLowBound) {LogicSim(seed, nFrame, DISTR_TYPE::UNIF); pSmlt0->GetMaxErrDist(*pSmlt1, false, maxErrLowBound);}
    inline int GetCnfVarIdOfIthPi(int iPi) const {assert(iPi >= 0 && iPi < static_cast<int>(cnfVarIdOfIthPi.size())); return cnfVarIdOfIthPi[iPi];}
};


/**
 * @brief Batch error estimator for multiple LACs
 */
class BatchErrEst {
private:
    const int ROUGH_SIM_FRAME = 1024; // number of simulation frames for rough error estimation
    METR_TYPE metrType;               // error metric type
    unsigned seed;                    // random seed
    int nFrame;                       // number of simulation frames

public:
    explicit BatchErrEst(METR_TYPE metr_type, unsigned _seed, int n_frame): metrType(metr_type), seed(_seed), nFrame(n_frame) {}
    ~BatchErrEst() = default;
    BatchErrEst(const BatchErrEst &) = delete;
    BatchErrEst(BatchErrEst &&) = delete;
    BatchErrEst & operator = (const BatchErrEst &) = delete;
    BatchErrEst & operator = (BatchErrEst &&) = delete;

    // std::shared_ptr<LAC> PickFirstValidLac(LACMan& lacMan, Simulator& accSmlt, NetMan& appNet, const NetMan& devNet, std::unordered_set<std::string>& lacBlackList);
    void CompLacErrsByEnumAndPruneBadLacs(LACMan& lacMan, Simulator& accSmlt, const NetMan& appNet, ll errUppBound);
    void CompLacErrsBySimAndPruneBadLacs(LACMan& lacMan, Simulator& accSmlt, const NetMan& appNet, ll errUppBound);
    void PruneLacsWithSim(LACMan& lacMan, Simulator& accSmlt, const NetMan& appNet, ll errUppBound, int nFramePrune, DISTR_TYPE distrType = DISTR_TYPE::UNIF);
    void CalcErrLooseUppBound(LACMan& lacMan, const NetMan& appNet);
};


// functions to generate deviation network
std::shared_ptr<NetMan> GenDevNet(METR_TYPE metrType, int outWidth);
std::shared_ptr<NetMan> GenDevCompNet(METR_TYPE metrType, int outWidth);
std::shared_ptr<NetMan> GenDevCompNetEmbedErrBound(std::shared_ptr<NetMan> pDevNet, int outWidth, ll errUppBound);
