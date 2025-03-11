/**
 * @file lac.h
 * @author Chang Meng (chang.meng@epfl.ch)
 * @brief Local approximate change (LAC) management
 *
 */
#pragma once


#include "my_abc.h"
#include "header.h"
#include "simulator.h"


/**
 * @brief Local approximate change type
 */
enum class LAC_TYPE {
    CONSTANT, SASIMI, RESUB
};
template <> 
struct fmt::formatter<LAC_TYPE> {
    constexpr auto parse(fmt::format_parse_context& ctx) {
        return ctx.begin();
    }
    template <typename FormatContext>
    auto format(const LAC_TYPE& value, FormatContext& ctx) {
        const std::string strs[] = {"CONST", "SASIMI", "RESUB"};
        return fmt::format_to(ctx.out(), "{}", strs[static_cast<int>(value)]);
    }
};
static inline LAC_TYPE Str2LacType(const std::string& str) {
    if (str == "CONSTANT")
        return LAC_TYPE::CONSTANT;
    else if (str == "SASIMI")
        return LAC_TYPE::SASIMI;
    else if (str == "RESUB")
        return LAC_TYPE::RESUB;
    else {
        fmt::print(stderr, "Error: unsupported LAC type {}\n", str);
        assert(0);
    }
}


/**
 * @brief local approximate change (LAC)
 */
class LAC {
private:
    int targId;        // target node ID
    int sizeGain;      // (estimated) size reduction after applying the LAC
    double err;        // (estimated) error after applying the LAC
    double err2;       // (estimated) error 2 after applying the LAC
    IntVect divs;      // divisors used to construct the new function
    std::string sop;   // new function in SOP form on the divisors

public:
    explicit LAC(): targId(-1), sizeGain(-1), err(DBL_MAX), err2(DBL_MAX) {}
    explicit LAC(int targ_node_id, int gain, const IntVect& _divs, const std::string & _sop): targId(targ_node_id), sizeGain(gain), err(DBL_MAX), err2(DBL_MAX), divs(_divs), sop(_sop) {}
    virtual ~LAC() = default;
    LAC(const LAC & oth_lac) = default;
    LAC(LAC &&) = default;
    LAC & operator = (const LAC & oth_lac) = default;
    LAC & operator = (LAC &&) = default;

    inline int GetTargId() const {return targId;}
    inline void SetTargId(ll targ_node_id) {targId = targ_node_id;}
    inline int GetSizeGain() const {return sizeGain;}
    inline double GetErr() const {return err;}
    inline void SetErr(double _err) {err = _err;}
    inline double GetErr2() const {return err2;}
    inline void SetErr2(double _err) {err2 = _err;}
    inline const IntVect& GetDivIds() const {return divs;}
    inline const std::string& GetSop() const {return sop;}
    inline std::string ToStrShort() const {return fmt::format("n{}d{}f{}", targId, fmt::join(divs, ","), sop);}
    inline std::string ToStr() const {auto _sop = sop; std::replace(_sop.begin(), _sop.end(), '\n', ';'); return fmt::format("node {}, sizeGain = {}, err = {}, err2 = {}, divs = [{}], sop = [{}]", targId, sizeGain, err, err2, fmt::join(divs, ", "), _sop);}
    inline bool IsConst0() const {return divs.empty() && sop == " 0\n";}
    inline bool IsConst1() const {return divs.empty() && sop == " 1\n";}
    inline bool IsConst() const {return divs.empty() && (sop == " 0\n" || sop == " 1\n");}

};
using LACPtr = std::shared_ptr<LAC>;
using LACPtrVect = std::vector<LACPtr>;


/**
 * @brief Compare two LACs; define the order of pLac0 < pLac1
 * @brief primary key: smaller error; secondary key: larger size gain; third key: smaller target node ID
 * 
 * @param pLac0     the first LAC
 * @param pLac1     the second LAC
 * @return bool     true if pLac0 < pLac1
 */
static inline bool Lac0BetterThanLac1(const std::shared_ptr<LAC>& pLac0, const std::shared_ptr<LAC>& pLac1) {
    if (DoubleLess(pLac0->GetErr(), pLac1->GetErr()))
        return true;
    if (DoubleEqual(pLac0->GetErr(), pLac1->GetErr())) {
        if (pLac0->GetSizeGain() > pLac1->GetSizeGain())
            return true;
        // if (pLac0->GetSizeGain() == pLac1->GetSizeGain())
        //     return pLac0->GetTargId() < pLac1->GetTargId();
    }
    return false;
}


/**
 * @brief Compare two LACs; define the order of pLac0 < pLac1
 * @brief primary key: smaller err2; secondary key: smaller err; third key: larger size gain; fourth key: smaller target node ID
 * 
 * @param pLac0     the first LAC
 * @param pLac1     the second LAC
 * @return bool     true if pLac0 < pLac1
 */
static inline bool LacBetterThanUsingErr2(const std::shared_ptr<LAC>& pLac0, const std::shared_ptr<LAC>& pLac1) {
    // primary key: smaller err2
    if (DoubleLess(pLac0->GetErr2(), pLac1->GetErr2()))
        return true;
    if (DoubleEqual(pLac0->GetErr2(), pLac1->GetErr2())) {
        // secondary key: smaller err
        if (DoubleLess(pLac0->GetErr(), pLac1->GetErr()))
            return true;
        if (DoubleEqual(pLac0->GetErr(), pLac1->GetErr())) {
            // third key: larger size gain
            if (pLac0->GetSizeGain() > pLac1->GetSizeGain())
                return true;
            // if (pLac0->GetSizeGain() == pLac1->GetSizeGain())
            //     // fourth key: smaller target node ID
            //     return pLac0->GetTargId() < pLac1->GetTargId();
        }
    }
    return false;
}


/**
 * @brief LAC manager
 */
class LACMan {
private:
    LACPtrVect pLacs;                        // Local approximate changes
    std::unordered_map<int, LACPtrVect> node2Lacs;    // LACs grouped by node ID

public:
    explicit LACMan() = default;
    ~LACMan() = default;
    LACMan(const LACMan &) = delete;
    LACMan(LACMan &&) = delete;
    LACMan & operator = (const LACMan &) = delete;
    LACMan & operator = (LACMan &&) = delete;
    void GenConstLACs(const NetMan& net);
    void GenSasimiLACs(const NetMan& net, int maxCandResub, bool inclConst = true);
    void GenResubLACs(const NetMan& net, unsigned seed, int maxLevelDiff, int nFrame4ResubGen, int maxCandResub, bool inclConst = true);
    void RegroupLACsByNode(bool forceUpd = false);
    void GetDivs(abc::Abc_Obj_t* pNode, int nLevDivMax, IntVect& divs);
    void PrintLACs(int firstK = -1) const;
    std::shared_ptr<LAC> GetBestLac(double errUppBound = DBL_MAX) const;
    void SortAndKeepTopKLACs(int k);
    void RemLacsFromBlackList(std::unordered_set<std::string>& blackList);
    void RemLargeErrLacs(double errUppBound);

    inline int GetLacNum() const {return static_cast<int> (pLacs.size());}
    inline std::shared_ptr<LAC> GetLac(int i) const {return pLacs[i];}
    inline const LACPtrVect& GetLacs() const {return pLacs;}
    inline LACPtrVect GetLacsCopy() const {return pLacs;}
    inline const std::unordered_map<int, LACPtrVect>& GetNode2Lacs() const {return node2Lacs;}
    inline void ReplLacs(const LACPtrVect& newlacs) {pLacs.clear(); pLacs.assign(newlacs.begin(), newlacs.end());}
};


// miscellaneous functions
void ApplyLac(NetMan& net, const LAC& lac);
int TempApplyLac(NetMan& net, const LAC& lac, IntVect& replTrace, bool fVerb);
void TempApplyLacWithController(NetMan& net, const LAC& lac, int ctrlId, IntVect& replTrace, bool fVerb = false);
void RecovNet(NetMan& net, const Int2DVect& replTraces, bool fVerb);