/**
 * @file my_abc.h
 * @author Chang Meng (chang.meng@epfl.ch)
 * @brief Wrapper for ABC, the open-source logic synthesis and verification system
 * @date 2025-01-25
 * 
 */
#pragma once


#include "header.h"
#include "my_util.hpp"
namespace abc {
#include "aig/aig/aig.h"
#include "aig/gia/gia.h"
#include "aig/hop/hop.h"
#include "base/main/main.h"
#include "base/main/mainInt.h"
#include "base/cmd/cmd.h"
#include "base/io/ioAbc.h"
#include "base/abc/abc.h"
#include "bool/bdc/bdc.h"
#include "bool/kit/kit.h"
#include "bool/dec/dec.h"
#include "misc/nm/nm.h"
#include "misc/nm/nmInt.h"
#include "misc/util/abc_global.h"
#include "misc/util/util_hack.h"
#include "map/mio/mio.h"
#include "map/mio/mioInt.h"
#include "map/mapper/mapper.h"
#include "map/mapper/mapperInt.h"
#include "map/scl/scl.h"
#include "map/scl/sclCon.h"
#include "map/scl/sclSize.h"
#include "opt/cut/cut.h"
#include "opt/cut/cutInt.h"
#include "opt/cut/cutList.h"
#include "opt/mfs/mfs.h"
#include "opt/mfs/mfsInt.h"
#include "opt/sim/sim.h"
#include "opt/rwr/rwr.h"
#include "proof/fraig/fraig.h"
#include "proof/ssw/ssw.h"


/**
 * @brief ABC time manager
 * @brief Copied from ABC source code
 * @brief Used to manage the timing information of the network
 * 
 */
struct Abc_ManTime_t_ {
    Abc_Time_t     tArrDef;
    Abc_Time_t     tReqDef;
    Vec_Ptr_t  *   vArrs;
    Vec_Ptr_t  *   vReqs;
    Abc_Time_t     tInDriveDef;
    Abc_Time_t     tOutLoadDef;
    Abc_Time_t *   tInDrive;
    Abc_Time_t *   tOutLoad;
};
}


/**
 * @brief Short-hand for ABC object
 * 
 */
using AbcObj = abc::Abc_Obj_t;
using AbcObjList = std::list<AbcObj*>;
using AbcObjSet = std::unordered_set<AbcObj*>;
using AbcObjVect = std::vector<AbcObj*>;
using AbcObjQueue = std::queue<AbcObj*>;
using AbcObjPair = std::pair<AbcObj*, AbcObj*>;


/**
 * @brief Network type
 * 
 */
enum class NET_TYPE {
    AIG, GATE, SOP, STRASH
};
template <> 
struct fmt::formatter<NET_TYPE> {
    constexpr auto parse(fmt::format_parse_context& ctx) {
        return ctx.begin();
    }
    template <typename FormatContext>
    auto format(const NET_TYPE& value, FormatContext& ctx) {
        const std::string strs[] = {"AIG", "GATE", "SOP", "STRASH"};
        assert(static_cast<int>(value) >= 0 && static_cast<int>(value) <= 3);
        return fmt::format_to(ctx.out(), "{}", strs[static_cast<int>(value)]);
    }
};


/**
 * @brief Optimization orientation
 * 
 */
enum class ORIENT {
    AREA, DELAY
};
static inline std::ostream& operator << (std::ostream& os, const ORIENT orient) {
    const std::string strs[2] = {"AREA", "DELAY"};
    os << strs[static_cast <int> (orient)];
    return os;
}


/**
 * @brief Mapping cell type
 * 
 */
enum class MAP_TYPE {
    LUT, SCL
};


/**
 * @brief Quality improvement flag
 * 
 */
enum class IMPR {
    GOOD, BAD, SAME
};


/**
 * @brief ABC manager, a wrapper of ABC command line interface
 */
class AbcMan {
public:
    explicit AbcMan() {
        if (abc::Abc_FrameGetGlobalFrame() == nullptr) {
            fmt::print(stderr, "Error: ABC global frame is not initialized\n");
            assert(0);
        }
    }
    ~AbcMan() = default;
    AbcMan(const AbcMan&) = delete;
    AbcMan(AbcMan&&) = delete;
    AbcMan& operator = (const AbcMan&) = delete;
    AbcMan& operator = (AbcMan&&) = delete;

    void Comm(const std::string& cmd, bool fVerb = false);
    void ReadNet(const std::string& fileName);
    void WriteNet(const std::string& fileName, bool fVerb = false);
    void ReadStandCell(const std::string& fileName);
    void ConvToAig();
    void ConvToGate();
    void ConvToSop();
    void Strash();
    void PrintStat();
    void TopoSort();
    void StatTimeAnal();
    void Synth(ORIENT orient, bool fVerb = false);
    void SynthAndMap(double maxDelay, bool fVerb = false);
    std::pair <double, double> Map(MAP_TYPE targ, ORIENT orient, bool fVerb = false);
    IMPR UpdNet(double oldArea, double oldDelay, abc::Abc_Ntk_t * oldNtk, double newArea, double newDelay, ORIENT orient);
    NET_TYPE GetNetType(abc::Abc_Ntk_t * pNtk) const;
    double GetArea(abc::Abc_Ntk_t * pNtk) const;
    double GetDelay(abc::Abc_Ntk_t * pNtk) const;
    bool CheckSCLNet(abc::Abc_Ntk_t * pNtk) const;
    AbcObj* GetTwinNode( AbcObj* pNode );
    void LoadAlias();

    inline abc::Abc_Frame_t* GetAbcFrame() const {return abc::Abc_FrameGetGlobalFrame();}
    inline abc::Abc_Ntk_t* GetNet() const {return abc::Abc_FrameReadNtk(GetAbcFrame());}
    inline NET_TYPE GetNetType() const {return GetNetType(GetNet());}
    inline double GetArea() const {return GetArea(GetNet());}
    inline double GetDelay() const {return GetDelay(GetNet());}
    inline bool IsLutNet() const {const int LUT_INP = 6; return GetNetType() != NET_TYPE::GATE && Abc_NtkGetFaninMax(GetNet()) <= LUT_INP;}
    inline void SetMainNet(abc::Abc_Ntk_t * pNtk) {assert(pNtk != nullptr); if (pNtk != GetNet()) Abc_FrameReplaceCurrentNetwork(GetAbcFrame(), pNtk);}
};


/**
 * @brief Network manager 
 * @brief A network means a logic circuit using the ABC-based network structure
 */
class NetMan: public AbcMan {
private:
    abc::Abc_Ntk_t* pNtk;         // ABC network
    bool isDupl;                  // who manages the memory of the ABC network pNtk; if true, NetMan manages; if false, ABC global frame manages

public:
    explicit NetMan();
    explicit NetMan(abc::Abc_Ntk_t * p_ntk, bool is_dupl = false);
    explicit NetMan(const std::string& fileName);
    ~NetMan();
    NetMan(const NetMan& net_man);
    NetMan(NetMan&& net_man);
    NetMan& operator = (const NetMan& net_man);
    NetMan& operator = (NetMan&& net_man);

    // print & check
    void Print(bool showFunct = false) const;
    void PrintObjBas(AbcObj* pObj, std::string&& endWith) const;
    void PrintObj(AbcObj* pObj, bool showFunct = false) const;
    bool IsPIOSame(const NetMan& oth_net) const;
    void WriteBlif(const std::string& fileName) const;
    void WriteDot(const std::string& fileName) const;

    // get sub-structure
    AbcObjVect CalcTopoOrd(bool inclConst = true) const;
    void CalcTopoOrdRec(AbcObj* pObj, AbcObjVect& nodes, bool inclConst) const; 
    IntVect CalcTopoOrdOfIds(bool inclConst = true) const;
    void CalcTopoOrdOfIdsRec(AbcObj* pObj, IntVect& nodes, bool inclConst) const; 
    AbcObjVect GetTFI(AbcObj* pObj) const;
    void GetTFIRec(AbcObj* pObj, AbcObjVect& nodes) const;
    IntVect GetTFI(AbcObj* pObj, const std::set <int>& critGraph) const;
    void GetTFIRec(AbcObj* pObj, IntVect& objs, const std::set <int>& critGraph) const;
    AbcObjVect GetTFO(AbcObj* pObj) const;
    void GetTFORec(AbcObj* pObj, AbcObjVect& nodes) const;
    IntVect GetTFO(AbcObj* pObj, const std::set <int>& critGraph) const;
    void GetTFORec(AbcObj* pObj, IntVect& objs, const std::set <int>& critGraph) const;
    int GetNodeMffcSize(AbcObj* pNode) const;
    AbcObjVect GetNodeMffc(AbcObj* pNode) const;

    // modify network
    IntPair GetConstIds(bool fVerb = false);
    IntPair CreateConstsIfNotExist(bool fVerb = false);
    void MergeConst(bool fVerb = false);
    void Comm(const std::string& cmd, bool fVerb = false);
    void Synth(ORIENT orient, bool fVerb = false);
    void SynthAndMap(double maxDelay, bool fVerb = false);
    void TempRepl_v2(int tsId, int ssId, IntVect& replTrace, bool fVerb = false);
    void Recov_v2(const IntVect& replTrace, bool fVerb = false);
    void PatchFanin(AbcObj* pObj, int iFanin, AbcObj* pFaninOld, AbcObj* pFaninNew);
    void Trunc(int truncBit);
    AbcObj* CreateNode(const AbcObjVect& pFanins, const std::string& sop);
    int CreateNode(const IntVect& faninIds, const std::string& sop);
    int CreateAIGStyleNodes(const IntVect& faninIds, const std::string& sop);
    AbcObj* CreateGate(AbcObjVect&& fanins, const std::string& gateName);
    AbcObj* DupObj(AbcObj* pObj, const char* pSuff);
    void LimFanout();
    void ReplaceByComplementedObj(int targId, int subId);
    void PropConst(int startId, bool fKeepDanglNodes, bool fVerb);
    void PropConst(bool fVerb = false);

    inline int CleanUp(bool fVerb = false) {return abc::Abc_NtkCleanup(pNtk, fVerb);}
    inline void Sweep(bool fVerb = false) {CleanUp(fVerb); PropConst(fVerb);}
    
    // get properties of network
    int NodeDeref_rec(AbcObj* pRoot, AbcObj* pNode, std::unordered_set<AbcObj*>& divSet) const;
    int NodeRef_rec(AbcObj* pRoot, AbcObj* pNode, std::unordered_set<AbcObj*>& divSet) const;
    int NodeDeref_rec_v2(int rootId, int nodeId, IntVect& delNodes) const;
    int NodeRef_rec_v2(int rootId, int nodeId) const;
    int GetSizeGain(int rootId, const IntVect& divIds) const;
    int GetSizeGain(const IntVect& targetIds, const IntVect& divIds) const;

    inline abc::Abc_Ntk_t* GetNet() const {return pNtk;}
    inline std::string GetNetName() const {if (pNtk->pName != nullptr) return std::string(pNtk->pName); else return "(null)";}
    inline NET_TYPE GetNetType() const {return AbcMan::GetNetType(pNtk);}
    inline bool IsStrash() const {return GetNetType() == NET_TYPE::STRASH;}
    inline int Check() const {return abc::Abc_NtkDoCheck(pNtk);} // return 1 if the network is okay
    inline int IsAcyclic() const {return abc::Abc_NtkIsAcyclic(pNtk);}
    inline double GetArea() const {return AbcMan::GetArea(pNtk);}
    inline double GetDelay() const {return AbcMan::GetDelay(pNtk);}
    inline void WriteNet(const std::string& fileName, bool fVerb = false) {AbcMan::SetMainNet(abc::Abc_NtkDup(pNtk)); AbcMan::WriteNet(fileName, fVerb);}
    inline void WriteNet(const std::string&& fileName, bool fVerb = false) {AbcMan::SetMainNet(abc::Abc_NtkDup(pNtk)); AbcMan::WriteNet(fileName, fVerb);}
    inline void PrintStat() {AbcMan::SetMainNet(abc::Abc_NtkDup(pNtk)); AbcMan::PrintStat();}
    inline bool IsInTopoOrd() const {auto type = GetNetType(); assert(type == NET_TYPE::AIG || type == NET_TYPE::GATE || type == NET_TYPE::SOP); return abc::Abc_SclCheckNtk(pNtk, 0);}
    inline bool CheckSCLNet() const {return AbcMan::CheckSCLNet(pNtk);}
    inline void CollMffc(int rootId, IntVect& mffcNodes) const {mffcNodes.clear(); int count0 = NodeDeref_rec_v2(rootId, rootId, mffcNodes); int count1 = NodeRef_rec_v2(rootId, rootId); assert(count0 == count1);}

    // convert network representation
    inline void ConvToSop() {abc::Abc_NtkToSop(pNtk, -1, ABC_INFINITY);}
    inline void Strash() {assert(isDupl); auto pNtkAig = abc::Abc_NtkStrash(pNtk, 0, 1, 0); abc::Abc_NtkDelete(pNtk); pNtk = pNtkAig;}

    // get PIs & POs & nodes
    inline int GetPiNum() const {return abc::Abc_NtkPiNum(pNtk);}
    inline int GetObjNumMax() const {return abc::Abc_NtkObjNumMax(pNtk);}
    inline int GetObjNum() const {return abc::Abc_NtkObjNum(pNtk);}
    inline int GetPoNum() const {return abc::Abc_NtkPoNum(pNtk);}
    inline int GetNodeNum() const {return abc::Abc_NtkNodeNum(pNtk);}
    inline AbcObj* GetPi(int i) const {return abc::Abc_NtkPi(pNtk, i);}
    inline AbcObj* GetObj(int i) const {return abc::Abc_NtkObj(pNtk, i);}
    inline AbcObj* GetPo(int i) const {return abc::Abc_NtkPo(pNtk, i);}
    inline int GetIdMaxPlus1() const {return abc::Abc_NtkObjNumMax(pNtk);}
    inline int GetIdMax() const {return abc::Abc_NtkObjNumMax(pNtk) - 1;}
    inline int GetId(AbcObj* pObj) const {return abc::Abc_ObjId(pObj);}
    inline int GetPiId(int i) const {return GetId(GetPi(i));}
    inline int GetPoId(int i) const {return GetId(GetPo(i));}
    inline AbcObj* GetPoDriv(int i) const {return abc::Abc_ObjFanin0(GetPo(i));}
    inline int GetPoDrivId(int i) const {return GetId(GetPoDriv(i));}
    inline AbcObj* GetNodeByName(const std::string& name) const {return abc::Abc_NtkFindNode(pNtk, const_cast<char *>(name.c_str()));}
    inline AbcObj* GetPiByName(const std::string& name) const {return abc::Abc_NtkFindCi(pNtk, const_cast<char *>(name.c_str()));}
    inline int GetConst1IdInStrashNet() const {assert(IsStrash()); return abc::Abc_AigConst1(pNtk)->Id;}

    // get/set properties of objects
    inline bool IsObj(AbcObj* pObj) const {return pObj != nullptr;}
    inline bool IsObj(int id) const {return IsObj(GetObj(id));}
    inline bool IsNode(AbcObj* pObj) const {return pObj != nullptr && abc::Abc_ObjIsNode(pObj);}
    inline bool IsNode(int id) const {return IsNode(GetObj(id));}
    inline bool IsObjPi(AbcObj* pObj) const {return abc::Abc_ObjIsPi(pObj);}
    inline bool IsObjPi(int id) const {return IsObjPi(GetObj(id));}
    inline bool IsObjPo(AbcObj* pObj) const {return abc::Abc_ObjIsPo(pObj);}
    inline bool IsObjPo(int id) const {return IsObjPo(GetObj(id));}
    inline bool IsConst(int id) const {if (!IsStrash()) return IsNode(GetObj(id)) && abc::Abc_NodeIsConst(GetObj(id)); else return id == GetConst1IdInStrashNet();}
    inline bool IsConst0(int id) const {assert(!IsStrash()); return IsNode(GetObj(id)) && abc::Abc_NodeIsConst0(GetObj(id));}
    inline bool IsConst1(int id) const {assert(!IsStrash()); return IsNode(GetObj(id)) && abc::Abc_NodeIsConst1(GetObj(id));}
    inline bool IsInv(AbcObj* pObj) const {return IsNode(pObj) && abc::Abc_NodeIsInv(pObj);}
    inline bool IsInv(int id) const {return IsNode(GetObj(id)) && abc::Abc_NodeIsInv(GetObj(id));}
    inline bool IsPoDriver(AbcObj* pObj) const {for (ll i = 0; i < GetFanoutNum(pObj); ++i) { auto fanout = GetFanout(pObj, i); if (IsObjPo(fanout)) return true; } return false;}
    inline bool IsPoDriver(ll id) const {return IsPoDriver(GetObj(id));}
    inline bool IsTheOnlyPoDriver(AbcObj* pObj) const {return GetFanoutNum(pObj) == 1 && IsObjPo(GetFanout(pObj, 0));}
    inline bool IsTheOnlyPoDriver(ll id) const {return IsTheOnlyPoDriver(GetObj(id));}
    inline std::string GetName(AbcObj* pObj) const {return std::string(abc::Abc_ObjName(pObj));}
    inline std::string GetName(int i) const {return GetName(GetObj(i));}
    inline std::string GetPiName(int i) const {return std::string(abc::Abc_ObjName(GetPi(i)));}
    inline std::string GetPoName(int i) const {return std::string(abc::Abc_ObjName(GetPo(i)));}
    inline std::string GetSop(AbcObj* pNode) const {return std::string(static_cast <char *> (pNode->pData));}
    inline std::string GetSop(int id) const {return GetSop(GetObj(id));}
    inline int GetNodeMffcSize(int i) const {return GetNodeMffcSize(GetObj(i));}
    inline void SetNetNotTrav() const {abc::Abc_NtkIncrementTravId(pNtk);}
    inline bool GetObjTrav(AbcObj* pObj) const {return abc::Abc_NodeIsTravIdCurrent(pObj);}
    inline void SetObjTrav(AbcObj* pObj) const {abc::Abc_NodeSetTravIdCurrent(pObj);}
    inline void PrintObj(int id, bool showFunct = false) const {PrintObj(GetObj(id), showFunct);}
    inline void PrintObjBas(int id, std::string && endWith = "") const {PrintObjBas(GetObj(id), std::move(endWith));}

    // get fanins & fanouts
    inline int GetFaninNum(AbcObj* pObj) const {return abc::Abc_ObjFaninNum(pObj);}
    inline int GetFaninNum(int id) const {return GetFaninNum(GetObj(id));}
    inline AbcObj* GetFanin(AbcObj* pObj, int i) const {return abc::Abc_ObjFanin(pObj, i);}
    inline AbcObj* GetFanin(int id, int i) const {return GetFanin(GetObj(id), i);}
    inline int GetFaninId(AbcObj* pObj, int i) const {return GetId(GetFanin(pObj, i));}
    inline int GetFaninId(int id, int i) const {return GetFaninId(GetObj(id), i);}
    inline int GetFaninCompl(int id, int i) const {return abc::Abc_ObjFaninC(GetObj(id), i);}
    inline int GetFanoutNum(AbcObj* pObj) const {return abc::Abc_ObjFanoutNum(pObj);}
    inline int GetFanoutNum(int id) const {return GetFanoutNum(GetObj(id));}
    inline AbcObj* GetFanout(AbcObj* pObj, int i) const {return abc::Abc_ObjFanout(pObj, i);}
    inline AbcObj* GetFanout(int id, int i) const {return GetFanout(GetObj(id), i);}
    inline int GetFanoutId(AbcObj* pObj, int i) const {return GetId(GetFanout(pObj, i));}
    inline int GetFanoutId(int id, int i) const {return GetFanoutId(GetObj(id), i);}
    inline AbcObjVect GetFanoutsThatArePos(AbcObj* pObj) {AbcObjVect pos; for (int i = 0; i < GetFanoutNum(pObj); ++i) { auto fanout = GetFanout(pObj, i); if (IsObjPo(fanout)) pos.emplace_back(fanout); } return pos;}
    inline AbcObjVect GetFanoutsThatArePos(int id) {return GetFanoutsThatArePos(GetObj(id));}

    // get timing information
    inline int GetLev() const {return abc::Abc_NtkLevel(pNtk);}
    inline int GetObjLev(AbcObj* pObj) const {if (IsObjPo(pObj)) { assert(GetFaninNum(pObj) == 1); return abc::Abc_ObjLevel(GetFanin(pObj, 0)) + 1; } else return abc::Abc_ObjLevel(pObj);}
    inline int GetObjLev(int i) const {return abc::Abc_ObjLevel(GetObj(i));}
    inline void SetObjLev(AbcObj* pObj, int lev) {abc::Abc_ObjSetLevel(pObj, lev);}
    inline void SetObjLev(int i, int lev) {SetObjLev(GetObj(i), lev);}
    inline double GetArrTime(AbcObj* pNode) const {if (static_cast <abc::SC_Lib *> (GetAbcFrame()->pLibScl) == nullptr) return static_cast <abc::Abc_Time_t *> (pNode->pNtk->pManTime->vArrs->pArray[pNode->Id])->Rise; else return pNode->dTemp;}
    inline double GetArrTime(int id) const {return GetArrTime(GetObj(id));}
    inline double GetGateDelay(AbcObj* pNode) const {return abc::Mio_GateReadDelayMax(static_cast <abc::Mio_Gate_t *> (pNode->pData));}
    inline double GetGateDelay(int id) const {return GetGateDelay(GetObj(id));}
    inline double GetInvDelay() const {return abc::Mio_LibraryReadDelayInvMax(static_cast <abc::Mio_Library_t *> (pNtk->pManFunc));}
    inline std::string GetGateName(AbcObj* pNode) const {assert(NetMan::GetNetType() == NET_TYPE::GATE); if (IsNode(pNode)) return std::string(abc::Mio_GateReadName(static_cast <abc::Mio_Gate_t *> (pNode->pData))); else return std::string("");}

    // modify network
    inline abc::Abc_Ntk_t* StartSopNet() {pNtk = abc::Abc_NtkAlloc(abc::ABC_NTK_LOGIC, abc::ABC_FUNC_SOP, 1); return pNtk;}
    inline abc::Abc_Ntk_t* StartStrashNet() {pNtk = abc::Abc_NtkAlloc(abc::ABC_NTK_STRASH, abc::ABC_FUNC_AIG, 1); return pNtk;}
    inline void AddFanin(AbcObj* pObj, AbcObj* pFanin) {abc::Abc_ObjAddFanin(pObj, pFanin);}
    inline void Replace(AbcObj* pTS, AbcObj* pSS) {abc::Abc_ObjReplace(pTS, pSS);}
    inline void Replace(int tsId, int ssId) {Replace(GetObj(tsId), GetObj(ssId));}
    inline void TransfFanout(AbcObj* pTS, AbcObj* pSS) {abc::Abc_ObjTransferFanout(pTS, pSS);}
    inline void TransfFanout(int tsId, int ssId) {TransfFanout(GetObj(tsId), GetObj(ssId));}
    inline void DelObj(AbcObj* pObj) {abc::Abc_NtkDeleteObj(pObj);}
    inline void DelObj(int id) {DelObj(GetObj(id));}
    inline void DelObjRec(AbcObj* pObj) {abc::Abc_NtkDeleteObj_rec(pObj, 1);}
    inline void DelObjRec(int id) {DelObjRec(GetObj(id));}
    inline AbcObj* CreateInv(AbcObj* pFanin) {assert(pFanin->pNtk == pNtk); return Abc_NtkCreateNodeInv(pNtk, pFanin);}
    inline int CreateInv(int faninId) {return GetId(CreateInv(GetObj(faninId)));}
    inline AbcObj* CreateBuf(AbcObj* pFanin) {assert(pFanin->pNtk == pNtk); return Abc_NtkCreateNodeBuf(pNtk, pFanin);}
    inline int CreateBuf(int faninId) {return GetId(CreateBuf(GetObj(faninId)));}
    inline AbcObj* CreateAnd(AbcObj* pA, AbcObj* pB) {assert(pA->pNtk == pNtk && pB->pNtk == pNtk); return CreateNode(AbcObjVect({pA, pB}), "11 1\n");}
    inline int CreateAnd(int a, int b) {return GetId(CreateAnd(GetObj(a), GetObj(b)));}
    inline AbcObj* CreateOr(AbcObj* pA, AbcObj* pB) {assert(pA->pNtk == pNtk && pB->pNtk == pNtk); return CreateNode(AbcObjVect({pA, pB}), "00 0\n");}
    inline int CreateOr(int a, int b) {return GetId(CreateOr(GetObj(a), GetObj(b)));}
    inline AbcObj* CreateXor(AbcObj* pA, AbcObj* pB) {assert(pA->pNtk == pNtk && pB->pNtk == pNtk); return CreateNode(AbcObjVect({pA, pB}), "01 1\n10 1\n");}
    inline int CreateXor(int a, int b) {return GetId(CreateXor(GetObj(a), GetObj(b)));}
    inline AbcObj* CreatePo(AbcObj* pFanin, const char* pName) {auto pPo = abc::Abc_NtkCreatePo(pNtk); AddFanin(pPo, pFanin); Abc_ObjAssignName(pPo, const_cast<char*>(pName), nullptr); return pPo;}
    inline AbcObj* CreatePi(const char* pName) {auto pPi = abc::Abc_NtkCreatePi(pNtk); Abc_ObjAssignName(pPi, const_cast<char*>(pName), nullptr); return pPi;}
    inline void RenameNet(const std::string& name) {if (pNtk->pName != nullptr) ABC_FREE(pNtk->pName); pNtk->pName = abc::Extra_UtilStrsav(name.c_str());}
    inline void Rename(AbcObj* pObj, const char* pName) {Abc_ObjAssignName(pObj, const_cast<char*>(pName), nullptr);}
    inline void Rename(int id, const char* pName) {Rename(GetObj(id), pName);}
    inline void RenameWithSuff(AbcObj* pObj, const char* pSuff) {Abc_ObjNameSuffix(pObj, const_cast<char*>(pSuff));}
    inline void RenameWithSuff(int id, const char* pSuff) {RenameWithSuff(GetObj(id), pSuff);}
};


/**
 * @brief Globally start ABC
 * 
 * @return void
 */
static inline void GlobStartAbc() {
    abc::Abc_Start();
    AbcMan abcMan;
    abcMan.LoadAlias();
}


/**
 * @brief Globally stop ABC
 * 
 * @return void
 */
static inline void GlobStopAbc() {
    abc::Abc_Stop();
}


/**
 * @brief Self-defined formatter for ABC object
 */
template <> 
struct fmt::formatter<AbcObj> {
    constexpr auto parse(fmt::format_parse_context& ctx) {
        return ctx.begin();
    }
    template <typename FormatContext>
    auto format(const AbcObj& obj, FormatContext& ctx) {
        auto pObj = const_cast<AbcObj*>(&obj);
        std::string str = fmt::format("{}({})", abc::Abc_ObjName(pObj), obj.Id);
        if (abc::Abc_ObjIsNode(pObj)) {
            abc::Abc_Ntk_t* pNtk = pObj->pNtk;
            if (abc::Abc_NtkIsLogic(pNtk) || abc::Abc_NtkIsNetlist(pNtk)) {
                if (abc::Abc_NodeIsConst0(pObj))
                    str += "(zero)";
                else if (abc::Abc_NodeIsConst1(pObj))
                    str += "(one)";
            }
        }
        return fmt::format_to(ctx.out(), "{}", str);
    }
};
template <> 
struct fmt::formatter<AbcObjVect> {
    constexpr auto parse(fmt::format_parse_context& ctx) {
        return ctx.begin();
    }
    template <typename FormatContext>
    auto format(const AbcObjVect& objs, FormatContext& ctx) {
        std::string str = "[";
        for (const auto& pObj: objs)
            str += fmt::format("{}, ", fmt::format("{}", *pObj));
        return fmt::format_to(ctx.out(), "{}]", str);
    }
};


/**
 * @brief Rename an ABC object
 * 
 * @param pObj    the object to be renamed
 * @param name    the new name
 * @return void
 */
static inline void RenameAbcObj(AbcObj* pObj, const std::string& name) {
    auto pNameMan = pObj->pNtk->pManName;
    auto pEntry = abc::Nm_ManTableLookupId(pNameMan, pObj->Id);
    if (pEntry != nullptr)
        abc::Nm_ManDeleteIdName(pNameMan, pObj->Id);
    abc::Abc_ObjAssignName(pObj, const_cast<char*>(name.c_str()), nullptr);
}


/**
 * @brief Compare PI names of the two networks
 * 
 * @param net0          the first network
 * @param net1          the second network
 * @return bool         true if all the net0's PI names appear in the first part of net1's PI names
 */
static inline bool ComparePi(const NetMan& net0, const NetMan& net1, bool allowExtraPi = false) {
    if (allowExtraPi) {
        if (net0.GetPiNum() > net1.GetPiNum())
            return false;
    }
    else {
        if (net0.GetPiNum() != net1.GetPiNum())
            return false;
    }
    // here, net0.GetPiNum() <= net1.GetPiNum()
    for (int i = 0; i < net0.GetPiNum(); ++i) {
        if (net0.GetPiName(i) != net1.GetPiName(i))
            return false;
    }
    return true;
}


/**
 * @brief Compare PO names of the two networks
 * 
 * @param net0          the first network
 * @param net1          the second network
 * @return bool         true if the PO names are the same
 */
static inline bool ComparePo(const NetMan& net0, const NetMan& net1) {
    if (net0.GetPoNum() != net1.GetPoNum())
        return false;
    for (int i = 0; i < net0.GetPoNum(); ++i) {
        if (net0.GetPoName(i) != net1.GetPoName(i))
            return false;
    }
    return true;
}

