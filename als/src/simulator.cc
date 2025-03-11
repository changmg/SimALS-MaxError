/**
 * @file simulator.cc
 * @author Chang Meng (chang.meng@epfl.ch)
 * @brief Circuit simulator
 * @date 2025-01-25
 * 
 */
#include "simulator.h"


using namespace abc;
using namespace std;
using namespace boost;
using namespace random;


/**
 * @brief Constructor of the simulator
 * 
 * @param net_man       the network to be simulated
 * @param _seed         the seed for random number generation
 * @param n_frame       the number of simulation frames
 * @param distr_type    the distribution type
 */
Simulator::Simulator(const NetMan& net_man, unsigned _seed, int n_frame, DISTR_TYPE distr_type): NetMan(net_man.GetNet(), false), seed(_seed), nFrame(n_frame), distrType(distr_type) {
    auto type = NetMan::GetNetType();
    assert(type == NET_TYPE::AIG || type == NET_TYPE::GATE || type == NET_TYPE::SOP || type == NET_TYPE::STRASH);
    if (distrType == DISTR_TYPE::ENUM) {
        int nPi = net_man.GetPiNum();
        assert(nPi < 30);
        nFrame = 1 << nPi;
    }
    dat.resize(NetMan::GetIdMaxPlus1(), BitVect(nFrame, 0));
}


/**
 * @brief Initialize the constant nodes
 * 
 * @retval dat  the simulation patterns for the constant nodes
 * @return void
 */
void Simulator::InitConstNodes() {
    if (!NetMan::IsStrash()) {
        for (int i = 0; i < NetMan::GetIdMaxPlus1(); ++i) {
            if (NetMan::IsConst0(i))
                dat[i].reset();
            else if (NetMan::IsConst1(i))
                dat[i].set();
        }
    }
    else
        dat[NetMan::GetConst1IdInStrashNet()].set();
}


/**
 * @brief Generate random input patterns using the uniform distribution
 * 
 * @retval dat  the simulation patterns for each node
 * @return void
 */
void Simulator::GenInpUnif() {
    // generate random input patterns for the PIs
    uniform_int<> unif01(0, 1);
    random::mt19937 eng(seed);
    variate_generator < random::mt19937, uniform_int<> > rand01(eng, unif01);
    assert(static_cast<int>(dat.size()) == NetMan::GetIdMaxPlus1());
    for (int i = 0; i < NetMan::GetPiNum(); ++i) {
        auto piId = NetMan::GetPiId(i);
        dat[piId].reset();
        for (int j = 0; j < nFrame; ++j) {
            if (rand01())
                dat[piId].set(j);
        }
    }
    // init the constant nodes
    InitConstNodes();
}


/**
 * @brief Generate random input patterns using the uniform distribution
 * @brief Fast version
 * 
 * @retval dat  the simulation patterns for each node
 * @return void
 */
void Simulator::GenInpUnifFast() {
    // if the number of frames is not a multiple of 64, use the normal version
    const int unitLength = 64;
    if ((nFrame & (unitLength - 1)) != 0) {
        GenInpUnif();
        return;
    }

    // prepare
    random::uniform_int_distribution <ull> unifUll;
    random::mt19937 eng(seed);
    variate_generator <random::mt19937, random::uniform_int_distribution <ull> > randUll(eng, unifUll);
    int nUnit = nFrame / unitLength;

    // generate random input patterns
    assert(static_cast<int>(dat.size()) == NetMan::GetIdMaxPlus1());
    for (int i = 0; i < NetMan::GetPiNum(); ++i) {
        auto piId = NetMan::GetPiId(i);
        dat[piId].resize(0);
        for (int j = 0; j < nUnit; ++j) {
            ull numb = randUll();
            dat[piId].append(numb);
        }
    }

    // init the constant nodes
    InitConstNodes();
}


/**
 * @brief Generate random input patterns using the enumeration method
 * 
 * @retval dat  the simulation patterns for each node
 * @return void
 */
void Simulator::GenInpEnum() {
    // Generate the input patterns for the PIs
    assert(GetPiNum() < 30);
    assert(1ll << GetPiNum() == nFrame);
    assert(static_cast<int>(dat.size()) == NetMan::GetIdMaxPlus1());
    for (int i = 0; i < GetPiNum(); ++i) {
        bool phase = 1;
        auto piId = GetPiId(i);
        dat[piId].reset();
        for (int j = 0; j < nFrame; ++j) {
            if (j % (1 << i) == 0)
                phase = !phase;
            if (phase)
                dat[piId].set(j);
        }
    }

    // init the constant nodes
    InitConstNodes();
}


/**
 * @brief Replace the iPatt-th simulation pattern of node objIds with the given values
 * 
 * @param iPatt    the pattern index
 * @param piVals   the PI values
 * @retval dat     the simulation patterns of the PI nodes
 * @return void
 */
void Simulator::ReplInp(int iPatt, const IntVect& piVals) {
    assert(NetMan::GetPiNum() == static_cast<int>(piVals.size()));
    assert(static_cast<int>(dat.size()) == NetMan::GetIdMaxPlus1());
    for (int iPi = 0; iPi < NetMan::GetPiNum(); ++iPi) {
        auto& piDat = dat[NetMan::GetPiId(iPi)];
        assert(iPatt < static_cast<int>(piDat.size()));
        assert(piVals[iPi] == 0 || piVals[iPi] == 1);
        dat[NetMan::GetPiId(iPi)][iPatt] = piVals[iPi];
    }
}


/**
 * @brief Append the given values to the simulation patterns of the PI nodes
 * 
 * @param piVals  the PI values
 * @retval dat     the simulation patterns of the PI nodes
 * @return void
 */
void Simulator::AppendInp(const IntVect& piVals) {
    assert(NetMan::GetPiNum() == static_cast<int>(piVals.size()));
    assert(static_cast<int>(dat.size()) == NetMan::GetIdMaxPlus1());
    ++nFrame;
    for (int iPi = 0; iPi < NetMan::GetPiNum(); ++iPi) {
        auto& piDat = dat[NetMan::GetPiId(iPi)];
        assert(piVals[iPi] == 0 || piVals[iPi] == 1);
        piDat.push_back(piVals[iPi]);
        assert(static_cast<int>(piDat.size()) == nFrame);
    }
}


/**
 * @brief Initialize the PI patterns using those of the other simulator
 * 
 * @param othSmlt    the other simulator
 * @retval dat       the simulation patterns of the PI objects
 * @return void
 */
void Simulator::GenInpFromOthSmlt(const Simulator& othSmlt) {
    // check
    assert(this->IsPIOSame(othSmlt));
    assert(static_cast<int>(dat.size()) == NetMan::GetIdMaxPlus1());
    // copy the PI patterns
    for (int i = 0; i < NetMan::GetPiNum(); ++i) {
        auto piId = NetMan::GetPiId(i);
        auto& piDat = othSmlt.GetDat(piId);
        assert(static_cast<int>(piDat.size()) == nFrame);
        dat[piId] = othSmlt.GetDat(piId);
    }
    // // resize the PI patterns
    // if (nFrame != othSmlt.GetFrameNumb()) {
    //     for (int i = 0; i < NetMan::GetPiNum(); ++i)
    //         dat[NetMan::GetPiId(i)].resize(nFrame);
    // }
    // init the constant nodes
    InitConstNodes();
}


/**
 * @brief Generate the input patterns for the circuit using the given bit vectors
 * 
 * @param piPatts[i]      the PI patterns (a bit vector) for the i-th PI node
 * @retval dat            the simulation patterns for the PI nodes
 * @return void 
 */
void Simulator::GenInpFromBitVects(const std::vector<BitVect>& piPatts) {
    // check
    int piPattSize = static_cast<int>(piPatts.size());
    assert(piPattSize <= NetMan::GetPiNum() && piPattSize > 0); // the number of PI patterns should be no more than the number of PIs
    assert(static_cast<int>(dat.size()) == NetMan::GetIdMaxPlus1());
    // update the frame number
    if (nFrame != static_cast<int>(piPatts[0].size()))
        nFrame = static_cast<int>(piPatts[0].size());
    // copy the PI patterns
    for (int i = 0; i < piPattSize; ++i) {
        assert(nFrame == static_cast<int>(piPatts[i].size()));
        dat[NetMan::GetPiId(i)] = piPatts[i];
    }
    // init the constant nodes
    InitConstNodes();
}


/**
 * @brief Logic simulation for the network NetMan::GetNet()
 * @retval dat  the simulation patterns for each node
 * @return void
 */
void Simulator::UpdNodeAndPoPatts() {
    auto type = GetNetType();
    auto nodes = CalcTopoOrd(false);
    for (const auto& pObj: nodes) {
        if (type == NET_TYPE::AIG)
            UpdAigNode(pObj);
        else if (type == NET_TYPE::SOP)
            UpdSopNode(pObj);
        else if (type == NET_TYPE::GATE)
            UpdGateNode(pObj);
        else if (type == NET_TYPE::STRASH)
            UpdStrashNode(pObj);
        else
            assert(0);
    }
    for (int i = 0; i < GetPoNum(); ++i) {
        auto pPo = GetPo(i);
        auto drivId = GetFaninId(pPo, 0);
        assert(!Abc_ObjIsComplement(pPo));
        if (type == NET_TYPE::AIG || type == NET_TYPE::GATE || type == NET_TYPE::SOP)
            dat[GetId(pPo)] = dat[drivId];
        else if (type == NET_TYPE::STRASH)
            dat[GetId(pPo)] = Abc_ObjFaninC0(pPo)? ~dat[drivId]: dat[drivId];
        else
            assert(0);
    }
}


/**
 * @brief Simulate a node with the given sum-of-products (SOP) function
 * 
 * @param faninIds      the fanin ids of the node
 * @param sop           the SOP of the node
 * @param res           the simulation result
 * @retval res          the simulation result
 * @return void
 */
void Simulator::SimSop(const IntVect& faninIds, const std::string& sop, BitVect& res) {
    if (sop == " 0\n") {
        res.reset();
        return;
    }
    if (sop == " 1\n") {
        res.set();
        return;
    }
    char* pSop = const_cast <char*> (sop.c_str());
    int nVars = Abc_SopGetVarNum(pSop);
    assert(nVars == static_cast<int>(faninIds.size()));

    BitVect product(nFrame, 0);
    for (char* pCube = pSop; *pCube; pCube += nVars + 3) {
        bool isFirst = true;
        for (int i = 0; pCube[i] != ' '; i++) {
            int faninId = faninIds[i];
            switch (pCube[i]) {
                case '-':
                    continue;
                    break;
                case '0':
                    if (isFirst) {
                        isFirst = false;
                        product = ~dat[faninId];
                    }
                    else
                        product &= ~dat[faninId];
                    break;
                case '1':
                    if (isFirst) {
                        isFirst = false;
                        product = dat[faninId];
                    }
                    else
                        product &= dat[faninId];
                    break;
                default:
                    assert(0);
            }
        }
        if (isFirst) {
            isFirst = false;
            product.set();
        }
        assert(!isFirst);
        if (pCube == pSop)
            res = product;
        else
            res |= product;
    }

    // complement
    if (Abc_SopIsComplement(pSop))
        res.flip();
}


void Simulator::UpdAigNode(AbcObj* pObj) {
    assert(Abc_ObjIsNode(pObj));
    auto pNtk = NetMan::GetNet();
    auto pMan = static_cast <Hop_Man_t *> (pNtk->pManFunc);
    auto pRoot = static_cast <Hop_Obj_t *> (pObj->pData);
    auto pRootR = Hop_Regular(pRoot);

    // skip constant node
    if (Hop_ObjIsConst1(pRootR))
        return;

    // get topological order of subnetwork in aig
    Vec_Ptr_t * vHopNodes = Hop_ManDfsNode(pMan, pRootR);

    // init internal hop nodes
    int maxHopId = -1;
    int i = 0;
    Hop_Obj_t * pHopObj = nullptr;
    Vec_PtrForEachEntry(Hop_Obj_t *, vHopNodes, pHopObj, i)
        maxHopId = max(maxHopId, pHopObj->Id);
    Vec_PtrForEachEntry( Hop_Obj_t *, pMan->vPis, pHopObj, i )
        maxHopId = max(maxHopId, pHopObj->Id);
    vector < BitVect > interData(maxHopId + 1, BitVect (nFrame, 0));
    unordered_map <int, BitVect *> hop2Data;
    AbcObj* pFanin = nullptr;
    Abc_ObjForEachFanin(pObj, pFanin, i)
        hop2Data[Hop_ManPi(pMan, i)->Id] = &dat[pFanin->Id];

    // special case for inverter or buffer
    if (pRootR->Type == AIG_PI) {
        pFanin = Abc_ObjFanin0(pObj);
        dat[pObj->Id] = dat[pFanin->Id];
    }

    // simulate
    Vec_PtrForEachEntry(Hop_Obj_t *, vHopNodes, pHopObj, i) {
        assert(Hop_ObjIsAnd(pHopObj));
        auto pHopFanin0 = Hop_ObjFanin0(pHopObj);
        auto pHopFanin1 = Hop_ObjFanin1(pHopObj);
        assert(!Hop_ObjIsConst1(pHopFanin0));
        assert(!Hop_ObjIsConst1(pHopFanin1));
        BitVect& data0 = Hop_ObjIsPi(pHopFanin0) ? *hop2Data[pHopFanin0->Id] : interData[pHopFanin0->Id];
        BitVect& data1 = Hop_ObjIsPi(pHopFanin1) ? *hop2Data[pHopFanin1->Id] : interData[pHopFanin1->Id];
        BitVect& out = (pHopObj == pRootR) ? dat[pObj->Id] : interData[pHopObj->Id];
        bool isFanin0C = Hop_ObjFaninC0(pHopObj);
        bool isFanin1C = Hop_ObjFaninC1(pHopObj);
        if (!isFanin0C && !isFanin1C)
            out = data0 & data1;
        else if (!isFanin0C && isFanin1C)
            out = data0 & ~data1;
        else if (isFanin0C && !isFanin1C)
            out = ~data0 & data1;
        else if (isFanin0C && isFanin1C)
            out = ~(data0 | data1);
    }

    // complement
    if (Hop_IsComplement(pRoot))
        dat[pObj->Id].flip();

    // recycle memory
    Vec_PtrFree(vHopNodes); 
}


void Simulator::UpdSopNode(AbcObj* pObj) {
    assert(Abc_ObjIsNode(pObj));
    // skip constant node
    if (Abc_NodeIsConst(pObj))
        return;
    // update sop
    char* pSop = static_cast <char*> (pObj->pData);
    UpdSop(pObj, pSop);
}


void Simulator::UpdGateNode(AbcObj* pObj) {
    assert(Abc_ObjIsNode(pObj));
    // skip constant node
    if (Abc_NodeIsConst(pObj))
        return;
    // update sop
    char* pSop = static_cast <char*> ((static_cast <Mio_Gate_t *> (pObj->pData))->pSop);
    UpdSop(pObj, pSop);
}


void Simulator::UpdStrashNode(AbcObj* pObj) {
    assert(Abc_ObjIsNode(pObj));
    // skip constant node
    auto pConst1 = Abc_AigConst1(NetMan::GetNet());
    if (Abc_ObjRegular(pObj) == pConst1) {
        assert(0);
        return;
    }
    // update node
    assert(!Abc_ObjIsComplement(pObj));
    auto pFanin0 = Abc_ObjFanin0(pObj), pFanin1 = Abc_ObjFanin1(pObj);
    int compl0 = Abc_ObjFaninC0(pObj), compl1 = Abc_ObjFaninC1(pObj);
    dynamic_bitset<ull> c0(nFrame, 0); if (compl0) c0.set();
    dynamic_bitset<ull> c1(nFrame, 0); if (compl1) c1.set();
    dat[pObj->Id] = (dat[pFanin0->Id] ^ c0) & (dat[pFanin1->Id] ^ c1);
}


void Simulator::UpdSop(AbcObj* pObj, char* pSop) {
    int nVars = Abc_SopGetVarNum(pSop);
    BitVect product(nFrame, 0);
    for (char* pCube = pSop; *pCube; pCube += nVars + 3) {
        bool isFirst = true;
        for (int i = 0; pCube[i] != ' '; i++) {
            AbcObj* pFanin = Abc_ObjFanin(pObj, i);
            switch (pCube[i]) {
                case '-':
                    continue;
                    break;
                case '0':
                    if (isFirst) {
                        isFirst = false;
                        product = ~dat[pFanin->Id];
                    }
                    else
                        product &= ~dat[pFanin->Id];
                    break;
                case '1':
                    if (isFirst) {
                        isFirst = false;
                        product = dat[pFanin->Id];
                    }
                    else
                        product &= dat[pFanin->Id];
                    break;
                default:
                    assert(0);
            }
        }
        if (isFirst) {
            isFirst = false;
            product.set();
        }
        assert(!isFirst);
        if (pCube == pSop)
            dat[pObj->Id] = product;
        else
            dat[pObj->Id] |= product;
    }

    // complement
    if (Abc_SopIsComplement(pSop))
        dat[pObj->Id].flip();
}


/**
 * @brief Get the circuit input encoded between lsb and msb for the iPatt-th pattern
 * 
 * @param iPatt      the pattern index
 * @param lsb        the least significant bit
 * @param msb        the most significant bit
 * @return BigInt    the input value
 */
BigInt Simulator::GetInput(int iPatt, int lsb, int msb) const {
    assert(lsb >= 0 && msb < NetMan::GetPiNum());
    assert(iPatt < nFrame);
    assert(lsb <= msb && msb - lsb < 512);
    BigInt ret(0);
    for (int k = msb; k >= lsb; --k) {
        ret <<= 1;
        if (dat[NetMan::GetPiId(k)][iPatt])
            ++ret;
    }
    return ret;
}


/**
 * @brief Get the PI vector for the iPatt-th pattern
 * 
 * @param iPatt   the pattern index
 * @param piVals  the PI values
 * @retval piVals the PI values
 * @return void
 */
void Simulator::GetInpVect(int iPatt, IntVect& piVals) const {
    piVals.resize(NetMan::GetPiNum());
    for (int i = 0; i < NetMan::GetPiNum(); ++i)
        piVals[i] = dat[NetMan::GetPiId(i)][iPatt];
}


/**
 * @brief Print the PI pattern of the circuit for the iPatt-th pattern
 * 
 * @param iPatt  the pattern index
 */
void Simulator::PrintInpStream(int iPatt) const {
    assert(iPatt < nFrame);
    for (int k = GetPiNum() - 1; k >= 0; --k)
        // cout << dat[GetPiId(k)][iPatt];
        fmt::print("{}", dat[GetPiId(k)][iPatt]? 1: 0);
    // cout << endl;
    fmt::print("\n");
}


/**
 * @brief Get the output of the circuit for the iPatt-th pattern
 * @brief The output is represented as a binary number
 * @brief The number of POs should be less than 500
 * @brief Assume the output value is an unsigned integer
 * 
 * @param  iPatt   the pattern index
 * @param  ret     the output value
 * @retval ret     the output value
 * @return void
 */
void Simulator::GetOutput(int iPatt, BigInt& ret) const {
    int msb = NetMan::GetPoNum() - 1;
    assert(iPatt < nFrame);
    assert(msb < 500);
    ret = 0;
    for (int k = msb; k >= 0; --k) {
        ret <<= 1;
        ret |= dat[NetMan::GetPoId(k)][iPatt];
    }
}


/**
 * @brief Get the output of the circuit for the iPatt-th pattern (fast version)
 * @brief The output is represented as a binary number
 * @brief The number of POs should be less than 63
 * @brief Assume the output value is an unsigned integer
 * 
 * @param iPatt  the pattern index
 * @return ll    the output value
 */
ll Simulator::GetOutputFast(int iPatt) const {
    int msb = NetMan::GetPoNum() - 1;
    assert(iPatt < nFrame);
    assert(msb < 63);
    ll ret = 0;
    for (int k = msb; k >= 0; --k) {
        ret <<= 1;
        ret |= dat[NetMan::GetPoId(k)][iPatt];
    }
    return ret;
}


/**
 * @brief Get the Temprorary output of the circuit for the iPatt-th pattern (fast version)
 * @brief The output is represented as a binary number
 * @brief The number of POs should be less than 64
 * 
 * @param iPatt  the pattern index
 * @return ll    the output value
 */
ll Simulator::GetTempOutputFast(int iPatt) const {
    int msb = NetMan::GetPoNum() - 1;
    assert(iPatt < nFrame);
    assert(msb < 64);
    int ret = 0;
    for (int k = msb; k >= 0; --k) {
        ret <<= 1;
        ret |= tempDat[NetMan::GetPoId(k)][iPatt];
    }
    return ret;
}


/**
 * @brief Print the PO pattern of the circuit for the iPatt-th pattern
 * 
 * @param iPatt  the pattern index
 */
void Simulator::PrintOutpStream(int iPatt) const {
    assert(iPatt < nFrame);
    for (int k = GetPoNum() - 1; k >= 0; --k)
        cout << dat[GetPoId(k)][iPatt];
    cout << endl;
}


double Simulator::GetSignalProb(int objId) const {
    assert(objId < NetMan::GetIdMaxPlus1());
    return dat[objId].count() / static_cast <double> (nFrame);
}


void Simulator::PrintSignalProb() const {
    for (int i = 0; i < NetMan::GetPoNum(); ++i) {
        cout << NetMan::GetName(NetMan::GetPo(i)) << " " << GetSignalProb(NetMan::GetPoId(i)) << endl;
    }
}


double Simulator::GetErrRate(const Simulator& oth_smlt, bool isCheck) const {
    if (isCheck)
        assert(NetMan::IsPIOSame(static_cast<const NetMan&>(oth_smlt)));
    BitVect temp(nFrame, 0);
    for (int i = 0; i < NetMan::GetPoNum(); ++i)
        temp |= (this->dat[NetMan::GetPoId(i)] ^ oth_smlt.dat[oth_smlt.GetPoId(i)]);
    return temp.count() / static_cast <double> (nFrame);
}


double Simulator::GetMeanErrDist(const Simulator& oth_smlt, bool isCheck) const {
    if (isCheck) {
        assert(NetMan::IsPIOSame(static_cast<const NetMan&>(oth_smlt)));
        assert(GetPoNum() < 63);
    }
    const int warnPoint = numeric_limits <int>::max() * 0.8;
    ll sed = 0;
    for (int i = 0; i < nFrame; ++i) {
        ll accOut = GetOutputFast(i);
        ll appOut = oth_smlt.GetOutputFast(i);
        sed += abs(accOut - appOut);
        if (sed >= warnPoint)
            assert(0);
    }
    // cout << "sed = " << sed << ", nFrame = " << nFrame << ", med = " << sed / static_cast <double> (nFrame) << endl;
    return sed / static_cast <double> (nFrame);
}


/**
 * @brief Get the maximum error distance the two simulators, this and oth_smlt (fast version)
 * 
 * @param  oth_smlt  the other simulator
 * @param  isCheck   whether to check the PI/PO information
 * @return ll        the maximum error distance
 */
ll Simulator::GetMaxErrDistFast(const Simulator& oth_smlt, bool isCheck) const {
    if (isCheck) {
        assert(NetMan::IsPIOSame(static_cast<const NetMan&>(oth_smlt)));
        assert(GetPoNum() < 63);
    }
    ll maxEd = 0;
    for (int i = 0; i < nFrame; ++i) {
        ll accOut = GetOutputFast(i);
        ll appOut = oth_smlt.GetOutputFast(i);
        maxEd = std::max(maxEd, abs(accOut - appOut));
    }
    return maxEd;
}


/**
 * @brief Get the maximum error distance the two simulators, this and oth_smlt
 * 
 * @param  oth_smlt        the other simulator
 * @param  isCheck         whether to check the PI/PO information
 * @param  maxErrLowBound  the maximum error distance (lower bound)
 * @retval maxErrLowBound  the maximum error distance (lower bound)
 * @return void
 */
void Simulator::GetMaxErrDist(const Simulator& oth_smlt, bool isCheck, BigInt& maxErrLowBound) const {
    if (isCheck) {
        assert(NetMan::IsPIOSame(static_cast<const NetMan&>(oth_smlt)));
        assert(GetPoNum() < 500);
    }
    maxErrLowBound = 0;
    BigInt accOut(0), appOut(0);
    for (int i = 0; i < nFrame; ++i) {
        GetOutput(i, accOut);
        oth_smlt.GetOutput(i, appOut);
        maxErrLowBound = std::max(maxErrLowBound, abs(accOut - appOut));
    }
}


/**
 * @brief Calculate the Boolean difference of the POs with respect to the node topoNodes[i]
 * 
 * @param  topoNodes             the nodes in topological order
 * @param  idx                   the index of the node to be flipped 
 * @param  bdPosWrtNode          the Boolean difference of the POs with respect to the node topoNodes[i]
 * @retval bdPosWrtNode[poIdx]   Boolean difference of the poIdx-th PO w.r.t. the node topoNodes[idx]
 * @return void
 */
void Simulator::CalcBoolDiff(const AbcObjVect& topoNodes, int idx, std::vector<BitVect>& bdPosWrtNode) {
    // check
    assert(GetNetType() == NET_TYPE::SOP);
    assert(idx < static_cast<int>(topoNodes.size()));

    // initizalize tempDat
    if (tempDat.size() != dat.size())
        tempDat.resize(dat.size(), BitVect(nFrame, 0));
    
    // prepare traversal mark
    auto pTarg = topoNodes[idx];
    SetNetNotTrav();   // mark all nodes as untraversed
    SetObjTrav(pTarg); // mark the i-th node as traversed

    // flip the node
    tempDat[pTarg->Id] = ~dat[pTarg->Id];
    // fmt::print("flip target {}: {}\n", *pTarg, tempDat[pTarg->Id]);

    // simulate from the i-th node
    for (int k = idx + 1; k < static_cast<int>(topoNodes.size()); ++k) {
        auto pObj = topoNodes[k];
        // get the traversing status of the fanins
        bool isOneFaninTrav = false;
        for (int iFanin = 0; iFanin < GetFaninNum(pObj); ++iFanin) {
            auto pFanin = GetFanin(pObj, iFanin);
            if (GetObjTrav(pFanin)) {
                isOneFaninTrav = true;
                break;
            }
        }
        // only if one of the fanins are traversed, then we need to update the node
        // otherwise, it means the node does not depend on the flipped node
        if (isOneFaninTrav) {
            UpdSopNodeForBoolDiff(pObj);
            // fmt::print("update node {}: {}\n", *pObj, tempDat[pObj->Id]);
        }
    }
    for (int i = 0; i < GetPoNum(); ++i) {
        UpdSopNodeForBoolDiff(GetPo(i));
        // fmt::print("update node {}: {}\n", *GetPo(i), tempDat[GetPoId(i)]);
    }

    // get Boolean difference
    bdPosWrtNode.resize(GetPoNum());
    for (int i = 0; i < GetPoNum(); ++i) {
        int poId = GetPoId(i);
        bdPosWrtNode[i] = dat[poId] ^ tempDat[poId];
        // fmt::print("po: {}, bdPosWrtNode[{}]: {}, dat[poId]: {}, tempDat[poId]: {}\n", *GetPo(i), i, bdPosWrtNode[i], dat[poId], tempDat[poId]);
    }
}


void Simulator::CalcLocBoolDiff(AbcObj* pObj, list <AbcObj*>& disjCut, vector <AbcObj*>& cutNtk, vector < BitVect >& bdCut2Node) {
    assert(pObj->pNtk == GetNet());
    if (tempDat.size() != dat.size())
        tempDat.resize(dat.size(), BitVect(nFrame, 0));
    // flip the node
    tempDat[pObj->Id] = ~dat[pObj->Id];
    // simulate
    Abc_NtkIncrementTravId(GetNet());
    Abc_NodeSetTravIdCurrent(pObj);
    for (auto& pInner: cutNtk)
        Abc_NodeSetTravIdCurrent(pInner);
    auto type = NetMan::GetNetType();
    if (type == NET_TYPE::AIG)
        assert(0);
        // UpdAigNodeForBoolDiff(pObj);
    else if (type == NET_TYPE::SOP) {
        for (auto& pInner: cutNtk)
            UpdSopNodeForBoolDiff(pInner);
    }
    else if (type == NET_TYPE::GATE) {
        for (auto& pInner: cutNtk)
            UpdGateNodeForBoolDiff(pInner);
    }
    else
        assert(0);
    // get boolean difference from the node to its disjoint cuts
    bdCut2Node.resize(disjCut.size(), BitVect(nFrame, 0));
    int i = 0;
    for (auto& pCut: disjCut) {
        bdCut2Node[i] = dat[pCut->Id] ^ tempDat[pCut->Id];
        ++i;
    }
}


/**
 * @brief Simulate the node pObj with its sum-of-products (SOP) function
 * @brief The simulation result is stored in tempDat
 * @brief Auxiliary function for computing Boolean difference
 * 
 * @param pObj   the node to be simulated
 */
void Simulator::UpdSopNodeForBoolDiff(AbcObj* pObj) {
    assert(!Abc_ObjIsPi(pObj));
    assert(!Abc_NodeIsConst(pObj));
    if (Abc_ObjIsPo(pObj)) {
        AbcObj* pDriver = Abc_ObjFanin0(pObj);
        assert(!Abc_ObjIsComplement(pObj));
        if (Abc_NodeIsTravIdCurrent(pDriver))
            tempDat[pObj->Id] = tempDat[pDriver->Id];
        else
            tempDat[pObj->Id] = dat[pDriver->Id];
        return;
    }
    UpdSopForBoolDiff(pObj, static_cast<char*>(pObj->pData));
}


void Simulator::UpdGateNodeForBoolDiff(AbcObj* pObj) {
    assert(!Abc_ObjIsPi(pObj));
    assert(!Abc_NodeIsConst(pObj));
    if (Abc_ObjIsPo(pObj)) {
        AbcObj* pDriver = Abc_ObjFanin0(pObj);
        if (Abc_NodeIsTravIdCurrent(pDriver))
            tempDat[pObj->Id] = tempDat[pDriver->Id];
        else
            tempDat[pObj->Id] = dat[pDriver->Id];
        return;
    }
    // update sop
    char* pSop = static_cast <char*> ((static_cast <Mio_Gate_t *> (pObj->pData))->pSop);
    UpdSopForBoolDiff(pObj, pSop);
}


void Simulator::UpdSopForBoolDiff(AbcObj* pObj, char* pSop) {
    int nVars = Abc_SopGetVarNum(pSop);
    BitVect product(nFrame, 0);
    for (char* pCube = pSop; *pCube; pCube += nVars + 3) {
        bool isFirst = true;
        for (int i = 0; pCube[i] != ' '; i++) {
            AbcObj* pFanin = Abc_ObjFanin(pObj, i);
            BitVect& datFi = Abc_NodeIsTravIdCurrent(pFanin)? tempDat[pFanin->Id]: dat[pFanin->Id];
            // fmt::print("pObj: {}, pFanin: {}, datFi: {}\n", *pObj, *pFanin, datFi);
            switch (pCube[i]) {
                case '-':
                    continue;
                    break;
                case '0':
                    if (isFirst) {
                        isFirst = false;
                        product = ~datFi;
                    }
                    else
                        product &= ~datFi;
                    break;
                case '1':
                    if (isFirst) {
                        isFirst = false;
                        product = datFi;
                    }
                    else
                        product &= datFi;
                    break;
                default:
                    assert(0);
            }
        }
        if (isFirst) {
            isFirst = false;
            product.set();
        }
        assert(!isFirst);
        if (pCube == pSop) {
            tempDat[pObj->Id] = product;
        }
        else
            tempDat[pObj->Id] |= product;
    }

    // complement
    if (Abc_SopIsComplement(pSop))
        tempDat[pObj->Id].flip();

    // mark the node as traversed
    Abc_NodeSetTravIdCurrent(pObj);
}


/**
 * @brief Print the simulation patterns of the circuit
 *
 * @return void 
 */
void Simulator::PrintDat() const {
    fmt::print("{}Simulation patterns{}\n", HALF_DASH_LINE, HALF_DASH_LINE);
    for (int i = 0; i < GetIdMaxPlus1(); ++i) {
        if (IsObj(i))
            fmt::print("{}:{}\n", *GetObj(i), dat[i]);
    }
    fmt::print("{}\n", DASH_LINE);
}