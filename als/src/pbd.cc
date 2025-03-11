/**
 * @file pbd.cc
 * @author Chang Meng (chang.meng@epfl.ch)
 * @brief Partial Boolean difference (PBD) computation
 * @date 2025-02-21
 * 
 */


#include "pbd.h"
#include "simulator.h"


using abc::Abc_Aig_t;
using abc::Abc_Obj_t;
using abc::Abc_Ntk_t;
using boost::dynamic_bitset;
using std::bitset;
using std::cout;
using std::endl;
using std::list;
using std::pair;
using std::string;
using std::thread;
using std::unordered_map;
using std::vector;


namespace mecals_v1 {


void PBDMan::BuildMit(const NetMan& accNet, NetMan& appNet, NetMan& mitNet) {
    // init
    assert(accNet.IsStrash() && appNet.IsStrash() && mitNet.IsStrash());
    assert(accNet.IsPIOSame(appNet));
    net.StartStrashNet();
    // copy accNet
    auto pNet = net.GetNet();
    auto pAccNet = accNet.GetNet();
    AbcObj* pObj = nullptr;
    int i = 0;
    Abc_NtkCleanCopy(pAccNet);
    Abc_AigConst1(pAccNet)->pCopy = Abc_AigConst1(pNet);
    Abc_NtkForEachPi(pAccNet, pObj, i) {
        Abc_NtkDupObj(pNet, pObj, 0);
        RenameAbcObj(pObj->pCopy, string(Abc_ObjName(pObj)) + "_acc");
    }
    Abc_AigForEachAnd(pAccNet, pObj, i) {
        pObj->pCopy = Abc_AigAnd((Abc_Aig_t *)pNet->pManFunc, Abc_ObjChild0Copy(pObj), Abc_ObjChild1Copy(pObj));
        RenameAbcObj(pObj->pCopy, string(Abc_ObjName(pObj)) + "_acc");
    }
    // copy appNet
    auto pAppNet = appNet.GetNet();
    Abc_NtkCleanCopy(pAppNet);
    Abc_AigConst1(pAppNet)->pCopy = Abc_AigConst1(pNet);
    Abc_NtkForEachPi(pAppNet, pObj, i) {
        Abc_NtkDupObj(pNet, pObj, 0);
        RenameAbcObj(pObj->pCopy, string(Abc_ObjName(pObj)) + "_app");
    }
    Abc_AigForEachAnd(pAppNet, pObj, i) {
        pObj->pCopy = Abc_AigAnd((Abc_Aig_t *)pNet->pManFunc, Abc_ObjChild0Copy(pObj), Abc_ObjChild1Copy(pObj));
        RenameAbcObj(pObj->pCopy, string(Abc_ObjName(pObj)) + "_app");
        auto pPo = Abc_NtkCreatePo(pNet);
        assert(!Abc_ObjIsComplement(pObj->pCopy));
        Abc_ObjAddFanin(pPo, pObj->pCopy);
        RenameAbcObj(pPo, string(Abc_ObjName(pObj)) + "_app");
    }
    Abc_NtkForEachPo(pAppNet, pObj, i) {
        Abc_NtkDupObj(pNet, pObj, 0);
        Abc_ObjAddFanin(pObj->pCopy, Abc_ObjChild0Copy(pObj));
        RenameAbcObj(pObj->pCopy, string(Abc_ObjName(pObj)) + "_po");
    }
    // copy miter
    auto pMitNet = mitNet.GetNet();
    int nPo = accNet.GetPoNum();
    assert(appNet.GetPoNum() == nPo && mitNet.GetPiNum() == nPo * 2);
    assert(mitNet.GetPoNum() == 1);
    Abc_NtkCleanCopy(pMitNet);
    Abc_AigConst1(pMitNet)->pCopy = Abc_AigConst1(pNet);
    Abc_NtkForEachPo(pAccNet, pObj, i)
        Abc_NtkPi(pMitNet, i)->pCopy = Abc_ObjChild0Copy(pObj);
    Abc_NtkForEachPo(pAppNet, pObj, i)
        Abc_NtkPi(pMitNet, i + nPo)->pCopy = Abc_ObjChild0Copy(pObj);
    Abc_AigForEachAnd(pMitNet, pObj, i) {
        pObj->pCopy = Abc_AigAnd((Abc_Aig_t *)pNet->pManFunc, Abc_ObjChild0Copy(pObj), Abc_ObjChild1Copy(pObj));
        // if (!Abc_ObjIsComplement(pObj->pCopy))
        //     RenameAbcObj(pObj->pCopy, string(Abc_ObjName(pObj)) + "_mit");
    }
    Abc_NtkForEachPo(pMitNet, pObj, i)
        Abc_NtkDupObj(pNet, pObj, 1);
    Abc_NtkForEachPo(pMitNet, pObj, i)
        Abc_ObjAddFanin(pObj->pCopy, Abc_ObjChild0Copy(pObj));
    // cout << "miter net" << endl;
    // net.PrintStat();
    // net.Print(true); cout << endl;
}


void PBDMan::BuildPBD(double exactPBDPerc) {
    assert(net.IsStrash());

    // update fanout information
    auto pNodes = net.CalcTopoOrd();
    oldFos.resize(net.GetIdMaxPlus1());
    appFos.resize(net.GetIdMaxPlus1());
    for (auto pNode: pNodes) {
        for (int i = 0; i < net.GetFanoutNum(pNode); ++i)
            oldFos[pNode->Id].emplace_back(net.GetFanout(pNode, i));
        if (net.GetName(pNode).ends_with("_app")) {
            for (int i = 0; i < net.GetFanoutNum(pNode); ++i) {
                auto pFanout = net.GetFanout(pNode, i);
                if (!net.IsObjPo(pFanout) && net.GetName(pFanout).ends_with("_app"))
                    appFos[pNode->Id].emplace_back(pFanout);
            }
        }
    }

    unordered_map<int, int> nFanouts;
    vector<pair<int, int>> nodeIdAndFoNumb; nodeIdAndFoNumb.reserve(pNodes.size());
    int countNodes = 0;
    for (auto it = pNodes.rbegin(); it != pNodes.rend(); ++it) {
        auto pNode = *it;
        if (!net.GetName(pNode).ends_with("_app"))
            continue;
        auto& pFanouts = appFos[pNode->Id];
        int nFo = pFanouts.size();
        if (!nFanouts.count(nFo))
            nFanouts[nFo] = 1;
        else
            nFanouts[nFo] += 1;
        if (nFo != 0)
            nodeIdAndFoNumb.emplace_back(pNode->Id, nFo);
        ++countNodes;
    }
    // cout << "#fanouts: " << endl;
    // for (auto& [key, value]: nFanouts)
    //     cout << key << "," << value << endl;
    std::sort(nodeIdAndFoNumb.begin(), nodeIdAndFoNumb.end(), [](const auto &a, const auto &b) {return a.second > b.second;});
    cout << "#nodes = " << countNodes << endl;
    cout << "#internal nodes = " << nodeIdAndFoNumb.size() << endl;
    cout << "#internal nodes using exact PBDs = " << nodeIdAndFoNumb.size() * exactPBDPerc << endl;

    IntVect useAccPBD(net.GetIdMaxPlus1(), 0);
    for (int i = 0; i < nodeIdAndFoNumb.size() * exactPBDPerc; ++i)
        useAccPBD[nodeIdAndFoNumb[i].first] = 1;

    // add PBD skeleton 
    auto pFPo = net.GetPo(net.GetPoNum() - 1);
    auto pF = net.GetFanin(pFPo, 0);
    assert(string(Abc_ObjName(pFPo)) == "f");
    auto pManFunc = (Abc_Aig_t *)(net.GetNet()->pManFunc);
    auto pConst1 = Abc_AigConst1(net.GetNet());
    auto pConst0 = Abc_ObjNot(pConst1);
    unordered_map<AbcObj*, AbcObj*> n2DFN;
    for (auto it = pNodes.rbegin(); it != pNodes.rend(); ++it) {
        auto pNode = *it;
        if (!net.GetName(pNode).ends_with("_app"))
            continue;
        auto& pFanouts = appFos[pNode->Id];
        int nFo = pFanouts.size();
        AbcObj* pDFN = nullptr;
        if (nFo == 0 || useAccPBD[pNode->Id]) { // if pNode is a PO; or pNode is an internal node, and pNode uses accurate PBD
            unordered_map<AbcObj*, AbcObj*> old2New;
            old2New[pNode] = Abc_ObjNot(pNode);
            auto tfo = GetTFO(pNode);
            for (auto pTfoNode: tfo) {
                auto fi0 = Abc_ObjFanin0(pTfoNode);
                if (old2New.count(fi0))
                    fi0 = old2New[fi0];
                auto fi1 = Abc_ObjFanin1(pTfoNode);
                if (old2New.count(fi1))
                    fi1 = old2New[fi1];
                old2New[pTfoNode] = Abc_AigAnd(pManFunc, Abc_ObjNotCond(fi0, Abc_ObjFaninC0(pTfoNode)), Abc_ObjNotCond(fi1, Abc_ObjFaninC1(pTfoNode)));
                if (!Abc_ObjIsComplement(old2New[pTfoNode]))
                    RenameAbcObj(old2New[pTfoNode], string(Abc_ObjName(pNode)) + "_tfo_" + string(Abc_ObjName(pTfoNode)) + "_tfo");
            }
            if (!old2New.count(pF)) {
                cout << Abc_ObjName(pNode) << endl;
                assert(0);
            }
            pDFN = Abc_AigXor(pManFunc, pF, old2New[pF]);
        }
        else if (nFo == 1)  { // if pNode is an internal node, and pNode uses approximate PBD, and #fanout = 1
            assert(pFanouts[0] != nullptr && !Abc_ObjIsPo(pFanouts[0]));
            assert(n2DFN.count(pFanouts[0]));
            pDFN = n2DFN[pFanouts[0]];
        }
        else { // if pNode is an internal node, and pNode uses approximate PBD, and #fanout > 1
            AbcObjVect dVNs(nFo, nullptr);
            for (int i = 0; i < nFo; ++i)
                dVNs[i] = GetLocPBD(pFanouts[i], pNode);
            AbcObjVect dFVs(nFo, nullptr);
            for (int i = 0; i < nFo; ++i) {
                assert(n2DFN.count(pFanouts[i]));
                dFVs[i] = n2DFN[pFanouts[i]];
            }
            auto pPart0 = pConst0;
            AbcObjVect bets(nFo, nullptr);
            for (int i = 0; i < nFo; ++i) {
                auto& bet = bets[i]; bet = pConst1;
                for (int j = 0; j < nFo; ++j)
                    bet = Abc_AigAnd(pManFunc, bet, (i == j)? dVNs[j]: Abc_ObjNot(dVNs[j]));
                pPart0 = Abc_AigOr(pManFunc, pPart0, Abc_AigAnd(pManFunc, bet, dFVs[i]));
            }
            auto alpha = pConst1;
            for (int j = 0; j < nFo; ++j)
                alpha = Abc_AigAnd(pManFunc, alpha, Abc_ObjNot(dVNs[j]));
            auto pPart1 = alpha;
            for (int i = 0; i < nFo; ++i)
                pPart1 = Abc_AigOr(pManFunc, pPart1, bets[i]);
            pDFN = Abc_AigOr(pManFunc, pPart0, Abc_ObjNot(pPart1));
        }
        n2DFN[pNode] = pDFN;
        net.CreatePo(pDFN, ("dF_" + net.GetName(pNode)).c_str());
    }
    Abc_NtkDeleteObj(pFPo);

    // merge PIs
    int nPo = net.GetPiNum() / 2;
    assert((net.GetPiNum() & 1) == 0);
    AbcObjVect pis;
    for (int i = 0; i < nPo; ++i) {
        auto pAccPi = net.GetPi(i);
        auto pAppPi = net.GetPi(nPo + i);
        Abc_ObjTransferFanout(pAccPi, pAppPi);
        assert(Abc_ObjFanoutNum(pAccPi) == 0);
        pis.emplace_back(pAccPi);
    }
    for (auto pi: pis)
        Abc_NtkDeleteObj(pi);

    // optimize
    cout << "start SAT sweeping the skeleton network" << endl;
    // net.SATSweep();
    net.Comm("ps; st; ps; ifraig; ps;");
    cout << "finish SAT sweeping the skeleton network" << endl;
}


int PBDMan::Synth(int fSASIMI) {
    {
    // test trivial cases, $\Delta f(n) \equiv 0$
    auto pNet = net.GetNet();
    auto pManFunc = (Abc_Aig_t *)(pNet->pManFunc);
    auto pConst1 = Abc_AigConst1(net.GetNet());
    auto pConst0 = Abc_ObjNot(pConst1);
    for (int i = 0; i < net.GetPoNum(); ++i) {
        auto pN = net.GetPo(i);
        auto name = net.GetName(pN);
        if (!name.ends_with("_app") || name.starts_with("dF_"))
            continue;
        auto pDFN = Abc_NtkFindCo(pNet, const_cast<char*>(("dF_" + net.GetName(pN)).c_str()));
        assert(pDFN != nullptr);
        auto pDFNDriv = Abc_ObjFanin0(pDFN); assert(!Abc_ObjIsComplement(pDFNDriv));
        if (pDFNDriv == pConst1 && Abc_ObjFaninC0(pDFN)) {
            auto pNDriv = Abc_ObjFanin0(pN); assert(!Abc_ObjIsComplement(pNDriv));
            cout << net.GetName(pN) << ", dFN \\equiv 0, replace its driver by constant 0" << endl;
            Abc_AigReplace(pManFunc, pNDriv, pConst0, 0);
            return 1;
        }
    }
    }

    {
    // get target nodes for constant
    auto pNet = net.GetNet();
    auto pManFunc = (Abc_Aig_t *)(pNet->pManFunc);
    vector<std::pair<AbcObj*, AbcObj*>> targNodes; targNodes.reserve(net.GetPoNum() / 2);
    for (int i = 0; i < net.GetPoNum(); ++i) {
        auto pN = net.GetPo(i);
        auto name = net.GetName(pN);
        if (!name.ends_with("_app") || name.starts_with("dF_"))
            continue;
        auto pDFN = Abc_NtkFindCo(pNet, const_cast<char*>(("dF_" + net.GetName(pN)).c_str())); assert(pDFN != nullptr);
        auto pDFNDriv = Abc_ObjFanin0(pDFN); assert(!Abc_ObjIsComplement(pDFNDriv));
        targNodes.emplace_back(pN, pDFN);
    }
    // test constant by simulation
    Simulator smlt(net, 19960822, 1 << 16);
    smlt.GenInpUnifFast();
    smlt.UpdNodeAndPoPatts();
    for (auto& targNode: targNodes) {
        auto pN = targNode.first, pDFN = targNode.second;
        auto pNDriv = Abc_ObjFanin0(pN); assert(!Abc_ObjIsComplement(pNDriv));
        auto pDFNDriv = Abc_ObjFanin0(pDFN); assert(!Abc_ObjIsComplement(pDFNDriv));
        const auto &nDat = smlt.GetDat(pN->Id), dFNDat = smlt.GetDat(pDFN->Id);
        // try constant
        if ((nDat & dFNDat).none()) {
            auto pFi0 = Abc_ObjFaninC0(pN)? Abc_ObjNot(pNDriv): pNDriv;
            auto pFi1 = Abc_ObjFaninC0(pDFN)? Abc_ObjNot(pDFNDriv): pDFNDriv;
            auto pDFS = Abc_AigAnd(pManFunc, pFi0, pFi1);
            cout << Abc_ObjName(pN) << ", const0 (simulation)" << endl;
            net.CreatePo(pDFS, ("ver_" + net.GetName(pDFN) + "*const0").c_str());
        }
        else if ((~nDat & dFNDat).none()) {
            auto pFi0 = Abc_ObjFaninC0(pN)? pNDriv: Abc_ObjNot(pNDriv);
            auto pFi1 = Abc_ObjFaninC0(pDFN)? Abc_ObjNot(pDFNDriv): pDFNDriv;
            auto pDFS = Abc_AigAnd(pManFunc, pFi0, pFi1);
            cout << Abc_ObjName(pN) << ", const1 (simulation)" << endl;
            net.CreatePo(pDFS, ("ver_" + net.GetName(pDFN) + "*const1").c_str());
        }
    }
    }
    {
    // test constant by SAT
    // net.SATSweep();
    net.Comm("ps; st; ps; ifraig; ps;");
    auto pNet = net.GetNet();
    auto pManFunc = (Abc_Aig_t *)(pNet->pManFunc);
    auto pConst1 = Abc_AigConst1(net.GetNet());
    auto pConst0 = Abc_ObjNot(pConst1);
    for (int i = 0; i < net.GetPoNum(); ++i) {
        auto pDFS = net.GetPo(i);
        auto name = net.GetName(pDFS);
        if (!name.starts_with("ver_"))
            continue;
        auto pDFSDriv = Abc_ObjFanin0(pDFS); assert(!Abc_ObjIsComplement(pDFSDriv));
        if (!(pDFSDriv == pConst1 && Abc_ObjFaninC0(pDFS)))
            continue;
        net.PrintObj(pDFS, true);
        auto pos0 = name.find("ver_dF_"), pos1 = name.find("*");
        assert(pos0 != std::string::npos && pos1 != std::string::npos);
        auto nName = name.substr(pos0 + 7, pos1 - pos0 - 7);
        auto pN = Abc_NtkFindCo(pNet, const_cast<char*>(nName.c_str())); assert(pN != nullptr);
        auto pNDriv = Abc_ObjFanin0(pN); assert(!Abc_ObjIsComplement(pNDriv));
        if (name.ends_with("const0")) {
            cout << Abc_ObjName(pN) << ", replace it by constant 0" << endl;
            if (!Abc_ObjFaninC0(pN))
                Abc_AigReplace(pManFunc, pNDriv, pConst0, 0);
            else
                Abc_AigReplace(pManFunc, pNDriv, pConst1, 0);
        }
        else if (name.ends_with("const1")) {
            cout << Abc_ObjName(pN) << ", replace it by constant 1" << endl;
            if (!Abc_ObjFaninC0(pN))
                Abc_AigReplace(pManFunc, pNDriv, pConst1, 0);
            else
                Abc_AigReplace(pManFunc, pNDriv, pConst0, 0);
        }
        return 1;
    }
    }

    if (!fSASIMI)
        return -1;
    {
    // get target nodes for SASIMI
    auto pNet = net.GetNet();
    auto pManFunc = (Abc_Aig_t *)(pNet->pManFunc);
    vector<std::pair<AbcObj*, AbcObj*>> targNodes; targNodes.reserve(net.GetPoNum() / 2);
    AbcObjVect subNodes; subNodes.reserve(net.GetPoNum() / 2);
    for (int i = 0; i < net.GetPoNum(); ++i) {
        auto pN = net.GetPo(i);
        auto name = net.GetName(pN);
        if (!name.ends_with("_app") || name.starts_with("dF_"))
            continue;
        auto pDFN = Abc_NtkFindCo(pNet, const_cast<char*>(("dF_" + net.GetName(pN)).c_str())); assert(pDFN != nullptr);
        auto pDFNDriv = Abc_ObjFanin0(pDFN); assert(!Abc_ObjIsComplement(pDFNDriv));
        targNodes.emplace_back(pN, pDFN);
        subNodes.emplace_back(pN);
    }
    for (int i = 0; i < net.GetPiNum(); ++i) {
        if (net.GetPiName(i).ends_with("_app"))
            subNodes.emplace_back(net.GetPi(i));
    }
    // test SASIMI by simulation
    net.GetLev();
    Simulator smlt(net, 19960822, 1 << 16);
    smlt.GenInpUnifFast();
    smlt.UpdNodeAndPoPatts();
    boost::timer::progress_display pd(targNodes.size());
    const int LAC_LIMIT_ON_NODE = 64;
    // const int LAC_LIMIT_ON_NODE = 16;
    int nLacs = 0;
    for (const auto& targNode: targNodes) {
        auto pN = targNode.first, pDFN = targNode.second;
        auto pNDriv = Abc_ObjFanin0(pN), pDFNDriv = Abc_ObjFanin0(pDFN);
        const auto &nDat = smlt.GetDat(pN->Id), dFNDat = smlt.GetDat(pDFN->Id);
        int lacCount = 0;
        for (const auto& pC: subNodes) {
            if (Abc_ObjIsPo(pC)) {
                auto pCDriv = Abc_ObjFanin0(pC);
                if (pCDriv->Level > pNDriv->Level || pCDriv == pNDriv)
                    continue;
            }
            else
                assert(Abc_ObjIsPi(pC));
            const auto &cDat = smlt.GetDat(pC->Id);
            if (((nDat ^ cDat) & dFNDat).none()) {
                // cout << pN << "," << pC << " buf (simulation)" << endl;
                AbcObj* pFi0 = Abc_ObjFaninC0(pN)? Abc_ObjNot(pNDriv): pNDriv;
                AbcObj* pFi1 = nullptr;
                if (Abc_ObjIsPo(pC)) {
                    auto pCDriv = Abc_ObjFanin0(pC);
                    pFi1 = Abc_ObjFaninC0(pC)? Abc_ObjNot(pCDriv): pCDriv;
                }
                else if (Abc_ObjIsPi(pC))
                    pFi1 = pC;
                else
                    assert(0);
                AbcObj* pFi2 = Abc_ObjFaninC0(pDFN)? Abc_ObjNot(pDFNDriv): pDFNDriv;
                auto pXor = Abc_AigXor(pManFunc, pFi0, pFi1);
                auto pDFS = Abc_AigAnd(pManFunc, pXor, pFi2);
                net.CreatePo(pDFS, ("sab_" + net.GetName(pDFN) + "*" + net.GetName(pC)).c_str());
                ++lacCount;
                ++nLacs;
            }
            else if (((nDat ^ ~cDat) & dFNDat).none()) {
                // cout << pN << "," << pC << " inv (simulation)" << endl;
                AbcObj* pFi0 = Abc_ObjFaninC0(pN)? Abc_ObjNot(pNDriv): pNDriv;
                AbcObj* pFi1 = nullptr;
                if (Abc_ObjIsPo(pC)) {
                    auto pCDriv = Abc_ObjFanin0(pC);
                    pFi1 = Abc_ObjFaninC0(pC)? pCDriv: Abc_ObjNot(pCDriv);
                }
                else if (Abc_ObjIsPi(pC))
                    pFi1 = Abc_ObjNot(pC);
                else
                    assert(0);
                AbcObj* pFi2 = Abc_ObjFaninC0(pDFN)? Abc_ObjNot(pDFNDriv): pDFNDriv;
                auto pXor = Abc_AigXor(pManFunc, pFi0, pFi1);
                auto pDFS = Abc_AigAnd(pManFunc, pXor, pFi2);
                net.CreatePo(pDFS, ("sai_" + net.GetName(pDFN) + "*" + net.GetName(pC)).c_str());
                ++lacCount;
                ++nLacs;
            }
            if (lacCount > LAC_LIMIT_ON_NODE)
                break;
        }
        ++pd;
    }
    cout << "nLacs = " << nLacs << endl;
    }
    {
    // test SASIMI by SAT
    // net.SATSweep();
    net.Comm("ps; st; ps; ifraig; ps;");
    auto pNet = net.GetNet();
    auto pManFunc = (Abc_Aig_t *)(pNet->pManFunc);
    auto pConst1 = Abc_AigConst1(net.GetNet());
    for (int i = 0; i < net.GetPoNum(); ++i) {
        auto pDFS = net.GetPo(i);
        auto name = net.GetName(pDFS);
        if (!name.starts_with("sa"))
            continue;
        auto pDFSDriv = Abc_ObjFanin0(pDFS); assert(!Abc_ObjIsComplement(pDFSDriv));
        if (!(pDFSDriv == pConst1 && Abc_ObjFaninC0(pDFS)))
            continue;
        net.PrintObj(pDFS, true);
        auto pos0 = name.find("sa"), pos1 = name.find("*");
        assert(pos0 != std::string::npos && pos1 != std::string::npos);
        auto nName = name.substr(pos0 + 7, pos1 - pos0 - 7);
        auto pN = Abc_NtkFindCo(pNet, const_cast<char*>(nName.c_str())); assert(pN != nullptr);
        auto pNDriv = Abc_ObjFanin0(pN); assert(!Abc_ObjIsComplement(pNDriv));
        auto cName = name.substr(pos1 + 1);
        auto pC = Abc_NtkFindCo(pNet, const_cast<char*>(cName.c_str()));
        if (pC == nullptr) pC = Abc_NtkFindCi(pNet, const_cast<char*>(cName.c_str())); 
        assert(pC != nullptr);
        cout << Abc_ObjName(pN) << "," << Abc_ObjName(pC) << endl;
        if (name.starts_with("sab")) {
            if (Abc_ObjIsPi(pC)) {
                if (!Abc_ObjFaninC0(pN))
                    Abc_AigReplace(pManFunc, pNDriv, pC, 0);
                else
                    Abc_AigReplace(pManFunc, pNDriv, Abc_ObjNot(pC), 0);
            }
            else if (Abc_ObjIsPo(pC)) {
                auto pCDriv = Abc_ObjFanin0(pC);
                auto pSub = Abc_ObjFaninC0(pC)? Abc_ObjNot(pCDriv): pCDriv;
                if (!Abc_ObjFaninC0(pN))
                    Abc_AigReplace(pManFunc, pNDriv, pSub, 0);
                else
                    Abc_AigReplace(pManFunc, pNDriv, Abc_ObjNot(pSub), 0);
            }
            else
                assert(0);
        }
        else if (name.starts_with("sai")) {
            if (Abc_ObjIsPi(pC)) {
                if (!Abc_ObjFaninC0(pN))
                    Abc_AigReplace(pManFunc, pNDriv, Abc_ObjNot(pC), 0);
                else
                    Abc_AigReplace(pManFunc, pNDriv, pC, 0);
            }
            else if (Abc_ObjIsPo(pC)) {
                auto pCDriv = Abc_ObjFanin0(pC);
                auto pSub = Abc_ObjFaninC0(pC)? pCDriv: Abc_ObjNot(pCDriv);
                if (!Abc_ObjFaninC0(pN))
                    Abc_AigReplace(pManFunc, pNDriv, pSub, 0);
                else
                    Abc_AigReplace(pManFunc, pNDriv, Abc_ObjNot(pSub), 0);
            }
            else
                assert(0);
        }
        else
            assert(0);
        return 1;
    }
    }

    return -1;
}


NetMan PBDMan::PostProc() {
    AbcObjVect delPos; delPos.reserve(net.GetPoNum());
    for (int i = 0; i < net.GetPoNum(); ++i) {
        auto name = net.GetPoName(i);
        if (!name.ends_with("_po") || name.starts_with("dF_") || name.starts_with("ver_") || name.starts_with("sa"))
            delPos.emplace_back(net.GetPo(i));
    }
    for (auto delPo: delPos)
        Abc_NtkDeleteObj(delPo);
    for (int i = 0; i < net.GetPoNum(); ++i) {
        auto name = net.GetPoName(i);
        auto loc = name.find("_po");
        assert(loc != name.npos);
        name = name.erase(loc);
        RenameAbcObj(net.GetPo(i), name);
    }
    for (int i = 0; i < net.GetPiNum(); ++i) {
        auto name = net.GetPiName(i);
        auto loc = name.find("_app");
        assert(loc != name.npos);
        name = name.erase(loc);
        RenameAbcObj(net.GetPi(i), name);
    }
    cout << "current approximate net" << endl;
    net.Synth(ORIENT::DELAY);
    // net.SATSweep();
    net.PrintStat();
    return net;
}


AbcObjVect PBDMan::GetTFO(AbcObj* pObj) const {
    AbcObjVect nodes;
    nodes.reserve(net.GetNodeNum());
    net.SetNetNotTrav();
    for (auto pFanout: oldFos[pObj->Id]) {
        if (!net.GetObjTrav(pFanout))
            GetTFORec(pFanout, nodes);
    }
    reverse(nodes.begin(), nodes.end());
    return nodes;
}


void PBDMan::GetTFORec(AbcObj* pObj, AbcObjVect& nodes) const {
    if (!net.IsNode(pObj))
        return;
    net.SetObjTrav(pObj);
    for (auto pFanout: oldFos[pObj->Id]) {
        if (!net.GetObjTrav(pFanout))
            GetTFORec(pFanout, nodes);
    }
    nodes.emplace_back(pObj);
}


AbcObj* PBDMan::GetLocPBD(AbcObj* pV, AbcObj* pU) {
    assert(!Abc_ObjIsComplement(pV));
    if (Abc_ObjFanin0(pV) == pU)
        return Abc_ObjChild1(pV);
    else {
        assert(Abc_ObjFanin1(pV) == pU);
        return Abc_ObjChild0(pV);
    }
}



} // namespace mecals_v1