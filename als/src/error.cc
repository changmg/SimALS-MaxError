/**
 * @file error.cc
 * @author Chang Meng
 * @brief Error analysis for approximate circuits
 * @date 2024-12-29
 * 
 */

#include "error.h"


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
using CMSat::SATSolver;
using CMSat::Lit;
using CMSat::lbool;



/**
 * @brief Constructor of the error manager
 * 
 * @param netMan0      the accurate network
 * @param netMan1      the approximate network
 */
ErrMan::ErrMan(const NetMan& netMan0, const NetMan& netMan1):
    net0(netMan0), net1(netMan1), pSmlt0(nullptr), pSmlt1(nullptr), pErrMit(nullptr), pSolver(nullptr) {
    // ensure the POs are the same
    if (!ComparePo(net0, net1)) {
        fmt::print(stderr, "Error: different POs\n");
        assert(0);
    }
    // ensure the PIs are the same
    if (!ComparePi(net0, net1, true)) {
        fmt::print(stderr, "Error: different PIs\n");
        assert(0);
    }
}


/**
 * @brief Constructor of the error manager
 * 
 * @param netMan0      the accurate network
 * @param netMan1      the approximate network
 * @param devNet       the deviation network
 */
ErrMan::ErrMan(const NetMan& netMan0, const NetMan& netMan1, const NetMan& devNet):
    net0(netMan0), net1(netMan1), pSmlt0(nullptr), pSmlt1(nullptr) {
    // ensure the POs are the same
    if (!ComparePo(net0, net1)) {
        fmt::print(stderr, "Error: different POs\n");
        assert(0);
    }
    // ensure the PIs are the same
    if (!ComparePi(net0, net1, true)) {
        fmt::print(stderr, "Error: different PIs\n");
        assert(0);
    }
    // initialize the error miter and SAT solver
    pErrMit = BuildErrMit(net0, net1, devNet);
    pSolver = BuildSatSolver_Abc(*pErrMit, cnfVarIdOfIthPi);
}


/**
 * @brief Perform logic simulation for the two networks net0 and net1
 * 
 * @param  seed        the seed for random number generation
 * @param  nFrame      the number of simulation frames
 * @param  distrType   the distribution type
 * @retval pSmlt0      the simulator for accurate network
 * @retval pSmlt1      the simulator for approximate network
 * @return void
 */
void ErrMan::LogicSim(unsigned seed, int nFrame, DISTR_TYPE distrType) {
    if (pSmlt0 != nullptr || pSmlt1 != nullptr) {
        fmt::print(stderr, "Error: simulators should not be initialized\n");
        assert(0);
    }
    assert(net0.IsPIOSame(net1));
    pSmlt0 = std::make_shared<Simulator>(net0, seed, nFrame, distrType);
    pSmlt1 = std::make_shared<Simulator>(net1, seed, nFrame, distrType); 
    pSmlt0->LogicSim();
    pSmlt1->LogicSim();
}


/**
 * @brief Compute the maximum error
 * 
 * @param metrType metric type
 * @return long long the maximum error
 */
BigInt ErrMan::ComputeMaxErr(METR_TYPE metrType) {
    // check the metric type
    if (metrType != METR_TYPE::MAXED && metrType != METR_TYPE::MAXHD) {
        fmt::print(stderr, "Error: unsupport metric type\n");
        assert(0);
    }

    // get the output width
    int outWidth = net0.GetPoNum();
    if (outWidth != net1.GetPoNum()) {
        fmt::print(stderr, "Error: different output width\n");
        assert(0);
    }

    // generate deviation network
    auto pDevCompNet = GenDevCompNet(metrType, outWidth);

    // build an error miter
    auto pErrMit = BuildErrMit(net0, net1, *pDevCompNet);
    
    // build a SAT solver
    auto pSolver = BuildSatSolver_Abc(*pErrMit, cnfVarIdOfIthPi);

    // return the maximum error
    int refErrWidth = outWidth;
    if (metrType == METR_TYPE::MAXHD)
        refErrWidth = static_cast<int>(log2(outWidth)) + 1;
    return SolveSatsForMaxErrBinSearch(*pErrMit, *pSolver, refErrWidth);
}


/**
 * @brief Build an error miter
 * 
 * @param accNet   the accurate network
 * @param appNet   the approximate network
 * @param devNet   the deviation network
 * @param errMit   the error miter
 * @retval errMit  the error miter
 * @return void
 */
std::shared_ptr<NetMan> ErrMan::BuildErrMit(const NetMan& accNet, const NetMan& appNet, const NetMan& devNet) {
    // check & prepare
    if (accNet.GetNetType() != NET_TYPE::SOP) {
        fmt::print(stderr, "Error: the accurate network should be in SOP\n");
        assert(0);
    }
    if (appNet.GetNetType() != NET_TYPE::SOP) {
        fmt::print(stderr, "Error: the approximate network should be in SOP\n");
        assert(0);
    }
    if (devNet.GetNetType() != NET_TYPE::SOP) {
        fmt::print(stderr, "Error: the deviation network should be in SOP\n");
        assert(0);
    }
    assert(ComparePi(accNet, appNet, true) && ComparePo(accNet, appNet));
    auto pErrMit = std::make_shared<NetMan>();
    auto& errMit = *pErrMit;
    errMit.StartSopNet();
    errMit.RenameNet("error_miter");

    // copy accNet
    auto pErrMitNet = errMit.GetNet();
    auto pAccNet = accNet.GetNet();
    AbcObj* pObj = nullptr, *pFanin = nullptr;
    int i = 0, k = 0;
    Abc_NtkCleanCopy(pAccNet);
    unordered_map <string, AbcObj*> name2PI;
    Abc_NtkForEachPi(pAccNet, pObj, i) {
        Abc_NtkDupObj(pErrMitNet, pObj, 0);
        RenameAbcObj(pObj->pCopy, string(Abc_ObjName(pObj)));
        name2PI[string(Abc_ObjName(pObj))] = pObj->pCopy;
    }
    Abc_NtkForEachNode(pAccNet, pObj, i) {
        Abc_NtkDupObj(pErrMitNet, pObj, 0);
        RenameAbcObj(pObj->pCopy, string(Abc_ObjName(pObj)) + "_acc");
    }
    Abc_NtkForEachNode(pAccNet, pObj, i) {
        Abc_ObjForEachFanin(pObj, pFanin, k)
            Abc_ObjAddFanin(pObj->pCopy, pFanin->pCopy);
    }

    // copy appNet
    auto pAppNet = appNet.GetNet();
    Abc_NtkCleanCopy(pAppNet);
    Abc_NtkForEachPi(pAppNet, pObj, i) {
        if (name2PI.count(string(Abc_ObjName(pObj)))) // the PI is in accNet
            pObj->pCopy = name2PI[string(Abc_ObjName(pObj))];
        else { // the PI is unique in appNet (a control signal)
            // fmt::print("Warning: control signal {}\n", Abc_ObjName(pObj));
            Abc_NtkDupObj(pErrMitNet, pObj, 0);
            RenameAbcObj(pObj->pCopy, string(Abc_ObjName(pObj)));
            name2PI[string(Abc_ObjName(pObj))] = pObj->pCopy;
        }
    }
    Abc_NtkForEachNode(pAppNet, pObj, i) {
        Abc_NtkDupObj(pErrMitNet, pObj, 0);
        RenameAbcObj(pObj->pCopy, string(Abc_ObjName(pObj)) + "_app");
    }
    Abc_NtkForEachNode(pAppNet, pObj, i) {
        Abc_ObjForEachFanin(pObj, pFanin, k)
            Abc_ObjAddFanin(pObj->pCopy, pFanin->pCopy);
    }

    // copy devNet
    auto pDevNet = devNet.GetNet();
    int nPo = accNet.GetPoNum();
    assert(appNet.GetPoNum() == nPo); 
    int hdRefErrWidth = static_cast<int>(log2(nPo)) + 1;
    // case 1: accNet PO + appNet PO + ref_err (for max error distance computation) 
    // case 2: accNet PO + appNet PO + hdRefErrWidth (for maximum hamming distance computation)
    // case 3: accNet PO + appNet PO (for error checking, no ref_err)
    assert(devNet.GetPiNum() == nPo * 3 || devNet.GetPiNum() == nPo * 2 + hdRefErrWidth || devNet.GetPiNum() == nPo * 2);
    assert(devNet.GetPoNum() >= 1);
    Abc_NtkCleanCopy(pDevNet);
    // link accNet POs to devNet PIs
    Abc_NtkForEachPo(pAccNet, pObj, i)
        Abc_NtkPi(pDevNet, i)->pCopy = Abc_ObjChild0Copy(pObj);
    // link appNet POs to devNet PIs
    Abc_NtkForEachPo(pAppNet, pObj, i)
        Abc_NtkPi(pDevNet, i + nPo)->pCopy = Abc_ObjChild0Copy(pObj);
    // copy ref_err signals (for maximum error distance computation & maximum hamming distance computation)
    if (devNet.GetPiNum() > nPo * 2) {
        int refErrWidth = Abc_NtkPiNum(pDevNet) - nPo * 2;
        for (int i = 0; i < refErrWidth; ++i) {
            auto pRefErrI = Abc_NtkPi(pDevNet, i + nPo * 2);
            Abc_NtkDupObj(pErrMitNet, pRefErrI, 0);
            RenameAbcObj(pRefErrI->pCopy, string(Abc_ObjName(pRefErrI)));
            assert(string(Abc_ObjName(pRefErrI)) == fmt::format("ref_err[{}]", i));
        }
    }
    // copy nodes
    Abc_NtkForEachNode(pDevNet, pObj, i) {
        Abc_NtkDupObj(pErrMitNet, pObj, 0);
        RenameAbcObj(pObj->pCopy, string(Abc_ObjName(pObj)) + "_mit");
    }
    Abc_NtkForEachNode(pDevNet, pObj, i) {
        Abc_ObjForEachFanin(pObj, pFanin, k)
            Abc_ObjAddFanin(pObj->pCopy, pFanin->pCopy);
    }
    // copy POs
    Abc_NtkForEachPo(pDevNet, pObj, i)
        Abc_NtkDupObj(pErrMitNet, pObj, 1);
    Abc_NtkForEachPo(pDevNet, pObj, i)
        Abc_ObjAddFanin(pObj->pCopy, Abc_ObjChild0Copy(pObj));

    return pErrMit;
}


/**
 * @brief Build SAT solver from the network
 * 
 * @param net      the network defining the SAT problem
 * @param solver   the SAT solver
 * @retval solver  the SAT solver
 * @return void
 */
std::shared_ptr<CMSat::SATSolver> ErrMan::BuildSatSolver_Naive(NetMan& net) {
    // convert the network to strashed AIG
    if (!net.IsStrash()) {
        // net.Comm("ps; st; ps;");
        // net.Comm("ps; st; ps; ifraig; ps;");
        net.Comm("st; ifraig;");
    }

    // start solver
    auto pSolver = std::make_shared<CMSat::SATSolver>();
    // solver.set_num_threads(nThreads);
    pSolver->set_max_confl(1ll << 20);

    // add variables: abc id = cmsat var id
    pSolver->new_vars(net.GetIdMaxPlus1()); 

    // add clauses
    // add clause for constant 1
    int constId = Abc_AigConst1(net.GetNet())->Id;
    vector<Lit> clause;
    clause.emplace_back(Lit(constId, false));
    pSolver->add_clause(clause);
    
    // add clause for each AND gate
    AbcObj* pObj = nullptr;
    int i = 0;
    Abc_AigForEachAnd(net.GetNet(), pObj, i) {
        // get fanin ids and their phases
        int fanin0 = Abc_ObjFanin0(pObj)->Id, faninC0 = Abc_ObjFaninC0(pObj);
        int fanin1 = Abc_ObjFanin1(pObj)->Id, faninC1 = Abc_ObjFaninC1(pObj);
        // function: pObj = (fanin0 ^ faninC0) & (fanin1 ^ faninC1)
        // clause 1: (faninC0? !fanin0: fanin0) | !pObj
        clause.clear();
        clause.emplace_back(Lit(fanin0, faninC0));
        clause.emplace_back(Lit(pObj->Id, true));
        pSolver->add_clause(clause);
        // clause 2: (faninC1? !fanin1: fanin1) | !pObj
        clause.clear();
        clause.emplace_back(Lit(fanin1, faninC1));
        clause.emplace_back(Lit(pObj->Id, true));
        pSolver->add_clause(clause);
        // clause 3: (faninC0? fanin0: !fanin0) | (faninC1? fanin1: !fanin1) | pObj
        clause.clear();
        clause.emplace_back(Lit(fanin0, !faninC0));
        clause.emplace_back(Lit(fanin1, !faninC1));
        clause.emplace_back(Lit(pObj->Id, false));
        pSolver->add_clause(clause);
    }

    // add clause for the PO: assert po == 1
    if (net.GetPoNum() != 1) {
        fmt::print(stderr, "Error: the network defining the SAT problem should have only one PO\n");
        assert(0);
    }
    AbcObj* pPo = Abc_NtkPo(net.GetNet(), 0);
    int driverId = Abc_ObjFanin0(pPo)->Id, driverC = Abc_ObjFaninC0(pPo);
    // fmt::print("PO driver: id = {}, phase = {}\n", driverId, driverC);
    clause.clear();
    clause.emplace_back(Lit(driverId, driverC));
    pSolver->add_clause(clause);

    // update the CNF variable ID of the i-th PI of the error miter
    // cnfVarId = abcId
    cnfVarIdOfIthPi.clear();
    cnfVarIdOfIthPi.resize(net.GetPiNum());
    for (int i = 0; i < net.GetPiNum(); ++i)
        cnfVarIdOfIthPi[i] = net.GetPiId(i);

    return pSolver;
}


/**
 * @brief Convert the literal to the variable
 * @brief Auxiliary function for BuildSatSolver_Abc
 * 
 * @param Lit     the literal
 * @return int    the variable
 */
static inline int Cnf_Lit2Var(int Lit) {return (Lit & 1)? -(Lit >> 1) - 1 : (Lit >> 1) + 1;}


/**
 * @brief Build SAT solver from the network
 * @brief Copied from ABC "&write_cnf" command
 * 
 * @param net      the network defining the SAT problem
 * @return a pointer to the SAT solver
 */
std::shared_ptr<CMSat::SATSolver> ErrMan::BuildSatSolver_Abc(NetMan& net, IntVect& cnfVarIdOfIthPi) {
    // preprocess the network
    if (net.GetPoNum() != 1) {
        fmt::print(stderr, "Error: the network defining the SAT problem should have only one PO\n");
        assert(0);
    }
    // net.Comm("ps; st; ps; ifraig; &get; &ps");
    net.Comm("st; ifraig; &get");
    auto pGia = AbcMan().GetAbcFrame()->pGia;

    // default parameters
    int nLutSize    = 8;
    int fCnfObjIds  = 0;
    int fAddOrCla   = 1;
    abc::Jf_Par_t pars;
    assert( nLutSize >= 3 && nLutSize <= 8 );
    abc::Mf_ManSetDefaultPars(&pars);
    pars.fGenCnf     = 1;
    pars.fCoarsen    = !fCnfObjIds;
    pars.nLutSize    = nLutSize;
    pars.fCnfObjIds  = fCnfObjIds;
    pars.fAddOrCla   = fAddOrCla;
    pars.fCnfMapping = 0;
    pars.fVerbose    = 0;

    // convert the GIA network to CNF
    // CNF variable mapping rules:
    // Assume CNF has N variables, with variable IDs running from 0 to N-1.
    // Variable number 0 is not used in the CNF.
    // Variables 1, 2, 3,... <nPOs> represent POs in their natural order.
    // Variables N-<nPIs>, N-<nPIs>+1, N-<nPIs>+2, ... N-1, represent PIs in their natural order.
    // The internal variables are ordered in a reverse topological order from outputs to inputs.
    // That is, smaller variable IDs tend to be closer to the outputs, while larger
    // variable IDs tend to be closer to the inputs. It was found that this ordering
    // leads to faster SAT solving for hard UNSAT CEC problems.
    auto pNew = abc::Mf_ManPerformMapping(pGia, &pars);
    abc::Gia_ManStopP(&pNew);
    auto pCnf = (abc::Cnf_Dat_t *)pGia->pData;

    // start solver
    auto pSolver = std::make_shared<CMSat::SATSolver>();
    // pSolver->set_max_confl(1ll << 20);
    // pSolver->set_max_confl(1ll << 16);
    pSolver->set_max_confl(1ll << 18);
    pSolver->new_vars(pCnf->nVars + 1); // in the ABC cnf, the variable id starts from 1

    // add clauses
    vector<Lit> clause;
    for (int i = 0; i < pCnf->nClauses; i++) {
        clause.clear();
        for (auto pLit = pCnf->pClauses[i], pStop = pCnf->pClauses[i+1]; pLit < pStop; pLit++) {
            int var = Cnf_Lit2Var(*pLit);
            clause.emplace_back(Lit(abs(var), var < 0));
        }
        pSolver->add_clause(clause);
    }

    // update the CNF variable ID of the i-th PI of the error miter
    assert(pSolver->nVars() - net.GetPiNum() >= 1);
    cnfVarIdOfIthPi.clear();
    cnfVarIdOfIthPi.resize(net.GetPiNum());
    for (int i = 0; i < net.GetPiNum(); ++i) {
        // Assume CNF has N variables, with variable IDs running from 0 to N-1; Variables N-<nPIs>, N-<nPIs>+1, N-<nPIs>+2, ... N-1, represent PIs in their natural order.
        cnfVarIdOfIthPi[i] = pSolver->nVars() - net.GetPiNum() + i;
    }
    
    // clean up
    Cnf_DataFree(pCnf);

    return pSolver;
}


/**
 * @brief Get the maximum error defined by the error miter and the corresponding SAT solver using binary search
 * 
 * @param net        the error miter
 * @param solver     the SAT solver corresponding to the error miter
 * @param refEdWidth the width of the "ref_err" signal
 * @return BigInt    the maximum error
 */
BigInt ErrMan::SolveSatsForMaxErrBinSearch(NetMan& net, CMSat::SATSolver& solver, int refEdWidth) {
    // init
    int netPiNum = net.GetPiNum();
    if (refEdWidth >= 500 || refEdWidth >= netPiNum) {
        fmt::print(stderr, "Error: the reference error width is too large\n");
        assert(0);
    }
    int refEdStartPi = netPiNum - refEdWidth;
    BigInt maxRefEd = (BigInt(1) << refEdWidth) - 1;
    vector<Lit> assumptions(refEdWidth, Lit(0, true));

    // solve
    BigInt left = 0, right = maxRefEd;
    assert(static_cast<int>(cnfVarIdOfIthPi.size()) == netPiNum);
    while (left <= right) {
        BigInt mid = left + (right - left) / 2;
        // set the assumptions
        BigInt refEdTmp = mid;
        for (int i = 0; i < refEdWidth; ++i) {
            auto piCnfVarId = cnfVarIdOfIthPi[refEdStartPi + i];
            assumptions[i] = Lit(piCnfVarId, !(refEdTmp & 1));
            refEdTmp >>= 1;
        }
        assert(refEdTmp == 0);
        // solve
        lbool res = solver.solve(&assumptions);
        // fmt::print("left = {}, right = {}, mid = {}, res = {}\n", left.str(), right.str(), mid.str(), res);
        // check the result
        if (res == CMSat::l_False) {
            right = mid - 1;
        } else if (res == CMSat::l_True) {
            left = mid + 1;
        } else {
            fmt::print("Error: SAT solver returns undefined\n");
            assert(0);
        }
    }
    return left;
}


/**
 * @brief Add a unit clause in the SAT solver for the i-th PI of the error miter
 * 
 * @param iPi           index of the PI
 * @param fVarCompl     the complement of the variable
 * @retval pSolver      the SAT solver
 * @return void
 */
void ErrMan::AddUnitClauseOfPi(int iPi, bool fVarCompl) { 
    assert(pSolver != nullptr); 
    int cnfVarId = GetCnfVarIdOfIthPi(iPi); 
    std::vector<CMSat::Lit> clause{CMSat::Lit(cnfVarId, fVarCompl)}; 
    pSolver->add_clause(clause);
}


// /**
//  * @brief Pick the first LAC whose real maximum error is no more than the given bound
//  * @brief Assume that the LACs are sorted: primary key = (smaller) error, secondary key = (larger) sizeGain
//  * 
//  * @param  lacMan         the LAC manager
//  * @param  accSmlt        the accurate network's simulator
//  * @param  appNet         the approximate network
//  * @param  devNet         the deviation network
//  * @param  lacBlackList   the black list of LACs, which cause the SAT solver to return undefined
//  * @return the first LAC whose real maximum error is no more than the given bound; nullptr if no such LAC
//  */
// std::shared_ptr<LAC> BatchErrEst::PickFirstValidLac(LACMan& lacMan, Simulator& accSmlt, NetMan& appNet, const NetMan& devNet, std::unordered_set<std::string>& lacBlackList) {
//     fmt::print("Check the maximum error for each LAC using SAT\n");
//     // prepare the counter example
//     auto& accNet = static_cast<const NetMan&>(accSmlt);
//     IntVect counterEx;
//     counterEx.resize(accNet.GetPiNum());
//     bool hasCounterEx = false;
//     // identify the first LAC whose real maximum error is no more than the given bound
//     std::shared_ptr<LAC> pRetLac = nullptr;
//     IntVect replTrace;
//     for (int i = 0; i < lacMan.GetLacNum(); ++i) {
//         auto pLac = lacMan.GetLac(i);
//         fmt::print("{}{} (estimated){}\n", HALF_DASH_LINE, pLac->ToStr(), HALF_DASH_LINE);
//         TempApplyLac(appNet, *pLac, replTrace, false);
//         ErrMan errMan(accNet, appNet, devNet);
//         auto res = errMan.SolveSat(counterEx, true);
//         if (res == CMSat::l_False) { // UNSAT, satisfy the error bound
//             pRetLac = pLac;
//             fmt::print("Satisfy the error bound\n");
//             RecovNet(appNet, {replTrace}, false);
//             break;
//         }
//         else if (res == CMSat::l_True) { // SAT, exceed the error bound
//             fmt::print("Exceed the error bound, save the {}-th counter example\n", countExNum);
//             // save the counter example in the countExNum-th PI pattern in the accSmlt
//             assert(ComparePi(accNet, errMan.GetErrMit(), false));
//             accSmlt.ReplInp(countExNum, counterEx);
//             hasCounterEx = true;
//             ++countExNum;
//             if (countExNum >= nFrame)
//                 countExNum = 0;
//         }
//         else {
//             fmt::print("Warning: SAT solver returns undefined, skip this LAC and add it to the black list\n");
//             lacBlackList.insert(pLac->ToStrShort());
//         }
//         RecovNet(appNet, {replTrace}, false);
//     }
//     // re-simulate the accurate network if counter examples are generated
//     if (hasCounterEx)
//         accSmlt.UpdNodeAndPoPatts();
//     return pRetLac;
// }


/**
 * @brief Error estimation for multiple LACs (using enumeration)
 * 
 * @param  lacMan       the LAC manager
 * @param  accNet       the accurate network
 * @param  appNet       the approximate network
 * @param  errUppBound  the upper bound of the error
 * @retval lacMan       the LAC manager with the error information
 * @return void
 */
void BatchErrEst::CompLacErrsByEnumAndPruneBadLacs(LACMan& lacMan, Simulator& accSmlt, const NetMan& appNet, ll errUppBound) {
    fmt::print("Max error computation using enumeration\n");
    assert(accSmlt.GetPiNum() < 30);
    int enumFrame = 1 << accSmlt.GetPiNum();
    // rough estimation
    if (ROUGH_SIM_FRAME < enumFrame) {
        auto startTime = std::chrono::high_resolution_clock::now();
        PruneLacsWithSim(lacMan, accSmlt, appNet, errUppBound, ROUGH_SIM_FRAME, DISTR_TYPE::UNIF);
        PrintRuntime(startTime, "rough simulation");
    }
    // enumeration
    auto startTime = std::chrono::high_resolution_clock::now();
    PruneLacsWithSim(lacMan, accSmlt, appNet, errUppBound, enumFrame, DISTR_TYPE::ENUM);
    PrintRuntime(startTime, "enumeration");
}


/**
 * @brief Compute the maximum error lower bound for each LAC using logic simulation and prune large-error LACs
 * 
 * @param lacMan       the LAC manager
 * @param accSmlt      the accurate network's simulator
 * @param appNet       the approximate network
 * @param errUppBound  the upper bound of the error
 */
void BatchErrEst::CompLacErrsBySimAndPruneBadLacs(LACMan& lacMan, Simulator& accSmlt, const NetMan& appNet, ll errUppBound) {
    // rough estimation
    auto startTime = std::chrono::high_resolution_clock::now();
    PruneLacsWithSim(lacMan, accSmlt, appNet, errUppBound, ROUGH_SIM_FRAME, DISTR_TYPE::UNIF);
    PrintRuntime(startTime, "rough simulation");
    // fine estimation
    if (ROUGH_SIM_FRAME < nFrame) {
        startTime = std::chrono::high_resolution_clock::now();
        PruneLacsWithSim(lacMan, accSmlt, appNet, errUppBound, nFrame, DISTR_TYPE::UNIF);
        PrintRuntime(startTime, "fine-grained simulation");
    }
}


/**
 * @brief Prune large-error LACs using logic simulation
 * 
 * @param lacMan       the LAC manager
 * @param accNet       the accurate network
 * @param appNet       the approximate network
 * @param errUppBound  the upper bound of the error
 */
void BatchErrEst::PruneLacsWithSim(LACMan& lacMan, Simulator& acc_smlt, const NetMan& appNet, ll errUppBound, int nFramePrune, DISTR_TYPE distrType) {
    // regroup LACs by node
    lacMan.RegroupLACsByNode(true);
    auto& node2Lacs = lacMan.GetNode2Lacs();

    // perform logic simulation
    auto& accNet = static_cast<const NetMan&>(acc_smlt);
    if (!accNet.IsPIOSame(appNet)) {
        fmt::print(stderr, "Error: different PI/PO\n");
        assert(0);
    }
    if (distrType == DISTR_TYPE::ENUM) {
        int nPi = appNet.GetPiNum();
        assert(nPi < 30);
        nFramePrune = 1 << nPi;
    }
    int nPo = appNet.GetPoNum();
    vector<BigInt> YAcc(nFramePrune, 0);
    vector<BitVect> accPos(nPo);
    Simulator appSmlt(appNet, seed, nFramePrune, distrType);
    if (nFramePrune == nFrame) {
        appSmlt.GenInpFromOthSmlt(acc_smlt);
        appSmlt.UpdNodeAndPoPatts();
        if (metrType == METR_TYPE::MAXED) {
            for (int iPatt = 0; iPatt < nFramePrune; ++iPatt)
                acc_smlt.GetOutput(iPatt, YAcc[iPatt]);
        }
        else if (metrType == METR_TYPE::MAXHD) {
            for (int iPo = 0; iPo < nPo; ++iPo)
                accPos[iPo] = acc_smlt.GetDat(acc_smlt.GetPoId(iPo));
        }
    }
    else {
        Simulator accSmltFewPatt(accNet, seed, nFramePrune, distrType);
        accSmltFewPatt.LogicSim();
        appSmlt.LogicSim();
        if (metrType == METR_TYPE::MAXED) {
            for (int iPatt = 0; iPatt < nFramePrune; ++iPatt)
                accSmltFewPatt.GetOutput(iPatt, YAcc[iPatt]);
        }
        else if (metrType == METR_TYPE::MAXHD) {
            for (int iPo = 0; iPo < nPo; ++iPo) {
                accPos[iPo] = accSmltFewPatt.GetDat(accSmltFewPatt.GetPoId(iPo));
                accPos[iPo].resize(nFramePrune);
            }
        }
    }

    // prune large-error LACs
    // prepare
    fmt::print("Compute maximum error lower bound for each of the {} LACs using {} simulation patterns\n", lacMan.GetLacNum(), nFramePrune);
    auto topoNodes = appNet.CalcTopoOrd(false);
    vector<BitVect> bdPosWrtNode;
    vector<BitVect> tempPos(nPo);
    vector<BitVect> poDiffs(nPo);
    BitVect nodePatt(nFramePrune, 0);
    bool useEnum = (distrType == DISTR_TYPE::ENUM);
    BigInt runMin(errUppBound);
    // iterate over the nodes
    boost::timer::progress_display pd(lacMan.GetLacNum());
    // fmt::print("#lacs = {}\n", lacMan.GetLacNum());
    // int count = 0;
    for (int iNode = 0; iNode < static_cast<int>(topoNodes.size()); ++iNode) {
        // skip the node if it is not associated with any LAC
        if (!node2Lacs.count(topoNodes[iNode]->Id))
            continue;
        // flip the i-th node and simulate the circuit
        appSmlt.CalcBoolDiff(topoNodes, iNode, bdPosWrtNode);
        // iterate over the LACs associated with the i-th node
        for (const auto& pLac: node2Lacs.at(topoNodes[iNode]->Id)) {
            ++pd;
            // ++count;
            int pTargId = pLac->GetTargId();
            if (pTargId != topoNodes[iNode]->Id) {
                fmt::print(stderr, "Error: inconsistent node id, pTargId = {}, topoNodes[iNode]->Id = {}\n", pTargId, topoNodes[iNode]->Id);
                assert(0);
            }
            // compute the simulation pattern of the target node of the LAC
            appSmlt.SimSop(pLac->GetDivIds(), pLac->GetSop(), nodePatt);
            // compute whether the node value change after applying the LAC
            auto& nodeChange = nodePatt;
            nodeChange ^= appSmlt.GetDat(pTargId);
            // get the PO patterns after applying the LAC
            for (int j = 0; j < nPo; ++j) {
                auto poId = appSmlt.GetPoId(j);
                tempPos[j] = appSmlt.GetDat(poId) ^ (nodeChange & bdPosWrtNode[j]); 
            }
            // get max error (obtained by simulation)
            if (metrType == METR_TYPE::MAXED) {
                BigInt maxErrSim(0), YNew(0);
                if (!useEnum) {
                    for (int iPatt = 0; iPatt < nFramePrune; ++iPatt) {
                        GetValue(tempPos, iPatt, YNew);
                        maxErrSim = std::max(maxErrSim, abs(YNew - YAcc[iPatt]));
                        if (maxErrSim > errUppBound)
                            break;
                    }
                }
                else { // enumeration
                    for (int iPatt = 0; iPatt < nFramePrune; ++iPatt) {
                        GetValue(tempPos, iPatt, YNew);
                        maxErrSim = std::max(maxErrSim, abs(YNew - YAcc[iPatt]));
                        // we only need to identify the LAC with the smallest error 
                        // so we can break if the error is already larger than the current minimum error
                        if (maxErrSim > runMin)
                            break;
                    }
                    runMin = std::min(runMin, maxErrSim);
                }
                pLac->SetErr(static_cast<double>(maxErrSim));
            }
            else if (metrType == METR_TYPE::MAXHD) {
                ll maxErrSim = 0;
                for (int iPo = 0; iPo < nPo; ++iPo)
                    poDiffs[iPo] = tempPos[iPo] ^ accPos[iPo];
                if (!useEnum) {
                    // BigInt YNew(0);
                    for (int iPatt = 0; iPatt < nFramePrune; ++iPatt) {
                        // GetValue(tempPos, iPatt, YNew);
                        ll hd = 0;
                        for (int iPo = 0; iPo < nPo; ++iPo)
                            hd += poDiffs[iPo][iPatt];
                        maxErrSim = std::max(maxErrSim, hd);
                        if (maxErrSim > errUppBound)
                            break;
                    }
                }
                else {
                    for (int iPatt = 0; iPatt < nFramePrune; ++iPatt) {
                        ll hd = 0;
                        for (int iPo = 0; iPo < nPo; ++iPo)
                            hd += poDiffs[iPo][iPatt];
                        maxErrSim = std::max(maxErrSim, hd);
                        if (maxErrSim > runMin)
                            break;
                    } 
                    runMin = std::min(runMin, static_cast<BigInt>(maxErrSim));
                }
                pLac->SetErr(static_cast<double>(maxErrSim));
            }
            else {
                fmt::print(stderr, "Error: unsupported metric type\n");
                assert(0);
            }
        }
    }
    // fmt::print("#LACs processed: {}\n", count);
    // assert(count == lacMan.GetLacNum());
    
    // update the LAC manager using the promising LACs
    lacMan.RemLargeErrLacs(static_cast<double>(errUppBound));
    fmt::print("#promising LACs after pruning: {}\n", lacMan.GetLacNum());
}


/**
 * @brief Sort the LACs by the upper bound of the maximum error
 * 
 * @param lacMan    the LAC manager
 * @param appNet    the approximate network
 * @retval lacMan   the LAC manager with the LACs sorted
 * @return void
 */
void BatchErrEst::CalcErrLooseUppBound(LACMan& lacMan, const NetMan& appNet) {
    fmt::print("Sort the LACs by the upper bound of the maximum error\n");
    // prepare
    auto topoNodes = appNet.CalcTopoOrdOfIds(false);
    DblVect uppBounds(appNet.GetIdMaxPlus1(), 0.0);
    lacMan.RegroupLACsByNode();

    // initialize the upper bound of the maximum error for each PO
    for (int i = 0; i < appNet.GetPoNum(); ++i) {
        auto poId = appNet.GetPoId(i);
        if (metrType == METR_TYPE::MAXED) // the increasing of MAXED caused by the change of the i-th appNet's PO is at most (1 << i)
            uppBounds[poId] = static_cast<double>(BigInt(1) << i);
        else if (metrType == METR_TYPE::MAXHD) // the increasing of MAXHD caused by the change of the i-th appNet's PO is at most 1
            uppBounds[poId] = 1.0; 
        else {
            fmt::print(stderr, "Error: unsupported metric type\n");
            assert(0);
        }
    }

    // traverse the nodes in reverse topological order
    for (auto itNode = topoNodes.rbegin(); itNode != topoNodes.rend(); ++itNode) {
        auto nodeId = *itNode;
        // upper bound of nodeId = the sum of the upper bounds of its fanouts
        for (int i = 0; i < appNet.GetFanoutNum(nodeId); ++i) {
            auto fanoutId = appNet.GetFanoutId(nodeId, i);
            uppBounds[nodeId] += uppBounds[fanoutId];
        }
        // update the upper bound of the LACs associated with the node
        if (lacMan.GetNode2Lacs().count(nodeId)) {
            for (const auto& pLac: lacMan.GetNode2Lacs().at(nodeId))
                pLac->SetErr2(uppBounds[nodeId]);
        }
    }
    
    // sort the LACs by the upper bound of the maximum error
    // lacMan.SortAndKeepTopKLACs(-1);
    // lacMan.PrintLACs(10);
}


/**
 * @brief Generate a deviation network 
 * @brief to compute the deviation between the accurate and approximate networks
 * 
 * @param metrType  metric type
 * @param outWidth  output width of the accurate and approximate networks
 * @return pDevNet  pointer to the generated network
 */
std::shared_ptr<NetMan> GenDevNet(METR_TYPE metrType, int outWidth) {
    // path to store the deviation network
    string folder = "./tmp";
    CreateDir(folder);

    // generate deviation network
    // prepare
    string devNetPath;
    string err;
    int refErrWidth = -1;
    assert(outWidth > 0);
    if (metrType == METR_TYPE::MAXED) {
        err = "error_distance";
        refErrWidth = outWidth;
    }
    else if (metrType == METR_TYPE::MAXHD) {
        err = "hamming_distance";
        refErrWidth = static_cast<int>(log2(outWidth)) + 1;
    }
    else {
        fmt::print(stderr, "Error: unsupport metric type\n");
        assert(0);
    }
    devNetPath = fmt::format("{}/{}_width{}", folder, err, outWidth);
    // if the deviation network file exists, use it
    if (IsPathExist(devNetPath + "_opt.blif"))
        fmt::print("Use the existing circuit in {} to compute {}\n", devNetPath + "_opt.blif", err);
    else { // otherwise, generate the deviation network
        fmt::print("Generating a circuit to compute {}\n", err);
        fmt::print("Currently, only support unsigned outputs\n");
        // generate the Verilog file
        auto fout = fmt::output_file(devNetPath + ".v");
        fout.print("module {}(a, b, {});\n", err, err);
        fout.print("parameter _bit = {};\n", outWidth);
        fout.print("input [_bit - 1: 0] a;\n");
        fout.print("input [_bit - 1: 0] b;\n");
        fout.print("output [{} - 1: 0] {};\n", refErrWidth, err);
        if (metrType == METR_TYPE::MAXED)
            fout.print("assign {} = (a > b)? (a - b): (b - a);\n", err);
        else if (metrType == METR_TYPE::MAXHD) {
            fout.print("wire [_bit - 1: 0] diff;\n");
            fout.print("assign diff = a ^ b;\n");
            fout.print("assign {} = ", err);
            for (int i = 0; i < outWidth; ++i) {
                if (i > 0)
                    fout.print(" + ");
                fout.print("diff[{}]", i);
            }
            fout.print(";\n");
        }
        else {
            fmt::print(stderr, "Error: unsupported metric type\n");
            assert(0);
        }
        fout.print("endmodule\n");
        fout.flush();
        // convert the Verilog file to BLIF using yosys
        string yosysComm = fmt::format("yosys -p \"read_verilog {}; synth; write_blif {}\" > {}", devNetPath + ".v", devNetPath + ".blif", devNetPath + "_yosys.log"); 
        fmt::print("Converting the Verilog file to BLIF\n");
        ExecSystComm(yosysComm);
        // convert the BLIF file to AIG
        AbcMan abc;
        string abcComm = fmt::format("r {}; st; resyn2rs; resyn2rs; resyn2rs; w {}", devNetPath + ".blif", devNetPath + "_opt.blif");
        fmt::print("Optimize the BLIF file\n");
        abc.Comm(abcComm, true);
    }

    // load the deviation network
    auto pDevNet = std::make_shared<NetMan>(devNetPath + "_opt.blif");
    pDevNet->Comm("ps");
    return pDevNet;
}


/**
 * @brief Generate a deviation & comparing network 
 * @brief to check whether the deviation between the accurate and approximate networks is no larger than the reference error
 * 
 * @param metrType  metric type
 * @param outWidth  output width of the accurate and approximate networks
 * @return pDevNet  pointer to the generated network
 */
std::shared_ptr<NetMan> GenDevCompNet(METR_TYPE metrType, int outWidth) {
    // path to store the deviation network
    string folder = "./tmp";
    CreateDir(folder);

    // generate deviation network
    // prepare
    string devCompNetPath;
    string err;
    int refErrWidth = -1;
    assert(outWidth > 0);
    if (metrType == METR_TYPE::MAXED) {
        err = "error_distance";
        refErrWidth = outWidth;
    }
    else if (metrType == METR_TYPE::MAXHD) {
        err = "hamming_distance";
        refErrWidth = static_cast<int>(log2(outWidth)) + 1;
    }
    else {
        fmt::print(stderr, "Error: unsupport metric type\n");
        assert(0);
    }
    devCompNetPath = fmt::format("{}/comp_{}_width{}", folder, err, outWidth);
    // if the deviation network file exists, use it
    if (IsPathExist(devCompNetPath + "_opt.blif"))
        fmt::print("Use the existing circuit in {} to compute {} and compare the error with the reference error\n", devCompNetPath + "_opt.blif", err);
    else { // otherwise, generate the deviation network
        fmt::print("Generating a circuit to compute {} and compare the error with the reference error\n", err);
        fmt::print("Currently, only support unsigned outputs\n");
        // generate the Verilog file
        auto fout = fmt::output_file(devCompNetPath + ".v");
        fout.print("module {}_devcomp(a, b, ref_err, f);\n", err);
        fout.print("parameter _bit = {};\n", outWidth);
        fout.print("input [_bit - 1: 0] a;\n");
        fout.print("input [_bit - 1: 0] b;\n");
        fout.print("input [{} - 1: 0] ref_err;\n", refErrWidth);
        fout.print("output f;\n");
        fout.print("wire [{} - 1: 0] {};\n", refErrWidth, err);
        if (metrType == METR_TYPE::MAXED)
            fout.print("assign {} = (a > b)? (a - b): (b - a);\n", err);
        else if (metrType == METR_TYPE::MAXHD) {
            fout.print("wire [_bit - 1: 0] diff;\n");
            fout.print("assign diff = a ^ b;\n");
            fout.print("assign {} = ", err);
            for (int i = 0; i < outWidth; ++i) {
                if (i > 0)
                    fout.print(" + ");
                fout.print("diff[{}]", i);
            }
            fout.print(";\n");
        }
        else {
            fmt::print(stderr, "Error: unsupported metric type\n");
            assert(0);
        }
        fout.print("assign f = ({} > ref_err);\n", err);
        fout.print("endmodule\n");
        fout.flush();
        // convert the Verilog file to BLIF using yosys
        string yosysComm = fmt::format("yosys -p \"read_verilog {}; synth; write_blif {}\" > {}", devCompNetPath + ".v", devCompNetPath + ".blif", devCompNetPath + "_yosys.log"); 
        fmt::print("Converting the Verilog file to BLIF\n");
        ExecSystComm(yosysComm);
        // convert the BLIF file to AIG
        AbcMan abc;
        string abcComm = fmt::format("r {}; st; resyn2rs; resyn2rs; resyn2rs; w {}", devCompNetPath + ".blif", devCompNetPath + "_opt.blif");
        fmt::print("Optimize the BLIF file\n");
        abc.Comm(abcComm, true);
    }

    // load the deviation network
    auto pDevCompNet = std::make_shared<NetMan>(devCompNetPath + "_opt.blif");
    pDevCompNet->Comm("ps");
    return pDevCompNet;
}


/**
 * @brief Generate a deviation network to check whether the maximum error is within the upper bound
 * 
 * @param pDevNet                the deviation network to measure the maximum error
 * @param outWidth               output width of the accurate and approximate networks
 * @param errUppBound            upper bound of the maximum error
 * @return pRetNet               the deviation network to check whether the maximum error is within the upper bound
 */
std::shared_ptr<NetMan> GenDevCompNetEmbedErrBound(std::shared_ptr<NetMan> pDevNet, int outWidth, ll errUppBound) {
    // copy the deviation network
    auto pRetNet = std::make_shared<NetMan>(*pDevNet);
    pRetNet->RenameNet(pDevNet->GetNetName() + "_embed_err_bound");
    // check the "ref_err" PIs
    int refErrWidth = -1;
    if (pDevNet->GetNetName() == "error_distance_devcomp")
        refErrWidth = outWidth;
    else if (pDevNet->GetNetName() == "hamming_distance_devcomp")
        refErrWidth = static_cast<int>(log2(outWidth)) + 1;
    else {
        fmt::print(stderr, "Error: unknown error metric\n");
        assert(0);
    }
    AbcObjVect pRefErrs;
    pRefErrs.reserve(refErrWidth);
    assert(pRetNet->GetPiNum() == outWidth * 2 + refErrWidth); // 0~outWidth-1: accOut; outWidth~2*outWidth-1: appOut; 2*outWidth~refErrWidth-1: ref_err
    for (int i = 0; i < refErrWidth; ++i) {
        auto pRefErrI = pRetNet->GetPi(outWidth * 2 + i);
        pRefErrs.emplace_back(pRefErrI);
        assert(pRetNet->GetName(pRefErrI) == fmt::format("ref_err[{}]", i));
    }
    // replace the "ref_err" input by the error upper bound
    auto constIds = pRetNet->CreateConstsIfNotExist();
    auto refErr = errUppBound;
    for (int i = 0; i < refErrWidth; ++i) {
        auto pRefErrIId = pRetNet->GetPiId(outWidth * 2 + i);
        if (refErr & 1) // replace by 1 
            pRetNet->TransfFanout(pRefErrIId, constIds.second); 
        else // replace by 0
            pRetNet->TransfFanout(pRefErrIId, constIds.first);
        refErr >>= 1;
    }
    assert(refErr == 0);
    // remove the "ref_err" PIs
    for (auto pRefErrI: pRefErrs)
        pRetNet->DelObj(pRefErrI);
    assert(pRetNet->GetPiNum() == outWidth * 2);
    // optimize the deviation network with embedded error upper bound
    fmt::print("{}\n", DASH_LINE);
    fmt::print("Optimizing the deviation network with embedded error upper bound {}\n", errUppBound);
    pRetNet->Comm("st; resyn2rs; resyn2rs; resyn2rs; logic; sop; ps;");
    fmt::print("{}\n", DASH_LINE);
    return pRetNet;
}

