/**
 * @file als.cc
 * @author Chang Meng (chang.meng@epfl.ch)
 * @brief Approximate logic synthesis
 * @date 2025-01-21
 * 
 */
#include "als.h"
#include "pbd.h"


/**
 * @brief Run the ALS algorithm, DATE'23 version
 * @brief In this version, the LAC type is fixed to SASIMI
 * 
 * @return void 
 */
void ALSMan::Run_v1() {
    fmt::print("Use MECALS 1.0 (DATE'23 version)\n");
    auto accNetStrash = accNet;
    accNetStrash.Comm("st");
    pDevCompNetEmbErr->Comm("st");
    auto appNet = accNetStrash;
    int round = 0;
    std::ostringstream oss("");
    std::string netName = std::string(accNetStrash.GetNet()->pName);
    oss << options.outpPath << round;
    appNet.WriteNet(oss.str() + ".blif", true); // appNet.WriteNet(oss.str() + ".dot", true);
    clock_t st = clock();
    for (round = 1; ; ++round) {
        std::cout << "---------- round " << round <<  "---------- " << std::endl;
        mecals_v1::PBDMan pbdMan;
        pbdMan.BuildMit(accNetStrash, appNet, *pDevCompNetEmbErr);
        pbdMan.BuildPBD(options.mecals1_exactPBDPerc);
        int ret = pbdMan.Synth(true); // enable SASIMI
        if (ret == -1)
            break;
        appNet = pbdMan.PostProc();
        // appNet.Print(true);
        std::ostringstream oss("");
        oss << options.outpPath << round;
        appNet.WriteNet(oss.str() + ".blif", true); // appNet.WriteNet(oss.str() + ".dot", true);
        std::cout << "current runtime = " << (clock() - st) / 1000000.0 << "s" << std::endl;
    }
    oss.str("");
    oss << options.outpPath << "final";
    appNet.WriteNet(oss.str() + ".blif", true);
    std::cout << "current runtime = " << (clock() - st) / 1000000.0 << "s" << std::endl;
    appNet.Comm("dch; amap; stime;");
}


/**
 * @brief Run the ALS algorithm, TCAD'25 version
 * 
 * @return void 
 */
void ALSMan::Run_v2() {
    // print exact circuit information
    fmt::print("Accurate circuit: size = {}, depth = {}, #PI = {}, #PO = {}\n", accNet.GetArea(), accNet.GetDelay(), accNet.GetPiNum(), accNet.GetPoNum());
    accNet.WriteBlif(fmt::format("{}r0_{}0_s{}_d{}.blif", options.outpPath, options.metrType, accNet.GetArea(), accNet.GetDelay()));

    // prepare the approximate circuit, counter example index, and LAC black list
    auto appNet = accNet;
    if (appNet.GetNetType() != NET_TYPE::SOP) {
        fmt::print(stderr, "Error: the network should be in SOP form\n");
        assert(0);
    }

    // main loop
    auto startTime = std::chrono::high_resolution_clock::now();
    int round = 1;
    bool fInclConst = 1; // when generating SASIMI LACs, include constant LACs
    bool fSimplify = 0;
    SimplifyWithSingleLac(LAC_TYPE::SASIMI, appNet, round, startTime, fInclConst, fSimplify);

    // accurate logic synthesis
    int resyn2rsTime = 3;
    std::string resyn2rsStr;
    for (int i = 0; i < resyn2rsTime; ++i)
        resyn2rsStr += "resyn2rs; ps;";
    appNet.Comm(fmt::format("ps; st; {}", resyn2rsStr));
    appNet.Comm("dch; amap; stime;");
    double area = appNet.GetArea(), delay = appNet.GetDelay();
    appNet.WriteNet(fmt::format("{}final_mapped_a{:.2f}_d{:.2f}.blif", options.outpPath, area, delay), true);
    PrintRuntime(startTime);
}


/**
 * @brief Run truncation-based ALS 
 * 
 * @return void 
 */
void ALSMan::Run_Trunc() {
    assert(options.metrType == METR_TYPE::MAXED);
    assert(accNet.GetNetType() == NET_TYPE::SOP);
    auto appNet = accNet;
    // create constant nodes
    auto constIds = appNet.CreateConstsIfNotExist();
    appNet.MergeConst();
    IntVect replTrace;
    for (int iBit = 0; iBit < appNet.GetPoNum(); ++iBit) {
        // set the iBit-th PO as 0
        auto drivId = appNet.GetPoDrivId(iBit);
        // replace the PO driver with a constant 0 node
        if (drivId != constIds.first)
            appNet.TempRepl_v2(drivId, constIds.first, replTrace, true);
        // measure the error
        ErrMan errMan(accNet, appNet, *pDevCompNetEmbErr);
        auto res = errMan.SolveSat(true);
        // check
        assert(res != CMSat::l_Undef);
        if (res == CMSat::l_True) {
            appNet.Recov_v2(replTrace, true);
            break;
        }
    }
    // accurate logic synthesis
    int resyn2rsTime = 3;
    std::string resyn2rsStr;
    for (int i = 0; i < resyn2rsTime; ++i)
        resyn2rsStr += "resyn2rs; ps;";
    appNet.Comm(fmt::format("ps; st; {}", resyn2rsStr));
    appNet.Comm("dch; amap; stime;");
    double area = appNet.GetArea(), delay = appNet.GetDelay();
    appNet.WriteNet(fmt::format("{}final_mapped_a{:.2f}_d{:.2f}.blif", options.outpPath, area, delay), true);
}


/**
 * @brief ALS using mixed LACs
 * 
 * @return void 
 */
void ALSMan::Run_FastFlow() {
    fmt::print("Use fast flow\n");
    auto startTime = std::chrono::high_resolution_clock::now();
    assert(accNet.GetNetType() == NET_TYPE::SOP);
    auto appNet = accNet;
    // truncation
    auto constIds = appNet.CreateConstsIfNotExist();
    appNet.MergeConst();
    IntVect replTrace;
    if (options.metrType == METR_TYPE::MAXED) {
    for (int iBit = 0; iBit < appNet.GetPoNum(); ++iBit) {
        auto drivId = appNet.GetPoDrivId(iBit);
        // try constant 0
        if (drivId != constIds.first)
            appNet.TempRepl_v2(drivId, constIds.first, replTrace, false);
        ErrMan errMan0(accNet, appNet);
        auto maxErr0 = errMan0.ComputeMaxErr(METR_TYPE::MAXED);
        appNet.Recov_v2(replTrace, false);
        // try constant 1
        if (drivId != constIds.second)
            appNet.TempRepl_v2(drivId, constIds.second, replTrace, false);
        ErrMan errMan1(accNet, appNet);
        auto maxErr1 = errMan1.ComputeMaxErr(METR_TYPE::MAXED);
        appNet.Recov_v2(replTrace, false);
        // apply the best one
        if (maxErr0 <= options.errUppBound || maxErr1 <= options.errUppBound) {
            if (maxErr0 <= maxErr1) {
                if (drivId != constIds.first)
                    appNet.TempRepl_v2(drivId, constIds.first, replTrace, true);
                fmt::print("current error = {}\n", maxErr0.str());
            }
            else {
                if (drivId != constIds.second)
                    appNet.TempRepl_v2(drivId, constIds.second, replTrace, true);
                fmt::print("current error = {}\n", maxErr1.str());
            }
        }
        else
            break;
    }
    }
    appNet.Sweep();
    appNet.WriteBlif(fmt::format("{}r0_{}xxx_s{}_d{}.blif", options.outpPath, options.metrType, appNet.GetArea(), appNet.GetDelay()));
    PrintRuntime(startTime);
    
    // try constant LAC
    int round = 1;
    bool fInclConst = 0;
    bool fSimplify = 1;
    SimplifyWithSingleLac(LAC_TYPE::CONSTANT, appNet, round, startTime, fInclConst, fSimplify);

    // try SASIMI LAC
    SimplifyWithSingleLac(LAC_TYPE::SASIMI, appNet, round, startTime, fInclConst, fSimplify);

    // accurate logic synthesis
    int resyn2rsTime = 3;
    std::string resyn2rsStr;
    for (int i = 0; i < resyn2rsTime; ++i)
        resyn2rsStr += "resyn2rs; ps;";
    appNet.Comm(fmt::format("ps; st; {}", resyn2rsStr));
    appNet.Comm("dch; amap; stime;");
    double area = appNet.GetArea(), delay = appNet.GetDelay();
    appNet.WriteNet(fmt::format("{}final_mapped_a{:.2f}_d{:.2f}.blif", options.outpPath, area, delay), true);
    PrintRuntime(startTime);
}


/**
 * @brief Simplify the network using a certain type of LAC
 * 
 * @param lacType 
 * @param appNet 
 * @param startTime 
 */
void ALSMan::SimplifyWithSingleLac(LAC_TYPE lacType, NetMan& appNet, int& round, const std::chrono::time_point<std::chrono::high_resolution_clock>& startTime, bool inclConst, bool fSimplify) {
    // zero error upper bound
    if (options.errUppBound == 0) {
        fmt::print("Early stop: the maximum error upper bound is 0\n");
        return;
    }

    // prepare
    fmt::print("{}\n", DASH_LINE);
    fmt::print("{}Using {} LAC{}\n", HALF_DASH_LINE, lacType, HALF_DASH_LINE);
    fmt::print("{}\n", DASH_LINE);
    static int countExNum = 0; // track the ending index of counter-example simulation patterns in accSmlt
    std::unordered_set<std::string> lacBlackList; // black list for LACs, which cause the SAT solver to return UNDEF
    int oldSize = static_cast<int>(appNet.GetArea()), oldDepth = static_cast<int>(appNet.GetDelay());

    // main loop
    for (; ; ++round) {
        fmt::print("{}round {}{}\n", HALF_DASH_LINE, round, HALF_DASH_LINE);
        // generate LACs
        LACMan lacMan;
        if (lacType == LAC_TYPE::CONSTANT)
            lacMan.GenConstLACs(appNet);
        else if (lacType == LAC_TYPE::SASIMI)
            lacMan.GenSasimiLACs(appNet, options.maxCandLacs, inclConst);
        else if (lacType == LAC_TYPE::RESUB)
            lacMan.GenResubLACs(appNet, options.seed, options.appResub_nFrame4ResubGen, options.appResub_maxLevelDiff, options.maxCandLacs, inclConst);
        else {
            fmt::print(stderr, "Error: unsupported LAC type\n");
            assert(0);
        }
        // remove LACs in the black list
        lacMan.RemLacsFromBlackList(lacBlackList);
        // use logic simulation to estimate maximum error lower bound and prune large-error LACs
        BatchErrEst errEst(options.metrType, options.seed, options.nFrame);
        errEst.CompLacErrsBySimAndPruneBadLacs(lacMan, *pAccSmlt, appNet, options.errUppBound);
        // if no valid LACs, break
        if (lacMan.GetLacNum() == 0)
            break;
        // sort LACs and keep top-K
        const int TOP_K_LAC = 100;
        lacMan.SortAndKeepTopKLACs(TOP_K_LAC);
        lacMan.PrintLACs(10);
        // use SAT to estimate error
        bool existValidLac = ApplyMultValidLacs(lacMan, appNet, lacBlackList, countExNum);
        // bool existValidLac = ApplyMultValidLacs_NoSimPrune(lacMan, appNet, lacBlackList, countExNum);
        if (!existValidLac)
            break;
        // post processing
        int size = static_cast<int>(appNet.GetArea()), depth = static_cast<int>(appNet.GetDelay());
        fmt::print("size = {}, depth = {}\n", size, depth);
        appNet.WriteBlif(fmt::format("{}r{}_{}xxx_s{}_d{}.blif", options.outpPath, round, options.metrType, size, depth));
        PrintRuntime(startTime);
        assert(size < oldSize || (size == oldSize && depth <= oldDepth));
        if (size == oldSize && depth == oldDepth) { // early stop due to no improvement
            fmt::print("Early stop: no size or depth improvement\n");
            break;
        }
        else {
            oldSize = size;
            oldDepth = depth;
        }
    }
    PrintRuntime(startTime);

    // final logic synthesis
    if (fSimplify)
        appNet.Comm("st; resyn2rs; ps; resyn2rs; ps; resyn2rs; ps; logic; sop;");
}


/**
 * @brief Apply multiple LACs, after which the real maximum error is no more than the given bound
 * @brief Assume that the LACs are sorted: primary key = (smaller) error, secondary key = (larger) sizeGain
 * 
 * @param  lacMan         the LAC manager
 * @param  accSmlt        the accurate network's simulator
 * @param  appNet         the approximate network
 * @param  devNet         the deviation network
 * @param  lacBlackList   the black list of LACs, which cause the SAT solver to return undefined
 * @param  countExNum     #counter examples
 * @retval accSmlt        the accurate network's simulator that includes the counter examples
 * @retval appNet         the approximate network after applying the LACs
 * @retval lacBlackList   the updated black list of LACs
 * @retval countExNum     the updated #counter examples
 * @return true if there exists at least one valid LAC; false otherwise
 */
bool ALSMan::ApplyMultValidLacs(const LACMan& lacMan, NetMan& appNet, std::unordered_set<std::string>& lacBlackList, int& countExNum) {
    auto startTime = std::chrono::high_resolution_clock::now();
    fmt::print("Check the maximum error for each LAC using SAT and apply multiple LACs\n");
    // prepare the counter example
    assert(pAccSmlt != nullptr);
    const auto& accNet = static_cast<const NetMan&>(*pAccSmlt);
    IntVect counterEx;
    int nPi = accNet.GetPiNum();
    assert(nPi > 0);
    counterEx.reserve(nPi);
    std::vector<BitVect> counterExPiPatts(nPi); // piPatts[i][j] is the i-th PI's j-th counter-example pattern
    for (auto& piPatt: counterExPiPatts)
        piPatt.resize(0);
    // apply the LACs in a heuristic way
    IntVect replTrace;
    IntSet frozTargNodes;
    bool existValidLac = false;
    // int nAppliedLac = 0;
    for (int iLacId = 0; iLacId < lacMan.GetLacNum(); ++iLacId) {
        auto pLac = lacMan.GetLac(iLacId);
        fmt::print("{}checking {}-th LAC: {}\n", HALF_DASH_LINE, iLacId, pLac->ToStr());
        // skip the LAC if the target node is frozen
        int targId = pLac->GetTargId();
        if (frozTargNodes.count(targId)) {
            fmt::print("Warning: the target node is frozen, skip this LAC\n");
            continue;
        }
        // check whether the target node is active or not
        if (appNet.GetFanoutNum(targId) == 0) {
            fmt::print("The target node is dangling, skip this LAC\n");
            continue;
        }
        // temporarily apply the LAC
        TempApplyLac(appNet, *pLac, replTrace, false);
        // if the network is cyclic, skip this LAC
        if (!appNet.IsAcyclic()) {
            fmt::print("Warning: the network is cyclic, skip this LAC");
            RecovNet(appNet, {replTrace}, false);
            continue;
        }
        // create an error manager
        assert(pDevCompNetEmbErr != nullptr);
        ErrMan errMan(accNet, appNet, *pDevCompNetEmbErr);
        const auto& errMit = errMan.GetErrMit();
        assert(ComparePi(accNet, errMit, false));
        // fast checking using counter examples
        if (counterExPiPatts[0].size()) {
            Simulator errMitSmlt(errMit, 0, counterExPiPatts[0].size(), DISTR_TYPE::UNIF);
            errMitSmlt.GenInpFromBitVects(counterExPiPatts);
            errMitSmlt.UpdNodeAndPoPatts();
            assert(errMit.GetPoNum() == 1);
            auto& outDat = errMitSmlt.GetDat(errMitSmlt.GetPoId(0));
            if (outDat.count()) { // if the only PO is 1, then the error constraint is not satisfied
                fmt::print("Fast checking: Exceed the error bound, skip this LAC\n");
                RecovNet(appNet, {replTrace}, false);
                continue;
            }
        }
        // solve the SAT problem
        auto res = errMan.SolveSat(counterEx, true);
        if (res == CMSat::l_False) { // UNSAT, satisfy the error bound, apply the LAC
            fmt::print("Satisfy the error bound, apply the LAC\n");
            existValidLac = true;
            // freeze nodes
            frozTargNodes.insert(targId);
        }
        else if (res == CMSat::l_True) { // SAT, exceed the error bound, skip the LAC
            fmt::print("Exceed the error bound, save the {}-th counter example\n", countExNum);
            // save the counter example in the countExNum-th PI pattern in the accSmlt
            pAccSmlt->ReplInp(countExNum, counterEx);
            ++countExNum;
            if (countExNum >= options.nFrame)
                countExNum = 0;
            for (int iPi = 0; iPi < nPi; ++iPi) {
                assert(counterEx[iPi] == 0 || counterEx[iPi] == 1);
                counterExPiPatts[iPi].push_back(counterEx[iPi]);
            }
            // recover the network
            RecovNet(appNet, {replTrace}, false);
        }
        else { // UNDEF, skip the LAC
            fmt::print("Warning: SAT solver returns undefined, skip this LAC and add it to the black list\n");
            // add the LAC to the black list
            lacBlackList.insert(pLac->ToStrShort());
            // recover the network
            RecovNet(appNet, {replTrace}, false);
        }
    }
    // re-simulate the accurate network if counter examples are generated
    if (counterExPiPatts[0].size())
        pAccSmlt->UpdNodeAndPoPatts();
    // clean up the appNet
    appNet.Sweep(false);
    assert(appNet.Check());
    PrintRuntime(startTime, "apply multiple LACs");
    return existValidLac;
}


bool ALSMan::ApplyMultValidLacs_NoSimPrune(const LACMan& lacMan, NetMan& appNet, std::unordered_set<std::string>& lacBlackList, int& countExNum) {
    auto startTime = std::chrono::high_resolution_clock::now();
    fmt::print("Check the maximum error for each LAC using SAT and apply multiple LACs\n");
    // prepare the counter example
    assert(pAccSmlt != nullptr);
    const auto& accNet = static_cast<const NetMan&>(*pAccSmlt);
    IntVect counterEx;
    int nPi = accNet.GetPiNum();
    assert(nPi > 0);
    counterEx.reserve(nPi);
    std::vector<BitVect> counterExPiPatts(nPi); // piPatts[i][j] is the i-th PI's j-th counter-example pattern
    for (auto& piPatt: counterExPiPatts)
        piPatt.resize(0);
    // apply the LACs in a heuristic way
    IntVect replTrace;
    IntSet frozTargNodes;
    // bool existValidLac = false;
    int nAppliedLac = 0;
    for (int iLacId = 0; iLacId < lacMan.GetLacNum(); ++iLacId) {
        auto pLac = lacMan.GetLac(iLacId);
        fmt::print("{}checking {}-th LAC: {}\n", HALF_DASH_LINE, iLacId, pLac->ToStr());
        // skip the LAC if the target node is frozen
        int targId = pLac->GetTargId();
        if (frozTargNodes.count(targId)) {
            fmt::print("Warning: the target node is frozen, skip this LAC\n");
            continue;
        }
        // check whether the target node is active or not
        if (appNet.GetFanoutNum(targId) == 0) {
            fmt::print("The target node is dangling, skip this LAC\n");
            continue;
        }
        // temporarily apply the LAC
        TempApplyLac(appNet, *pLac, replTrace, false);
        // if the network is cyclic, skip this LAC
        if (!appNet.IsAcyclic()) {
            fmt::print("Warning: the network is cyclic, skip this LAC");
            RecovNet(appNet, {replTrace}, false);
            continue;
        }
        // create an error manager
        assert(pDevCompNetEmbErr != nullptr);
        ErrMan errMan(accNet, appNet, *pDevCompNetEmbErr);
        const auto& errMit = errMan.GetErrMit();
        assert(ComparePi(accNet, errMit, false));
        // fast checking using counter examples
        if (counterExPiPatts[0].size()) {
            Simulator errMitSmlt(errMit, 0, counterExPiPatts[0].size(), DISTR_TYPE::UNIF);
            errMitSmlt.GenInpFromBitVects(counterExPiPatts);
            errMitSmlt.UpdNodeAndPoPatts();
            assert(errMit.GetPoNum() == 1);
            auto& outDat = errMitSmlt.GetDat(errMitSmlt.GetPoId(0));
            if (outDat.count()) { // if the only PO is 1, then the error constraint is not satisfied
                fmt::print("Fast checking: Exceed the error bound, skip this LAC\n");
                RecovNet(appNet, {replTrace}, false);
                continue;
            }
        }
        // solve the SAT problem
        auto res = errMan.SolveSat(counterEx, true);
        if (res == CMSat::l_False) { // UNSAT, satisfy the error bound, apply the LAC
            fmt::print("Satisfy the error bound, apply the LAC\n");
            // existValidLac = true;
            ++nAppliedLac;
            const int MAX_APPLY_NUM = 100;
            if (nAppliedLac >= MAX_APPLY_NUM)
                break;
            // freeze nodes
            frozTargNodes.insert(targId);
        }
        else if (res == CMSat::l_True) { // SAT, exceed the error bound, skip the LAC
            fmt::print("Exceed the error bound, save the {}-th counter example\n", countExNum);
            // save the counter example in the countExNum-th PI pattern in the accSmlt
            pAccSmlt->ReplInp(countExNum, counterEx);
            ++countExNum;
            if (countExNum >= options.nFrame)
                countExNum = 0;
            for (int iPi = 0; iPi < nPi; ++iPi) {
                assert(counterEx[iPi] == 0 || counterEx[iPi] == 1);
                counterExPiPatts[iPi].push_back(counterEx[iPi]);
            }
            // recover the network
            RecovNet(appNet, {replTrace}, false);
        }
        else { // UNDEF, skip the LAC
            fmt::print("Warning: SAT solver returns undefined, skip this LAC and add it to the black list\n");
            // add the LAC to the black list
            lacBlackList.insert(pLac->ToStrShort());
            // recover the network
            RecovNet(appNet, {replTrace}, false);
        }
    }
    // re-simulate the accurate network if counter examples are generated
    if (counterExPiPatts[0].size())
        pAccSmlt->UpdNodeAndPoPatts();
    // clean up the appNet
    appNet.Sweep(false);
    assert(appNet.Check());
    PrintRuntime(startTime, "apply multiple LACs");
    // return existValidLac;
    return (nAppliedLac > 0);
}


// /**
//  * @brief Apply multiple LACs using incremental SAT, after which the real maximum error is no more than the given bound
//  * @brief Assume that the LACs are sorted: primary key = (smaller) error, secondary key = (larger) sizeGain
//  * 
//  * @param  lacMan         the LAC manager
//  * @param  appNet         the approximate network
//  * @param  lacBlackList   the black list of LACs, which cause the SAT solver to return undefined
//  * @param  countExNum     #counter examples
//  * @retval pAccSmlt       the accurate network's simulator that includes the counter examples
//  * @retval appNet         the approximate network after applying the LACs
//  * @retval lacBlackList   the updated black list of LACs
//  * @retval countExNum     the updated #counter examples
//  * @return true if there exists at least one valid LAC; false otherwise
//  */
// bool ALSMan::ApplyMultValidLacs_IncSAT(const LACMan& lacMan, NetMan& appNet, std::unordered_set<std::string>& lacBlackList, int& counterExNum) {
//     auto startTime = std::chrono::high_resolution_clock::now();
//     fmt::print("Check the maximum error for each LAC using SAT and apply multiple LACs\n");

//     // prepare counter example patterns
//     assert(pAccSmlt != nullptr);
//     const auto& accNet = static_cast<const NetMan&>(*pAccSmlt);
//     int accPiNum = accNet.GetPiNum();
//     assert(accPiNum > 0);
//     std::vector<BitVect> counterExPiPatts(accPiNum); // piPatts[i][j] is the i-th PI's j-th counter-example pattern
//     for (auto& piPatt: counterExPiPatts)
//         piPatt.resize(0);

//     // incremental SAT solving for each group
//     const int GROUP_SIZE = 20;
//     bool fUpd = false;
//     IntSet appearedNodes;
//     LACPtrVect candidateLacs;
//     candidateLacs.reserve(GROUP_SIZE);
//     int iGroup = 0, procLacIdx = 0;
// while (1) {
//     fmt::print("{}group {}{}\n", HALF_DASH_LINE, iGroup++, HALF_DASH_LINE);
//     // make sure that the candidate LACs are applied to different nodes
//     candidateLacs.clear();
//     while (static_cast<int>(candidateLacs.size()) < GROUP_SIZE && procLacIdx < lacMan.GetLacNum()) {
//         auto pLac = lacMan.GetLac(procLacIdx++);
//         int nodeId = pLac->GetTargId();
//         if (appearedNodes.count(nodeId))
//             fmt::print("ignore LAC: {}\n", pLac->ToStr());
//         else {
//             appearedNodes.insert(nodeId);
//             candidateLacs.emplace_back(pLac);
//         }
//     }
//     fmt::print("#candidate LACs: {}\n", candidateLacs.size());
//     if (candidateLacs.empty())
//         break;

//     // integrate candidate LACs into the approximate network
//     Int2DVect replTraces;
//     int nCandLacs = candidateLacs.size();
//     replTraces.resize(nCandLacs);
//     for (int iLac = 0; iLac < nCandLacs; ++iLac) {
//         auto pLac = candidateLacs[iLac];
//         TempApplyLacWithController(appNet, *pLac, iLac, replTraces[iLac], false);
//     }
//     assert(appNet.IsAcyclic());

//     // build the error miter
//     assert(pDevCompNetEmbErr != nullptr);
//     ErrMan errMan(accNet, appNet, *pDevCompNetEmbErr);
//     const auto& errMit = errMan.GetErrMit();
//     assert(ComparePi(accNet, errMit, true));
//     assert(errMit.GetPiNum() == accPiNum + nCandLacs); // errMit PIs = accNet PIs + control signals (#control signals = #candidate LACs)
//     assert(errMit.GetPoNum() == 1);

//     // prepare the counter example simulator
//     Simulator errMitSmlt(errMit, 0, counterExPiPatts[0].size(), DISTR_TYPE::UNIF);
//     errMitSmlt.GenInpFromBitVects(counterExPiPatts);

//     // incremental SAT solving
//     std::vector<CMSat::Lit> assumpts(nCandLacs);
//     IntVect applLacIds;
//     IntVect counterEx;
//     for (int iLac = 0; iLac < nCandLacs; ++iLac) {
//         auto pLac = candidateLacs[iLac];
//         fmt::print("{}checking {}-th LAC: {}\n", HALF_DASH_LINE, iLac, pLac->ToStr());
//         // fast checking using simulation
//         if (errMitSmlt.GetFrameNumb() > 0) {
//             for (int iCtrl = 0; iCtrl < nCandLacs; ++iCtrl) // reset the control signals
//                 errMitSmlt.SetPiConst(accPiNum + iCtrl, 0);
//             errMitSmlt.SetPiConst(accPiNum + iLac, 1); // the i-th LAC is applied
//             for (auto applLacId: applLacIds) // already applied LACs
//                 errMitSmlt.SetPiConst(accPiNum + applLacId, 1);
//             errMitSmlt.UpdNodeAndPoPatts();
//             auto simRes = (errMitSmlt.GetDat(errMitSmlt.GetPoId(0))).count();
//             if (simRes) { // if the only PO is 1, then the error constraint is not satisfied
//                 fmt::print("Fast checking: Exceed the error bound, skip this LAC\n");
//                 errMan.AddUnitClauseOfPi(accPiNum + iLac, true); // add unit clause for the i-th LAC (ctrl_i = 0)
//                 continue;
//             }
//         }
//         // here, LAC 0~iLac-1 have already been processed by adding unit clauses
//         // we only need to add assumptions for LAC iLac~nCandLacs-1
//         assumpts.clear();
//         assumpts.emplace_back(errMan.GetCnfVarIdOfIthPi(accPiNum + iLac), false); // ctrl_{iLac} = 1
//         for (int iCtrl = iLac + 1; iCtrl < nCandLacs; ++iCtrl)
//             assumpts.emplace_back(errMan.GetCnfVarIdOfIthPi(accPiNum + iCtrl), true); // ctrl_{iCtrl} = 0
//         // solve the SAT problem
//         auto res = errMan.SolveSat(assumpts, counterEx, true);
//         if (res == CMSat::l_False) { // UNSAT, satisfy the error bound
//             fmt::print("Satisfy the error bound, apply the LAC\n");
//             applLacIds.emplace_back(iLac);
//             // add unit clause for the i-th LAC (ctrl_i = 1)
//             errMan.AddUnitClauseOfPi(accPiNum + iLac, false);
//         }
//         else if (res == CMSat::l_True) { // SAT, exceed the error bound
//             fmt::print("Exceed the error bound, save the {}-th counter example\n", counterExNum);
//             // save the counter example in the countExNum-th PI pattern in the accSmlt
//             assert(static_cast<int>(counterEx.size()) == accPiNum + nCandLacs);
//             errMitSmlt.AppendInp(counterEx);
//             counterEx.resize(accPiNum);
//             pAccSmlt->ReplInp(counterExNum, counterEx);
//             ++counterExNum;
//             if (counterExNum >= options.nFrame)
//                 counterExNum = 0;
//             for (int iPi = 0; iPi < accPiNum; ++iPi) {
//                 assert(counterEx[iPi] == 0 || counterEx[iPi] == 1);
//                 counterExPiPatts[iPi].push_back(counterEx[iPi]);
//             }
//             // add unit clause for the i-th LAC (ctrl_i = 0)
//             errMan.AddUnitClauseOfPi(accPiNum + iLac, true);
//         }
//         else {
//             fmt::print("Warning: SAT solver returns undefined, skip this LAC and add it to the black list\n");
//             // add the LAC to the black list
//             lacBlackList.insert(pLac->ToStrShort());
//             // add unit clause for the i-th LAC (ctrl_i = 0)
//             errMan.AddUnitClauseOfPi(accPiNum + iLac, true);
//         }
//     }

//     // recover the appNet
//     RecovNet(appNet, replTraces, false);
    
//     // apply the valid LACs
//     fmt::print("apply the valid LACs\n");
//     IntVect tmpReplTrace;
//     for (auto appliLacId: applLacIds) {
//         auto pLac = candidateLacs[appliLacId];
//         fmt::print("apply the {}-th LAC: {}\n", appliLacId, pLac->ToStr());
//         TempApplyLac(appNet, *pLac, tmpReplTrace, false);
//     }
//     assert(appNet.Check());
//     fUpd |= (!applLacIds.empty());
// }

//     // sweep the appNet
//     appNet.Sweep(false);

//     // re-simulate the accurate network because counter examples may be generated
//     pAccSmlt->UpdNodeAndPoPatts();
 
//     PrintRuntime(startTime, "apply multiple LACs");
//     return fUpd;
// }


// /**
//  * @brief Apply multiple LACs, after which the real maximum error is no more than the given bound
//  * @brief Assume that the LACs are sorted: primary key = (smaller) error, secondary key = (larger) sizeGain
//  * 
//  * @param  lacMan         the LAC manager
//  * @param  appNet         the approximate network
//  * @param  lacBlackList   the black list of LACs, which cause the SAT solver to return undefined
//  * @param  countExNum     #counter examples
//  * @retval accSmlt        the accurate network's simulator that includes the counter examples
//  * @retval appNet         the approximate network after applying the LACs
//  * @retval lacBlackList   the updated black list of LACs
//  * @retval countExNum     the updated #counter examples
//  * @return true if there exists at least one valid LAC; false otherwise
//  */
// bool ALSMan::ApplyMultValidLacsUsingBaseErr(const LACMan& lacMan, NetMan& appNet, std::unordered_set<std::string>& lacBlackList, int& countExNum) {
//     auto startTime = std::chrono::high_resolution_clock::now();
//     fmt::print("Apply multiple LACs using base error\n");
//     // get real maximum error
//     ErrMan errManReal(accNet, appNet);
//     auto baseErr = errManReal.ComputeMaxErr(options.metrType);
//     auto errMarg = BigInt(options.errUppBound) - baseErr;
//     fmt::print("Base error = {}, error margin = {}\n", baseErr.str(), errMarg.str());
//     assert(errMarg >= 0 && errMarg <= LL_MAX);
//     // if (errMarg == 0) {
//     //     fmt::print("Error margin is 0, no need to apply LACs\n");
//     //     return false;
//     // }
//     // copy appNet & prepare deviation+comparator network
//     auto appNetCopy = appNet;
//     auto pCurrDevCompNet = GenDevCompNetEmbedErrBound(pDevCompNet, accNet.GetPoNum(), static_cast<ll>(errMarg));
//     // prepare the counter example
//     assert(pAccSmlt != nullptr);
//     int nPi = pAccSmlt->GetPiNum();
//     assert(nPi > 0);
//     IntVect counterEx;
//     counterEx.reserve(nPi);
//     std::vector<BitVect> counterExPiPatts(nPi); // piPatts[i][j] is the i-th PI's j-th counter-example pattern
//     for (auto& piPatt: counterExPiPatts)
//         piPatt.resize(0);
//     // iteratively apply the LACs (sorted by error upper bound) until exceeding the error bound
//     IntVect replTrace;
//     IntSet frozTargNodes;
//     bool existValidLac = false;
//     for (int iLacId = 0; iLacId < lacMan.GetLacNum(); ++iLacId) {
//         auto pLac = lacMan.GetLac(iLacId);
//         fmt::print("{}checking {}-th LAC: {}\n", HALF_DASH_LINE, iLacId, pLac->ToStr());
//         // skip the LAC if the target node is frozen
//         int targId = pLac->GetTargId();
//         if (frozTargNodes.count(targId)) {
//             fmt::print("Warning: the target node is frozen, skip this LAC\n");
//             continue;
//         }
//         // temporarily apply the LAC
//         TempApplyLac(appNet, *pLac, replTrace, false);
//         // if the network is cyclic, skip this LAC
//         if (!appNet.IsAcyclic()) {
//             fmt::print("Warning: the network is cyclic, skip this LAC");
//             RecovNet(appNet, {replTrace}, false);
//             continue;
//         }
//         // create an error manager
//         ErrMan errMan(appNetCopy, appNet, *pCurrDevCompNet);
//         const auto& errMit = errMan.GetErrMit();
//         assert(ComparePi(appNetCopy, errMit, false));
//         // fast checking using counter examples
//         if (counterExPiPatts[0].size()) {
//             Simulator errMitSmlt(errMit, 0, counterExPiPatts[0].size(), DISTR_TYPE::UNIF);
//             errMitSmlt.GenInpFromBitVects(counterExPiPatts);
//             errMitSmlt.UpdNodeAndPoPatts();
//             assert(errMit.GetPoNum() == 1);
//             auto& outDat = errMitSmlt.GetDat(errMitSmlt.GetPoId(0));
//             if (outDat.count()) { // if the only PO is 1, then the error constraint is not satisfied
//                 fmt::print("Fast checking: Exceed the error bound, skip this LAC\n");
//                 RecovNet(appNet, {replTrace}, false);
//                 continue;
//             }
//         }
//         // solve the SAT problem
//         auto res = errMan.SolveSat(counterEx, true);
//         if (res == CMSat::l_False) { // UNSAT, satisfy the error bound
//             fmt::print("Satisfy the error bound, apply the LAC\n");
//             existValidLac = true;
//             // freeze nodes
//             frozTargNodes.insert(targId);
//         }
//         else if (res == CMSat::l_True) { // SAT, exceed the error bound
//             fmt::print("Exceed the error bound, save the {}-th counter example\n", countExNum);
//             // save the counter example in the countExNum-th PI pattern in the accSmlt
//             pAccSmlt->ReplInp(countExNum, counterEx);
//             ++countExNum;
//             if (countExNum >= options.nFrame)
//                 countExNum = 0;
//             for (int iPi = 0; iPi < nPi; ++iPi) {
//                 assert(counterEx[iPi] == 0 || counterEx[iPi] == 1);
//                 counterExPiPatts[iPi].push_back(counterEx[iPi]);
//             }
//             // recover the network
//             RecovNet(appNet, {replTrace}, false);
//         }
//         else {
//             fmt::print("Warning: SAT solver returns undefined, skip this LAC and add it to the black list\n");
//             // add the LAC to the black list
//             lacBlackList.insert(pLac->ToStrShort());
//             // recover the network
//             RecovNet(appNet, {replTrace}, false);
//         }
//     }
//     // re-simulate the accurate network if counter examples are generated
//     if (counterExPiPatts[0].size())
//         pAccSmlt->UpdNodeAndPoPatts();
//     // clean up the appNet, check, and print runtime
//     appNet.Sweep(false);
//     assert(appNet.Check());
//     PrintRuntime(startTime, "apply multiple LACs");
//     // return 
//     return existValidLac;
// }