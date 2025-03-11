/**
 * @file my_abc.cc
 * @author Chang Meng (chang.meng@epfl.ch)
 * @brief Wrapper for ABC, the open-source logic synthesis and verification system
 * @date 2025-01-25
 * 
 */
#include "my_abc.h"


using namespace abc;
using namespace std;


void AbcMan::Comm(const string& cmd, bool fVerb) {
    if (fVerb)
        fmt::print("Execute abc command: {}\n", cmd);
    if (Cmd_CommandExecute(GetAbcFrame(), cmd.c_str())) {
        fmt::print(stderr, "run ABC command {} failed\n", cmd);
        assert(0);
    }
}


void AbcMan::ReadNet(const std::string& fileName) {
    assert(IsPathExist(fileName));
    Comm("r " + fileName);
}


void AbcMan::WriteNet(const std::string& fileName, bool fVerb) {
    Comm("w " + fileName, fVerb);
}


void AbcMan::ReadStandCell(const std::string& fileName) {
    assert(IsPathExist(fileName));
    Comm("r " + fileName);
}


void AbcMan::ConvToAig() {
    Comm("aig");
}


void AbcMan::ConvToGate() {
    Map(MAP_TYPE::SCL, ORIENT::AREA);
}


void AbcMan::ConvToSop() {
    if (GetNetType() == NET_TYPE::STRASH)
        Comm("logic;");
    Comm("sop");
}


void AbcMan::Strash() {
    Comm("st");
}


void AbcMan::PrintStat() {
    Comm("ps");
}


void AbcMan::TopoSort() {
    auto type = GetNetType();
    assert(type == NET_TYPE::AIG || type == NET_TYPE::SOP || type == NET_TYPE::GATE);
    Comm("topo");

    // fix twin nodes
    auto pNtk = GetNet();
    if (Abc_NtkHasMapping(pNtk)) {
        Abc_Ntk_t * pNtkNew; 
        AbcObj* pObj, * pFanin;
        int i, k;
        assert(pNtk != nullptr);
        // start the network
        pNtkNew = Abc_NtkStartFrom( pNtk, pNtk->ntkType, pNtk->ntkFunc );
        // copy the internal nodes
        assert(!Abc_NtkIsStrash(pNtk));
        // duplicate the nets and nodes (CIs/COs/latches already dupped)
        set <int> skip;
        Abc_NtkForEachObj( pNtk, pObj, i ) {
            if ( pObj->pCopy == NULL && skip.count(pObj->Id) == 0 ) {
                Abc_NtkDupObj(pNtkNew, pObj, Abc_NtkHasBlackbox(pNtk) && Abc_ObjIsNet(pObj));
                auto pTwin = GetTwinNode(pObj);
                if (pTwin != nullptr) {
                    Abc_NtkDupObj(pNtkNew, pTwin, Abc_NtkHasBlackbox(pNtk) && Abc_ObjIsNet(pTwin));
                    skip.insert(pTwin->Id);
                }
            }
        }
        // reconnect all objects (no need to transfer attributes on edges)
        Abc_NtkForEachObj( pNtk, pObj, i )
            if ( !Abc_ObjIsBox(pObj) && !Abc_ObjIsBo(pObj) )
                Abc_ObjForEachFanin( pObj, pFanin, k )
                    Abc_ObjAddFanin( pObj->pCopy, pFanin->pCopy );
        // duplicate the EXDC Ntk
        if ( pNtk->pExdc )
            pNtkNew->pExdc = Abc_NtkDup( pNtk->pExdc );
        if ( pNtk->pExcare )
            pNtkNew->pExcare = Abc_NtkDup( (Abc_Ntk_t *)pNtk->pExcare );
        // duplicate timing manager
        if ( pNtk->pManTime )
            Abc_NtkTimeInitialize( pNtkNew, pNtk );
        if ( pNtk->vPhases )
            Abc_NtkTransferPhases( pNtkNew, pNtk );
        if ( pNtk->pWLoadUsed )
            pNtkNew->pWLoadUsed = Abc_UtilStrsav( pNtk->pWLoadUsed );
        // check correctness
        if ( !Abc_NtkCheck( pNtkNew ) )
            fprintf( stdout, "Abc_NtkDup(): Network check has failed.\n" );
        pNtk->pCopy = pNtkNew;
        // return pNtkNew;
        SetMainNet(pNtkNew);
    }
}


void AbcMan::StatTimeAnal() {
    assert(GetNetType() == NET_TYPE::GATE);
    assert(GetAbcFrame()->pLibScl != nullptr);
    TopoSort();
    Comm("stime");
}


void AbcMan::Synth(ORIENT orient, bool fVerb) {
    Comm("st", fVerb);
    if (fVerb)
        cout << orient << "-oriented synthesis" << endl;
    const string areaComm = {"st; compress2rs"};
    // const string delayComm = {"st; resyn2"};
    const string delayComm2 = {"st; ifraig; resyn2"};
    double oldArea = GetArea();
    double oldDelay = GetDelay();
    bool isCont = true;
    while (isCont) {
        isCont = false;
        auto oldNtk = Abc_NtkDup(GetNet());
        if (orient == ORIENT::AREA)
            Comm(areaComm, fVerb);
        else if (orient == ORIENT::DELAY)
            Comm(delayComm2, fVerb);
        else
            assert(0);
        auto res = make_pair <double, double> (GetArea(), GetDelay());
        if (fVerb)
            PrintStat();
        double newArea = res.first;
        double newDelay = res.second;
        IMPR impr = UpdNet(oldArea, oldDelay, oldNtk, newArea, newDelay, orient);
        if (impr == IMPR::GOOD) {
            oldArea = newArea;
            oldDelay = newDelay;
            isCont = true;
        }
        if (fVerb)
            cout << (impr == IMPR::GOOD? "accept": "cancel") << endl;
    }
    if (fVerb)
        PrintStat();
}


void AbcMan::SynthAndMap(double maxDelay, bool fVerb) {
    bool cont = true;
    if (fVerb)
        cout << "maxDelay = " << maxDelay << endl;
    TopoSort();
    while (cont) {
        double oldArea = numeric_limits <double>::max(), oldDelay = numeric_limits <double>::max();
        if (GetNetType(GetNet()) == NET_TYPE::GATE)
            oldArea = GetArea(), oldDelay = GetDelay();
        auto pOldNtk = Abc_NtkDup(GetNet());
        if (fVerb)
            cout << "oldArea = " << oldArea << ", " << "oldDelay = " << oldDelay << endl;
        Comm("st; compress2rs; dch; amap;", fVerb);
        TopoSort();
        double newArea = GetArea(), newDelay = GetDelay();
        if (fVerb)
            cout << "newArea = " << newArea << ", " << "newDelay = " << newDelay << endl;
        if (newDelay <= maxDelay) {
            auto impr = UpdNet(oldArea, oldDelay, pOldNtk, newArea, newDelay, ORIENT::AREA);
            if (impr != IMPR::GOOD) {
                cont = false;
                if (fVerb)
                    cout << "reject" << endl;
            }
            else {
                if (fVerb)
                    cout << "accept" << endl;
            }
        }
        else {
            SetMainNet(pOldNtk);
                cont = false;
            if (fVerb)
                cout << "reject" << endl;
        }
    }
    PrintStat();
}


pair <double, double> AbcMan::Map(MAP_TYPE cell, ORIENT orient, bool fVerb) {
    double oldArea = numeric_limits <double>::max();
    double oldDelay = numeric_limits <double>::max();
    const int LUT_INP = 6;
    ostringstream LutInpStr("");
    LutInpStr << LUT_INP;
    if ((cell == MAP_TYPE::SCL && GetNetType() == NET_TYPE::GATE) ||
        (cell == MAP_TYPE::LUT && IsLutNet())) {
        oldArea = GetArea();
        oldDelay = GetDelay();
    }
    bool isFirst = true;
    bool isCont = true;
    while (isCont) {
        auto oldNtk = Abc_NtkDup(GetNet());
        if (isFirst) {
            Comm("st; dch;", fVerb);
            isFirst = false;
        }
        else
            Comm("st; b;", fVerb);
        if (cell == MAP_TYPE::SCL) {
            if (orient == ORIENT::AREA)
                Comm("amap", fVerb);
            else if (orient == ORIENT::DELAY)
                Comm("map", fVerb);
            else
                assert(0);
        }
        else if (cell == MAP_TYPE::LUT) {
            if (orient == ORIENT::AREA)
                Comm("if -a -K " + LutInpStr.str(), fVerb);
            else if (orient == ORIENT::DELAY)
                Comm("if -K " + LutInpStr.str(), fVerb);
            else
                assert(0);
        }
        else
            assert(0);
        double newArea = GetArea();
        double newDelay = GetDelay();
        IMPR impr = UpdNet(oldArea, oldDelay, oldNtk, newArea, newDelay, orient);
        if (impr == IMPR::GOOD) {
            oldArea = newArea;
            oldDelay = newDelay;
        }
        else
            isCont = false;
        // PrintStat();
    }
    return make_pair(oldArea, oldDelay);
}


// pair <double, double> AbcMan::Map2(double maxDelay, bool fVerb) {
//     double oldArea = numeric_limits <double>::max();
//     double oldDelay = numeric_limits <double>::max();
//     ostringstream LutInpStr("");
//     LutInpStr << LUT_INP;
//     assert(GetNetType() == NET_TYPE::STRASH);
//     bool isFirst = true;
//     bool isCont = true;
//     while (isCont) {
//         auto oldNtk = Abc_NtkDup(GetNet());
//         if (isFirst) {
//             Comm("st; dch;", fVerb);
//             isFirst = false;
//         }
//         else
//             Comm("st; b;", fVerb);
//         ostringstream oss("");
//         oss << "map -D " << maxDelay;
//         Comm(oss.str(), fVerb);
//         double newArea = GetArea();
//         double newDelay = GetDelay();
//         IMPR impr = UpdNet(oldArea, oldDelay, oldNtk, newArea, newDelay, ORIENT::AREA);
//         if (impr == IMPR::GOOD) {
//             oldArea = newArea;
//             oldDelay = newDelay;
//         }
//         else
//             isCont = false;
//         // PrintStat();
//     }
//     return make_pair(oldArea, oldDelay);
// }


IMPR AbcMan::UpdNet(double oldArea, double oldDelay, Abc_Ntk_t * oldNtk, double newArea, double newDelay, ORIENT orient) {
    IMPR impr = IMPR::SAME;
    if (orient == ORIENT::AREA) {
        if (DoubleGreat(newArea, oldArea) || (DoubleEqual(newArea, oldArea) && DoubleGreat(newDelay, oldDelay)))
            impr = IMPR::BAD;
        else if (DoubleEqual(newArea, oldArea) && DoubleEqual(newDelay, oldDelay))
            impr = IMPR::SAME;
        else
            impr = IMPR::GOOD;
    }
    else if (orient == ORIENT::DELAY) {
        if (DoubleGreat(newDelay, oldDelay) || (DoubleEqual(newDelay, oldDelay) && DoubleGreat(newArea, oldArea)))
            impr = IMPR::BAD;
        else if (DoubleEqual(newDelay, oldDelay) && DoubleEqual(newArea, oldArea))
            impr = IMPR::SAME;
        else
            impr = IMPR::GOOD;
    }
    else
        assert(0);
    if (impr == IMPR::BAD) {
        assert(oldArea != numeric_limits <double>::max() && oldDelay != numeric_limits <double>::max());
        assert(oldNtk != nullptr);
        // cout << "Cancel the last abc command" << endl;
        SetMainNet(oldNtk);
    }
    else
        Abc_NtkDelete(oldNtk);
    return impr;
}


/**
 * @brief Get the network type
 * 
 * @param pNtk            the network
 * @return NET_TYPE       the network type
 */
NET_TYPE AbcMan::GetNetType(Abc_Ntk_t * pNtk) const {
    if (Abc_NtkIsAigLogic(pNtk))
        return NET_TYPE::AIG;
    else if (Abc_NtkIsMappedLogic(pNtk))
        return NET_TYPE::GATE;
    else if (Abc_NtkIsSopLogic(pNtk))
        return NET_TYPE::SOP;
    else if (Abc_NtkIsStrash(pNtk))
        return NET_TYPE::STRASH;
    else {
        fmt::print(stderr, "invalid network type\n");
        assert(0);
    }
}


double AbcMan::GetArea(Abc_Ntk_t * pNtk) const {
    auto type = GetNetType(pNtk);
    if (type == NET_TYPE::AIG || type == NET_TYPE::STRASH)
        return Abc_NtkNodeNum(pNtk);
    else if (type == NET_TYPE::SOP) {
        AbcObj* pObj = nullptr;
        int i = 0;
        int ret = Abc_NtkNodeNum(pNtk);
        Abc_NtkForEachNode(pNtk, pObj, i) {
            if (Abc_NodeIsConst(pObj))
                --ret;
        }
        return ret;
    }
    else if (type == NET_TYPE::GATE) {
        auto pLibScl = static_cast <SC_Lib *> (GetAbcFrame()->pLibScl);
        if (pLibScl == nullptr)
            return Abc_NtkGetMappedArea(pNtk);
        else {
            return Abc_NtkGetMappedArea(pNtk);
            // assert(pNtk->nBarBufs2 == 0);
            // assert(CheckSCLNet(pNtk));
            // SC_Man * p = Abc_SclManStart(pLibScl, pNtk, 0, 1, 0, 0);
            // double area = Abc_SclGetTotalArea(p->pNtk);
            // Abc_SclManFree(p);
            // return area;
        }
    }
    else
        assert(0);
}


double AbcMan::GetDelay(Abc_Ntk_t * pNtk) const {
    auto type = GetNetType(pNtk);
    if (type == NET_TYPE::AIG || type == NET_TYPE::SOP || type == NET_TYPE::STRASH)
        return Abc_NtkLevel(pNtk);
    else if (type == NET_TYPE::GATE) {
        auto pLibScl = static_cast <SC_Lib *> (GetAbcFrame()->pLibScl);
        if (pLibScl == nullptr) {
            return Abc_NtkDelayTrace(pNtk, nullptr, nullptr, 0);
        }
        else {
            assert(pNtk->nBarBufs2 == 0);
            assert(CheckSCLNet(pNtk));
            SC_Man * p = Abc_SclManStart(pLibScl, pNtk, 0, 1, 0, 0);
            int fRise = 0;
            AbcObj* pPivot = Abc_SclFindCriticalCo(p, &fRise); 
            double delay = Abc_SclObjTimeOne(p, pPivot, fRise);
            AbcObj* pObj = nullptr;
            int i = 0;
            Abc_NtkForEachObj(pNtk, pObj, i)
                pObj->dTemp = Abc_SclObjTimeMax(p, pObj);
            Abc_SclManFree(p);
            return delay;
        }
    }
    else
        assert(0);
}


bool AbcMan::CheckSCLNet(abc::Abc_Ntk_t * pNtk) const {
    AbcObj* pObj, * pFanin;
    int i, k, fFlag = 1;
    Abc_NtkIncrementTravId( pNtk );        
    Abc_NtkForEachCi( pNtk, pObj, i )
        Abc_NodeSetTravIdCurrent( pObj );
    Abc_NtkForEachNode( pNtk, pObj, i )
    {
        Abc_ObjForEachFanin( pObj, pFanin, k )
            if ( !Abc_NodeIsTravIdCurrent( pFanin ) )
                printf( "obj %d and its fanin %d are not in the topo order\n", Abc_ObjId(pObj), Abc_ObjId(pFanin) ), fFlag = 0;
        Abc_NodeSetTravIdCurrent( pObj );
        if ( Abc_ObjIsBarBuf(pObj) )
            continue;
        // if ( Abc_ObjFanoutNum(pObj) == 0 )
        //     printf( "node %d has no fanout\n", Abc_ObjId(pObj) ), fFlag = 0;
        if ( !fFlag )
            break;
    }
    // if ( fFlag && fVerbose )
    //     printf( "The network is in topo order and no dangling nodes.\n" );
    return fFlag;
}


AbcObj* AbcMan::GetTwinNode( AbcObj* pNode ) {
    assert( Abc_NtkHasMapping(pNode->pNtk) );
    Mio_Gate_t * pGate = (Mio_Gate_t *)pNode->pData;
    if ( pGate == nullptr || Mio_GateReadTwin(pGate) == nullptr )
        return nullptr;
    AbcObj* pNode2 = nullptr;
    int id = 0;
    AbcObj* pTwin = nullptr;
    int count = 0;
    Abc_NtkForEachNode(pNode->pNtk, pNode2, id) {
        if ( Abc_ObjFaninNum(pNode) != Abc_ObjFaninNum(pNode2) )
            continue;
        bool sameFanin = true;
        for (int faninId = 0; faninId < Abc_ObjFaninNum(pNode); ++faninId) {
            if (Abc_ObjFanin(pNode, faninId) != Abc_ObjFanin(pNode2, faninId)) {
                sameFanin = false;
                break;
            }
        }
        if (!sameFanin)
            continue;
        if ( Mio_GateReadTwin(pGate) != (Mio_Gate_t *)pNode2->pData )
            continue;
        pTwin = pNode2;
        ++count;
        if (count > 1)
            assert(0);
    }
    return pTwin;
}


void AbcMan::LoadAlias() {
    Comm("alias hi history", false);
    Comm("alias b balance", false);
    Comm("alias cg clockgate", false);
    Comm("alias cl cleanup", false);
    Comm("alias clp collapse", false);
    Comm("alias cs care_set", false);
    Comm("alias el eliminate", false);
    Comm("alias esd ext_seq_dcs", false);
    Comm("alias f fraig", false);
    Comm("alias fs fraig_sweep", false);
    Comm("alias fsto fraig_store", false);
    Comm("alias fres fraig_restore", false);
    Comm("alias fr fretime", false);
    Comm("alias ft fraig_trust", false);
    Comm("alias ic indcut", false);
    Comm("alias lp lutpack", false);
    Comm("alias pcon print_cone", false);
    Comm("alias pd print_dsd", false);
    Comm("alias pex print_exdc -d", false);
    Comm("alias pf print_factor", false);
    Comm("alias pfan print_fanio", false);
    Comm("alias pg print_gates", false);
    Comm("alias pl print_level", false);
    Comm("alias plat print_latch", false);
    Comm("alias pio print_io", false);
    Comm("alias pk print_kmap", false);
    Comm("alias pm print_miter", false);
    Comm("alias ps print_stats ", false);
    Comm("alias psb print_stats -b", false);
    Comm("alias psu print_supp", false);
    Comm("alias psy print_symm", false);
    Comm("alias pun print_unate", false);
    Comm("alias q quit", false);
    Comm("alias r read", false);
    Comm("alias ra read_aiger", false);
    Comm("alias r3 retime -M 3", false);
    Comm("alias r3f retime -M 3 -f", false);
    Comm("alias r3b retime -M 3 -b", false);
    Comm("alias ren renode", false);
    Comm("alias rh read_hie", false);
    Comm("alias ri read_init", false);
    Comm("alias rl read_blif", false);
    Comm("alias rb read_bench", false);
    Comm("alias ret retime", false);
    Comm("alias dret dretime", false);
    Comm("alias rp read_pla", false);
    Comm("alias rt read_truth", false);
    Comm("alias rv read_verilog", false);
    Comm("alias rvl read_verlib", false);
    Comm("alias rsup read_super mcnc5_old.super", false);
    Comm("alias rlib read_library", false);
    Comm("alias rlibc read_library cadence.genlib", false);
    Comm("alias rty read_liberty", false);
    Comm("alias rlut read_lut", false);
    Comm("alias rw rewrite", false);
    Comm("alias rwz rewrite -z", false);
    Comm("alias rf refactor", false);
    Comm("alias rfz refactor -z", false);
    Comm("alias re restructure", false);
    Comm("alias rez restructure -z", false);
    Comm("alias rs resub", false);
    Comm("alias rsz resub -z", false);
    Comm("alias sa set autoexec ps", false);
    Comm("alias scl scleanup", false);
    Comm("alias sif if -s", false);
    Comm("alias so source -x", false);
    Comm("alias st strash", false);
    Comm("alias sw sweep", false);
    Comm("alias ssw ssweep", false);
    Comm("alias tr0 trace_start", false);
    Comm("alias tr1 trace_check", false);
    Comm("alias trt \"r c.blif; st; tr0; b; tr1\"", false);
    Comm("alias u undo", false);
    Comm("alias w write", false);
    Comm("alias wa write_aiger", false);
    Comm("alias wb write_bench", false);
    Comm("alias wc write_cnf", false);
    Comm("alias wh write_hie", false);
    Comm("alias wl write_blif", false);
    Comm("alias wp write_pla", false);
    Comm("alias wv write_verilog", false);
    Comm("alias resyn       \"b; rw; rwz; b; rwz; b\"", false);
    Comm("alias resyn2      \"b; rw; rf; b; rw; rwz; b; rfz; rwz; b\"", false);
    Comm("alias resyn2a     \"b; rw; b; rw; rwz; b; rwz; b\"", false);
    Comm("alias resyn3      \"b; rs; rs -K 6; b; rsz; rsz -K 6; b; rsz -K 5; b\"", false);
    Comm("alias compress    \"b -l; rw -l; rwz -l; b -l; rwz -l; b -l\"", false);
    Comm("alias compress2   \"b -l; rw -l; rf -l; b -l; rw -l; rwz -l; b -l; rfz -l; rwz -l; b -l\"", false);
    Comm("alias choice      \"fraig_store; resyn; fraig_store; resyn2; fraig_store; fraig_restore\"", false);
    Comm("alias choice2     \"fraig_store; balance; fraig_store; resyn; fraig_store; resyn2; fraig_store; resyn2; fraig_store; fraig_restore\"", false);
    Comm("alias rwsat       \"st; rw -l; b -l; rw -l; rf -l\"", false);
    Comm("alias drwsat2     \"st; drw; b -l; drw; drf; ifraig -C 20; drw; b -l; drw; drf\"", false);
    Comm("alias share       \"st; multi -m; sop; fx; resyn2\"", false);
    Comm("alias addinit     \"read_init; undc; strash; zero\"", false);
    Comm("alias blif2aig    \"undc; strash; zero\"", false);
    Comm("alias v2p         \"&vta_gla; &ps; &gla_derive; &put; w 1.aig; pdr -v\"", false);
    Comm("alias g2p         \"&ps; &gla_derive; &put; w 2.aig; pdr -v\"", false);
    Comm("alias &sw_        \"&put; sweep; st; &get\"", false);
    Comm("alias &fx_        \"&put; sweep; sop; fx; st; &get\"", false);
    Comm("alias &dc3        \"&b; &jf -K 6; &b; &jf -K 4; &b\"", false);
    Comm("alias &dc4        \"&b; &jf -K 7; &fx; &b; &jf -K 5; &fx; &b\"", false);
    Comm("alias src_rw      \"st; rw -l; rwz -l; rwz -l\"", false);
    Comm("alias src_rs      \"st; rs -K 6 -N 2 -l; rs -K 9 -N 2 -l; rs -K 12 -N 2 -l\"", false);
    Comm("alias src_rws     \"st; rw -l; rs -K 6 -N 2 -l; rwz -l; rs -K 9 -N 2 -l; rwz -l; rs -K 12 -N 2 -l\"", false);
    Comm("alias resyn2rs    \"b; rs -K 6; rw; rs -K 6 -N 2; rf; rs -K 8; b; rs -K 8 -N 2; rw; rs -K 10; rwz; rs -K 10 -N 2; b; rs -K 12; rfz; rs -K 12 -N 2; rwz; b\"", false);
    Comm("alias compress2rs \"b -l; rs -K 6 -l; rw -l; rs -K 6 -N 2 -l; rf -l; rs -K 8 -l; b -l; rs -K 8 -N 2 -l; rw -l; rs -K 10 -l; rwz -l; rs -K 10 -N 2 -l; b -l; rs -K 12 -l; rfz -l; rs -K 12 -N 2 -l; rwz -l; b -l\"", false);
    Comm("alias fix_aig     \"logic; undc; strash; zero\"", false);
    Comm("alias fix_blif    \"undc; strash; zero\"", false);
    Comm("alias recadd3     \"st; rec_add3; b; rec_add3; dc2; rec_add3; if -K 8; bidec; st; rec_add3; dc2; rec_add3; if -g -K 6; st; rec_add3\"", false);
}


static AbcObj* Abc_NtkDupObj_KeepName( Abc_Ntk_t * pNtkNew, AbcObj* pObj, int fCopyName ) {
    AbcObj* pObjNew;
    // create the new object
    pObjNew = Abc_NtkCreateObj( pNtkNew, (Abc_ObjType_t)pObj->Type );
    // transfer names of the terminal objects
    if ( fCopyName )
    {
        if ( Abc_ObjIsCi(pObj) )
        {
            if ( !Abc_NtkIsNetlist(pNtkNew) )
                Abc_ObjAssignName( pObjNew, Abc_ObjName(Abc_ObjFanout0Ntk(pObj)), NULL );
        }
        else if ( Abc_ObjIsCo(pObj) )
        {
            if ( !Abc_NtkIsNetlist(pNtkNew) )
            {
                if ( Abc_ObjIsPo(pObj) )
                    Abc_ObjAssignName( pObjNew, Abc_ObjName(Abc_ObjFanin0Ntk(pObj)), NULL );
                else
                {
                    assert( Abc_ObjIsLatch(Abc_ObjFanout0(pObj)) );
                    Abc_ObjAssignName( pObjNew, Abc_ObjName(pObj), NULL );
                }
            }
        }
        else if ( Abc_ObjIsBox(pObj) || Abc_ObjIsNet(pObj) || Abc_ObjIsNode(pObj) )
            Abc_ObjAssignName( pObjNew, Abc_ObjName(pObj), NULL );
    }
    // copy functionality/names
    if ( Abc_ObjIsNode(pObj) ) // copy the function if functionality is compatible
    {
        if ( pNtkNew->ntkFunc == pObj->pNtk->ntkFunc ) 
        {
            if ( Abc_NtkIsStrash(pNtkNew) ) 
            {}
            else if ( Abc_NtkHasSop(pNtkNew) || Abc_NtkHasBlifMv(pNtkNew) )
                pObjNew->pData = Abc_SopRegister( (Mem_Flex_t *)pNtkNew->pManFunc, (char *)pObj->pData );
#ifdef ABC_USE_CUDD
            else if ( Abc_NtkHasBdd(pNtkNew) )
                pObjNew->pData = Cudd_bddTransfer((DdManager *)pObj->pNtk->pManFunc, (DdManager *)pNtkNew->pManFunc, (DdNode *)pObj->pData), Cudd_Ref((DdNode *)pObjNew->pData);
#endif
            else if ( Abc_NtkHasAig(pNtkNew) )
                pObjNew->pData = Hop_Transfer((Hop_Man_t *)pObj->pNtk->pManFunc, (Hop_Man_t *)pNtkNew->pManFunc, (Hop_Obj_t *)pObj->pData, Abc_ObjFaninNum(pObj));
            else if ( Abc_NtkHasMapping(pNtkNew) )
                pObjNew->pData = pObj->pData, pNtkNew->nBarBufs2 += !pObj->pData;
            else assert( 0 );
        }
    }
    else if ( Abc_ObjIsNet(pObj) ) // copy the name
    {
    }
    else if ( Abc_ObjIsLatch(pObj) ) // copy the reset value
        pObjNew->pData = pObj->pData;
    pObjNew->fPersist = pObj->fPersist;
    // transfer HAIG
//    pObjNew->pEquiv = pObj->pEquiv;
    // remember the new node in the old node
    pObj->pCopy = pObjNew;
    return pObjNew;
}


static Abc_Ntk_t * Abc_NtkDup_KeepName( Abc_Ntk_t * pNtk ) {
    Abc_Ntk_t * pNtkNew; 
    AbcObj* pObj, * pFanin;
    int i, k;
    if ( pNtk == NULL )
        return NULL;
    // start the network
    pNtkNew = Abc_NtkStartFrom( pNtk, pNtk->ntkType, pNtk->ntkFunc );
    // copy the internal nodes
    if ( Abc_NtkIsStrash(pNtk) )
    {
        // copy the AND gates
        Abc_AigForEachAnd( pNtk, pObj, i )
            pObj->pCopy = Abc_AigAnd( (Abc_Aig_t *)pNtkNew->pManFunc, Abc_ObjChild0Copy(pObj), Abc_ObjChild1Copy(pObj) );
        // relink the choice nodes
        Abc_AigForEachAnd( pNtk, pObj, i )
            if ( pObj->pData )
                pObj->pCopy->pData = ((AbcObj*)pObj->pData)->pCopy;
        // relink the CO nodes
        Abc_NtkForEachCo( pNtk, pObj, i )
            Abc_ObjAddFanin( pObj->pCopy, Abc_ObjChild0Copy(pObj) );
        // get the number of nodes before and after
        if ( Abc_NtkNodeNum(pNtk) != Abc_NtkNodeNum(pNtkNew) )
            printf( "Warning: Structural hashing during duplication reduced %d nodes (this is a minor bug).\n",
                Abc_NtkNodeNum(pNtk) - Abc_NtkNodeNum(pNtkNew) );
    }
    else
    {
        // duplicate the nets and nodes (CIs/COs/latches already dupped)
        Abc_NtkForEachObj( pNtk, pObj, i )
            if ( pObj->pCopy == NULL )
                Abc_NtkDupObj_KeepName(pNtkNew, pObj, 1);
        // reconnect all objects (no need to transfer attributes on edges)
        Abc_NtkForEachObj( pNtk, pObj, i )
            if ( !Abc_ObjIsBox(pObj) && !Abc_ObjIsBo(pObj) )
                Abc_ObjForEachFanin( pObj, pFanin, k )
                    Abc_ObjAddFanin( pObj->pCopy, pFanin->pCopy );
    }
    // duplicate the EXDC Ntk
    if ( pNtk->pExdc )
        pNtkNew->pExdc = Abc_NtkDup( pNtk->pExdc );
    if ( pNtk->pExcare )
        pNtkNew->pExcare = Abc_NtkDup( (Abc_Ntk_t *)pNtk->pExcare );
    // duplicate timing manager
    if ( pNtk->pManTime )
        Abc_NtkTimeInitialize( pNtkNew, pNtk );
    if ( pNtk->vPhases )
        Abc_NtkTransferPhases( pNtkNew, pNtk );
    if ( pNtk->pWLoadUsed )
        pNtkNew->pWLoadUsed = Abc_UtilStrsav( pNtk->pWLoadUsed );
    // check correctness
    if ( !Abc_NtkCheck( pNtkNew ) )
        fprintf( stdout, "Abc_NtkDup(): Network check has failed.\n" );
    pNtk->pCopy = pNtkNew;
    return pNtkNew;
}


NetMan::NetMan(): AbcMan(), pNtk(nullptr), isDupl(true) {
}


NetMan::NetMan(Abc_Ntk_t * p_ntk, bool is_dupl): AbcMan(), isDupl(is_dupl) {
    if (is_dupl)
        pNtk = Abc_NtkDup_KeepName(p_ntk);
    else
        pNtk = p_ntk;
}


NetMan::NetMan(const std::string& fileName): AbcMan(), isDupl(true) {
    AbcMan::ReadNet(fileName);
    pNtk = Abc_NtkDup_KeepName(AbcMan::GetNet());
}


NetMan::~NetMan() {
    if (isDupl && pNtk != AbcMan::GetNet()) {
        if (pNtk != nullptr) {
            Abc_NtkDelete(pNtk);
            pNtk = nullptr;
        }
    }
}


NetMan::NetMan(const NetMan& net_man): AbcMan(), isDupl(true) {
    pNtk = Abc_NtkDup_KeepName(net_man.pNtk);
}


NetMan::NetMan(NetMan&& net_man): AbcMan(), pNtk(net_man.pNtk), isDupl(net_man.isDupl) {
    net_man.isDupl = false;
    net_man.pNtk = nullptr;
}


NetMan& NetMan::operator = (const NetMan& net_man) {
    if (this == &net_man)
        return *this;
    if (isDupl && pNtk != nullptr && pNtk != AbcMan::GetNet() && pNtk != net_man.GetNet())
        Abc_NtkDelete(pNtk);
    pNtk = Abc_NtkDup_KeepName(net_man.GetNet());
    isDupl = true;
    return *this;
}


NetMan& NetMan::operator = (NetMan&& net_man) {
    // cout << "move assign netman" << endl;
    if (this == &net_man)
        return *this;
    if (isDupl && pNtk != nullptr && pNtk != AbcMan::GetNet() && pNtk != net_man.GetNet())
        Abc_NtkDelete(pNtk);
    pNtk = net_man.pNtk;
    isDupl = net_man.isDupl;
    net_man.isDupl = false;
    net_man.pNtk = nullptr;
    return *this;
}


/**
 * @brief Get constant node IDs in the network.
 * 
 * @param fVerb           Whether to print the information.
 * @return The constant node IDs.
 */
IntPair NetMan::GetConstIds(bool fVerb) {
    IntPair ret(-1, -1);
    auto type = GetNetType();
    AbcObj* pObj = nullptr;
    int i = 0;
    Abc_NtkForEachNode(GetNet(), pObj, i) {
        if (type == NET_TYPE::GATE || type == NET_TYPE::SOP) {
            if (Abc_NodeIsConst0(pObj)) {
                if (fVerb)
                    fmt::print("find const 0: {}\n", *pObj);
                if (ret.first == -1)
                    ret.first = GetId(pObj);
            }
            else if (Abc_NodeIsConst1(pObj)) {
                if (fVerb)
                    fmt::print("find const 1: {}\n", *pObj);
                if (ret.second == -1)
                    ret.second = GetId(pObj);
            }
        }
        else if (type == NET_TYPE::AIG) {
            auto pHopObj = static_cast <Hop_Obj_t *> (pObj->pData);
            auto pHopObjR = Hop_Regular(pHopObj);
            if (Hop_ObjIsConst1(pHopObjR)) {
                assert(Hop_ObjFanin0(pHopObjR) == nullptr);
                assert(Hop_ObjFanin1(pHopObjR) == nullptr);
                if (!Hop_IsComplement(pHopObj))
                    ret.second = GetId(pObj);
                else 
                    ret.first = GetId(pObj);
            }
        }
        else if (type == NET_TYPE::STRASH) {
            ret.first = -1;
            ret.second = Abc_AigConst1(NetMan::GetNet())->Id;
        }
        else
            assert(0);
    }
    return ret;
}


/**
 * @brief Create constant nodes in the network.
 * @brief If the constant nodes exist, return their IDs.
 * @brief Otherwise, create the constant nodes and return their IDs.
 * 
 * @param fVerb         Whether to print the information.
 * @return The constant node IDs.
 */
IntPair NetMan::CreateConstsIfNotExist(bool fVerb) {
    auto consts = GetConstIds(fVerb);
    IntPair ret(consts);
    if (ret.first == -1) {
        auto pObj = Abc_NtkCreateNodeConst0(GetNet());
        Rename(pObj, "zero");
        ret.first = GetId(pObj);
    }
    if (ret.second == -1) {
        auto pObj = Abc_NtkCreateNodeConst1(GetNet());
        Rename(pObj, "one");
        ret.second = GetId(pObj);
    }
    return ret;
}


/**
 * @brief Merge the constant nodes in the network.
 * 
 * @retval pNtk    The network after merging the constant nodes.
 * @return void
 */
void NetMan::MergeConst(bool fVerb) {
    IntPair ret(-1, -1);
    auto type = GetNetType();
    AbcObj* pObj = nullptr;
    int i = 0;
    Abc_NtkForEachNode(GetNet(), pObj, i) {
        if (type == NET_TYPE::GATE || type == NET_TYPE::SOP) {
            if (Abc_NodeIsConst0(pObj)) {
                if (ret.first == -1) {
                    if (fVerb)
                        fmt::print("find const 0: {}\n", *pObj);
                    ret.first = GetId(pObj);
                }
                else {
                    if (fVerb)
                        fmt::print("merge const 0: {} -> {}\n", *pObj, *GetObj(ret.first));
                    Abc_ObjReplace(pObj, GetObj(ret.first));
                }
            }
            else if (Abc_NodeIsConst1(pObj)) {
                if (ret.second == -1) {
                    if (fVerb)
                        fmt::print("find const 1: {}\n", *pObj);
                    ret.second = GetId(pObj);
                }
                else {
                    if (fVerb)
                        fmt::print("merge const 1: {} -> {}\n", *pObj, *GetObj(ret.second));
                    Abc_ObjReplace(pObj, GetObj(ret.second));
                }
            }
        }
        else
            assert(0);
    }
}


/**
 * @brief Calculate the topological order of the nodes in the network.
 * @param inclConst Whether to include the constant nodes.
 * 
 * @return AbcObjVect   The topological order of the nodes in the network.
 */
AbcObjVect NetMan::CalcTopoOrd(bool inclConst) const {
    AbcObjVect nodes;
    nodes.reserve(GetNodeNum());
    SetNetNotTrav();
    for (int i = 0; i < GetPoNum(); ++i) {
        auto pDriver = GetFanin(GetPo(i), 0);
        if (!GetObjTrav(pDriver))
            CalcTopoOrdRec(pDriver, nodes, inclConst);
    }
    return nodes;
}


/**
 * @brief Recursively calculate the topological order of the nodes in the network.
 * 
 * @param pObj      The current node.
 * @param nodes     A vector to store the topological order of the nodes.
 * @param inclConst Whether to include the constant nodes.
 * @retval nodes    A vector to store the topological order of the nodes.
 * @return void
 */
void NetMan::CalcTopoOrdRec(AbcObj* pObj, AbcObjVect& nodes, bool inclConst) const {
    if (!IsNode(pObj) || IsConst(pObj->Id))
        return;
    SetObjTrav(pObj);
    for (int i = 0; i < GetFaninNum(pObj); ++i) {
        auto pFanin = GetFanin(pObj, i);
        if (!GetObjTrav(pFanin))
            CalcTopoOrdRec(pFanin, nodes, inclConst);
    }
    nodes.emplace_back(pObj);
}


/**
 * @brief Calculate the topological order of the node ids in the network.

 * @param inclConst Whether to include the constant nodes.
 * @return IntVect   The topological order of the node ids in the network.
 */
IntVect NetMan::CalcTopoOrdOfIds(bool inclConst) const {
    IntVect nodes;
    nodes.reserve(GetNodeNum());
    SetNetNotTrav();
    for (int i = 0; i < GetPoNum(); ++i) {
        auto pDriver = GetFanin(GetPo(i), 0);
        if (!GetObjTrav(pDriver))
            CalcTopoOrdOfIdsRec(pDriver, nodes, inclConst);
    }
    return nodes;
}


/**
 * @brief Recursively calculate the topological order of the node ids in the network.
 * 
 * @param pObj      The current node.
 * @param nodes     A vector to store the topological order of the node ids.
 * @param inclConst Whether to include the constant nodes.
 * @retval nodes    A vector to store the topological order of the node ids.
 * @return void
 */
void NetMan::CalcTopoOrdOfIdsRec(AbcObj* pObj, IntVect& nodes, bool inclConst) const {
    if (!IsNode(pObj) || IsConst(pObj->Id))
        return;
    SetObjTrav(pObj);
    for (int i = 0; i < GetFaninNum(pObj); ++i) {
        auto pFanin = GetFanin(pObj, i);
        if (!GetObjTrav(pFanin))
            CalcTopoOrdOfIdsRec(pFanin, nodes, inclConst);
    }
    nodes.emplace_back(GetId(pObj));
}


AbcObjVect NetMan::GetTFI(AbcObj* pObj) const {
    AbcObjVect nodes;
    nodes.reserve(GetNodeNum());
    SetNetNotTrav();
    for (int i = 0; i < GetFaninNum(pObj); ++i) {
        auto pFanin = GetFanin(pObj, i);
        if (!GetObjTrav(pFanin))
            GetTFIRec(pFanin, nodes);
    }
    return nodes;
}


void NetMan::GetTFIRec(AbcObj* pObj, AbcObjVect& nodes) const {
    if (!IsNode(pObj))
        return;
    SetObjTrav(pObj);
    for (int i = 0; i < GetFaninNum(pObj); ++i) {
        auto pFanin = GetFanin(pObj, i);
        if (!GetObjTrav(pFanin))
            GetTFIRec(pFanin, nodes);
    }
    nodes.emplace_back(pObj);
}


IntVect NetMan::GetTFI(AbcObj* pObj, const set <int>& critGraph) const {
    IntVect objs;
    objs.reserve(GetNodeNum());
    SetNetNotTrav();
    for (int i = 0; i < GetFaninNum(pObj); ++i) {
        auto pFanin = GetFanin(pObj, i);
        if (!GetObjTrav(pFanin))
            GetTFIRec(pFanin, objs, critGraph);
    }
    return objs;
}


void NetMan::GetTFIRec(AbcObj* pObj, IntVect& objs, const set <int>& critGraph) const {
    if (critGraph.count(pObj->Id) == 0)
        return;
    // if (!IsNode(pObj))
    //     return;
    SetObjTrav(pObj);
    for (int i = 0; i < GetFaninNum(pObj); ++i) {
        auto pFanin = GetFanin(pObj, i);
        if (!GetObjTrav(pFanin))
            GetTFIRec(pFanin, objs, critGraph);
    }
    objs.emplace_back(pObj->Id);
}


AbcObjVect NetMan::GetTFO(AbcObj* pObj) const {
    AbcObjVect nodes;
    nodes.reserve(GetNodeNum());
    SetNetNotTrav();
    for (int i = 0; i < GetFanoutNum(pObj); ++i) {
        auto pFanout = GetFanout(pObj, i);
        if (!GetObjTrav(pFanout))
            GetTFORec(pFanout, nodes);
    }
    reverse(nodes.begin(), nodes.end());
    return nodes;
}


void NetMan::GetTFORec(AbcObj* pObj, AbcObjVect& nodes) const {
    if (!IsNode(pObj))
        return;
    SetObjTrav(pObj);
    for (int i = 0; i < GetFanoutNum(pObj); ++i) {
        auto pFanout = GetFanout(pObj, i);
        if (!GetObjTrav(pFanout))
            GetTFORec(pFanout, nodes);
    }
    nodes.emplace_back(pObj);
}


IntVect NetMan::GetTFO(AbcObj* pObj, const set <int>& critGraph) const {
    IntVect objs;
    objs.reserve(GetNodeNum());
    SetNetNotTrav();
    for (int i = 0; i < GetFanoutNum(pObj); ++i) {
        auto pFanout = GetFanout(pObj, i);
        if (!GetObjTrav(pFanout))
            GetTFORec(pFanout, objs, critGraph);
    }
    reverse(objs.begin(), objs.end());
    return objs;
}


void NetMan::GetTFORec(AbcObj* pObj, IntVect& objs, const set <int>& critGraph) const {
    if (critGraph.count(pObj->Id) == 0)
        return;
    // if (!IsNode(pObj))
    //     return;
    SetObjTrav(pObj);
    for (int i = 0; i < GetFanoutNum(pObj); ++i) {
        auto pFanout = GetFanout(pObj, i);
        if (!GetObjTrav(pFanout))
            GetTFORec(pFanout, objs, critGraph);
    }
    objs.emplace_back(pObj->Id);
}


void NetMan::Comm(const std::string& cmd, bool fVerb) {
    assert(isDupl == true);
    AbcMan::SetMainNet(pNtk); // abc manage the memory of the old network
    AbcMan::Comm(cmd, fVerb);
    pNtk = Abc_NtkDup(AbcMan::GetNet()); // NetMan manage the memory of the duplicated network
}


void NetMan::Synth(ORIENT orient, bool fVerb) {
    assert(isDupl == true);
    AbcMan::SetMainNet(pNtk); // abc manage the memory of the old network
    AbcMan::Synth(orient, fVerb);
    pNtk = Abc_NtkDup(AbcMan::GetNet()); // NetMan manage the memory of the duplicated network
}


void NetMan::SynthAndMap(double maxDelay, bool fVerb) {
    assert(isDupl == true);
    AbcMan::SetMainNet(pNtk); // abc manage the memory of the old network
    AbcMan::SynthAndMap(maxDelay, fVerb);
    pNtk = Abc_NtkDup(AbcMan::GetNet()); // NetMan manage the memory of the duplicated network
}


/**
 * @brief Print the network.
 * 
 * @param showFunct     Whether to show the functionality of the nodes.
 * @return void
 */
void NetMan::Print(bool showFunct) const {
    fmt::print("{}\n", GetNetName());
    for (int i = 0; i < GetIdMaxPlus1(); ++i) {
        if (IsObj(i))
            PrintObj(i, showFunct);
    }
}


void NetMan::PrintObjBas(AbcObj* pObj, string&& endWith) const {
    cout << GetName(pObj) << "(" << GetId(pObj) << ")" << endWith;
}


void NetMan::PrintObj(AbcObj* pObj, bool showFunct) const {
    PrintObjBas(pObj, ":");
    for (int i = 0; i < GetFaninNum(pObj); ++i)
        PrintObjBas(GetFanin(pObj, i), ",");
    if (showFunct) {
        if (NetMan::GetNetType() == NET_TYPE::SOP) {
            if (IsNode(pObj)) {
                auto pSop = static_cast <char *> (pObj->pData);
                for (auto pCh = pSop; *pCh != '\0'; ++pCh) {
                    if (*pCh != '\n')
                        cout << *pCh;
                    else
                        cout << "\\n";
                }
                cout << endl;
            }
            else
                cout << endl;
        }
        else if (NetMan::GetNetType() == NET_TYPE::GATE) {
            if (IsNode(pObj))
                cout << Mio_GateReadName(static_cast <Mio_Gate_t *> (pObj->pData)) << endl;
            else
                cout << endl;
        }
        else if (NetMan::GetNetType() == NET_TYPE::STRASH) {
            if (Abc_AigNodeIsAnd(pObj)) {
                assert(!Abc_ObjIsComplement(pObj));
                cout << !Abc_ObjFaninC0(pObj) << !Abc_ObjFaninC1(pObj) << " " << !Abc_ObjIsComplement(pObj) << endl;
            }
            else if (Abc_AigNodeIsConst(pObj)) {
                assert(pObj == Abc_AigConst1(GetNet()));
                cout << " 1" << endl;
            }
            else if (Abc_ObjIsPi(pObj))
                cout << endl;
            else if (Abc_ObjIsPo(pObj)) {
                if (Abc_ObjFaninC0(pObj))
                    cout << "0 1" << endl;
                else
                    cout << "1 1" << endl;
            }
            else
                assert(0);
        }
        else
            assert(0);
    }
    else
        cout << endl;
}


/**
 * @brief Check whether the PI names & PO names of this network are the same as the other network.
 * 
 * @param oth_net   The other network
 * @return bool     whether the PIs & POs of the two networks are the same; true if they are the same
 */
bool NetMan::IsPIOSame(const NetMan& oth_net) const {
    if (this->GetPiNum() != oth_net.GetPiNum())
        return false;
    for (int i = 0; i < this->GetPiNum(); ++i) {
        if (this->GetPiName(i) != oth_net.GetPiName(i))
            return false;
    }
    if (this->GetPoNum() != oth_net.GetPoNum())
        return false;
    for (int i = 0; i < this->GetPoNum(); ++i) {
        if (this->GetPoName(i) != oth_net.GetPoName(i))
            return false;
    }
    return true;
}


void NetMan::WriteBlif(const string& fileName) const {
    cout << "write blif to " << fileName << endl;
    FILE * fp = fopen(fileName.c_str(), "w");
    assert(fp != nullptr);
    fprintf(fp, ".model %s\n", GetNet()->pName);
    // dump inputs
    fprintf(fp, ".inputs ");
    for (int i = 0; i < GetPiNum(); ++i)
        fprintf(fp, "%s ", const_cast<char*>(GetPiName(i).c_str()));
    fprintf(fp, "\n");
    // dump outputs
    fprintf(fp, ".outputs ");
    unordered_map<string, int> isPoPrinted;
    for (int i = 0; i < GetPoNum(); ++i) {
        fprintf(fp, "%s ", const_cast<char*>(GetPoName(i).c_str()));
        isPoPrinted[GetPoName(i)] = 0;
    }
    fprintf(fp, "\n");
    // dump nodes
    auto netType = GetNetType();
    assert(netType == NET_TYPE::SOP || netType == NET_TYPE::GATE);
    for (int id = 0; id < GetIdMaxPlus1(); ++id) {
        auto pObj = GetObj(id);
        if (IsNode(pObj)) {
            fprintf(fp, ".names ");
            for (int i = 0; i < GetFaninNum(pObj); ++i)
                fprintf(fp, "%s ", const_cast<char*>(GetName(GetFanin(pObj, i)).c_str()));
            fprintf(fp, "%s\n", const_cast<char*>(GetName(pObj).c_str()));
            if (netType == NET_TYPE::SOP)
                fprintf(fp, "%s\n", static_cast <char *> (pObj->pData));
            else if (netType == NET_TYPE::GATE) {
                fprintf(fp, "# %s\n", Mio_GateReadName((Mio_Gate_t *)pObj->pData));
                fprintf(fp, "%s\n", Mio_GateReadSop((Mio_Gate_t *)pObj->pData));
            }
            else
                assert(0);
            if (isPoPrinted.count(GetName(pObj)))
                isPoPrinted[GetName(pObj)] = 1;
        }
    }
    for (int id = 0; id < GetPoNum(); ++id) {
        auto pObj = GetPo(id);
        if (isPoPrinted[GetName(pObj)] == 0) {
            fprintf(fp, ".names ");
            fprintf(fp, "%s ", const_cast<char*>(GetName(GetFanin(pObj, 0)).c_str()));
            fprintf(fp, "%s\n", const_cast<char*>(GetName(pObj).c_str()));
            assert(Abc_ObjIsComplement(pObj) == 0);
            fprintf(fp, "1 1\n");
        }
    }

    fprintf(fp, ".end\n");
    fclose(fp);
}


static char * Abc_NtkPrintSop( char * pSop ) {
    static char Buffer[1000];
    char * pGet, * pSet;
    pSet = Buffer;
    for ( pGet = pSop; *pGet; pGet++ )
    {        
        if ( *pGet == '\n' )
        {
            *pSet++ = '\\';
            *pSet++ = 'n';
        }
        else
            *pSet++ = *pGet;
    }
    *(pSet-2) = 0;
    return Buffer;
}


static int Abc_NtkCountLogicNodes( Vec_Ptr_t * vNodes ) {
    AbcObj* pObj;
    int i, Counter = 0;
    Vec_PtrForEachEntry( AbcObj*, vNodes, pObj, i )
    {
        if ( !Abc_ObjIsNode(pObj) )
            continue;
        if ( Abc_ObjFaninNum(pObj) == 0 && Abc_ObjFanoutNum(pObj) == 0 )
            continue;
        Counter ++;
    }
    return Counter;
}


static void Net_WriteDotNtk( Abc_Ntk_t * pNtk, Vec_Ptr_t * vNodes, const char * pFileName, int fGateNames, int fUseReverse ) {
    FILE * pFile;
    AbcObj* pNode, * pFanin;
    char * pSopString;
    int LevelMin, LevelMax, fHasCos, Level, i, k, fHasBdds, fCompl, Prev;
    // int Limit = 500;
    // int Limit = 1500;
    int Limit = 2500;

    assert( Abc_NtkIsStrash(pNtk) || Abc_NtkIsLogic(pNtk) );

    if ( vNodes->nSize < 1 )
    {
        printf( "The set has no nodes. DOT file is not written.\n" );
        return;
    }

    if ( vNodes->nSize > Limit )
    {
        printf( "The set has more than %d nodes. DOT file is not written.\n", Limit );
        return;
    }

    // start the stream
    if ( (pFile = fopen( pFileName, "w" )) == NULL )
    {
        fprintf( stdout, "Cannot open the intermediate file \"%s\".\n", pFileName );
        return;
    }

    // transform logic functions from BDD to SOP
    if ( (fHasBdds = Abc_NtkIsBddLogic(pNtk)) )
    {
        if ( !Abc_NtkBddToSop(pNtk, -1, ABC_INFINITY, 1) )
        {
            printf( "Io_WriteDotNtk(): Converting to SOPs has failed.\n" );
            return;
        }
    }

    // get the levels of nodes
    LevelMax = Abc_NtkLevel( pNtk );
    if ( fUseReverse )
    {
        LevelMin = Abc_NtkLevelReverse( pNtk );
        assert( LevelMax == LevelMin );
        Vec_PtrForEachEntry( AbcObj*, vNodes, pNode, i )
            if ( Abc_ObjIsNode(pNode) )
                pNode->Level = LevelMax - pNode->Level + 1;
    }

    // find the largest and the smallest levels
    LevelMin = 10000;
    LevelMax = -1;
    fHasCos  = 0;
    Vec_PtrForEachEntry( AbcObj*, vNodes, pNode, i )
    {
        if ( Abc_ObjIsCo(pNode) )
        {
            fHasCos = 1;
            continue;
        }
        if ( LevelMin > (int)pNode->Level )
            LevelMin = pNode->Level;
        if ( LevelMax < (int)pNode->Level )
            LevelMax = pNode->Level;
    }

    // set the level of the CO nodes
    if ( fHasCos )
    {
        LevelMax++;
        Vec_PtrForEachEntry( AbcObj*, vNodes, pNode, i )
        {
            if ( Abc_ObjIsCo(pNode) )
                pNode->Level = LevelMax;
        }
    }

    // write the DOT header
    fprintf( pFile, "# %s\n",  "Network structure generated by ABC" );
    fprintf( pFile, "\n" );
    fprintf( pFile, "digraph network {\n" );
    fprintf( pFile, "size = \"7.5,10\";\n" );
//    fprintf( pFile, "size = \"10,8.5\";\n" );
//    fprintf( pFile, "size = \"14,11\";\n" );
//    fprintf( pFile, "page = \"8,11\";\n" );
//  fprintf( pFile, "ranksep = 0.5;\n" );
//  fprintf( pFile, "nodesep = 0.5;\n" );
    fprintf( pFile, "center = true;\n" );
//    fprintf( pFile, "orientation = landscape;\n" );
//  fprintf( pFile, "edge [fontsize = 10];\n" );
//  fprintf( pFile, "edge [dir = none];\n" );
    fprintf( pFile, "edge [dir = back];\n" );
    fprintf( pFile, "\n" );

    // labels on the left of the picture
    fprintf( pFile, "{\n" );
    fprintf( pFile, "  node [shape = plaintext];\n" );
    fprintf( pFile, "  edge [style = invis];\n" );
    fprintf( pFile, "  LevelTitle1 [label=\"\"];\n" );
    fprintf( pFile, "  LevelTitle2 [label=\"\"];\n" );
    // generate node names with labels
    for ( Level = LevelMax; Level >= LevelMin; Level-- )
    {
        // the visible node name
        fprintf( pFile, "  Level%d", Level );
        fprintf( pFile, " [label = " );
        // label name
        fprintf( pFile, "\"" );
        fprintf( pFile, "\"" );
        fprintf( pFile, "];\n" );
    }

    // genetate the sequence of visible/invisible nodes to mark levels
    fprintf( pFile, "  LevelTitle1 ->  LevelTitle2 ->" );
    for ( Level = LevelMax; Level >= LevelMin; Level-- )
    {
        // the visible node name
        fprintf( pFile, "  Level%d",  Level );
        // the connector
        if ( Level != LevelMin )
            fprintf( pFile, " ->" );
        else
            fprintf( pFile, ";" );
    }
    fprintf( pFile, "\n" );
    fprintf( pFile, "}" );
    fprintf( pFile, "\n" );
    fprintf( pFile, "\n" );

    // generate title box on top
    fprintf( pFile, "{\n" );
    fprintf( pFile, "  rank = same;\n" );
    fprintf( pFile, "  LevelTitle1;\n" );
    fprintf( pFile, "  title1 [shape=plaintext,\n" );
    fprintf( pFile, "          fontsize=20,\n" );
    fprintf( pFile, "          fontname = \"Times-Roman\",\n" );
    fprintf( pFile, "          label=\"" );
    fprintf( pFile, "%s", "Network structure visualized by ABC" );
    fprintf( pFile, "\\n" );
    fprintf( pFile, "Benchmark \\\"%s\\\". ", pNtk->pName );
    fprintf( pFile, "Time was %s. ",  Extra_TimeStamp() );
    fprintf( pFile, "\"\n" );
    fprintf( pFile, "         ];\n" );
    fprintf( pFile, "}" );
    fprintf( pFile, "\n" );
    fprintf( pFile, "\n" );

    // generate statistics box
    fprintf( pFile, "{\n" );
    fprintf( pFile, "  rank = same;\n" );
    fprintf( pFile, "  LevelTitle2;\n" );
    fprintf( pFile, "  title2 [shape=plaintext,\n" );
    fprintf( pFile, "          fontsize=18,\n" );
    fprintf( pFile, "          fontname = \"Times-Roman\",\n" );
    fprintf( pFile, "          label=\"" );
    if ( Abc_NtkObjNum(pNtk) == Vec_PtrSize(vNodes) )
        fprintf( pFile, "The network contains %d logic nodes and %d latches.", Abc_NtkNodeNum(pNtk), Abc_NtkLatchNum(pNtk) );
    else
        fprintf( pFile, "The set contains %d logic nodes and spans %d levels.", Abc_NtkCountLogicNodes(vNodes), LevelMax - LevelMin + 1 );
    fprintf( pFile, "\\n" );
    fprintf( pFile, "\"\n" );
    fprintf( pFile, "         ];\n" );
    fprintf( pFile, "}" );
    fprintf( pFile, "\n" );
    fprintf( pFile, "\n" );

    // generate the POs
    if ( fHasCos )
    {
        fprintf( pFile, "{\n" );
        fprintf( pFile, "  rank = same;\n" );
        // the labeling node of this level
        fprintf( pFile, "  Level%d;\n",  LevelMax );
        // generate the PO nodes
        Vec_PtrForEachEntry( AbcObj*, vNodes, pNode, i )
        {
            if ( !Abc_ObjIsCo(pNode) )
                continue;
            {
                fprintf( pFile, "  Node%d [label = \"%s%s\"", 
                    pNode->Id, 
                    (Abc_ObjIsBi(pNode)? Abc_ObjName(Abc_ObjFanout0(pNode)):Abc_ObjName(pNode)), 
                    (Abc_ObjIsBi(pNode)? "_in":"") );
            }
            fprintf( pFile, ", shape = %s", (Abc_ObjIsBi(pNode)? "box":"invtriangle") );
            {
                if ( pNode->fMarkB )
                    fprintf( pFile, ", style = filled" );
            }
            fprintf( pFile, ", color = coral, fillcolor = coral" );
            fprintf( pFile, "];\n" );
        }
        fprintf( pFile, "}" );
        fprintf( pFile, "\n" );
        fprintf( pFile, "\n" );
    }

    // generate nodes of each rank
    for ( Level = LevelMax - fHasCos; Level >= LevelMin && Level > 0; Level-- )
    {
        fprintf( pFile, "{\n" );
        fprintf( pFile, "  rank = same;\n" );
        // the labeling node of this level
        fprintf( pFile, "  Level%d;\n",  Level );
        Vec_PtrForEachEntry( AbcObj*, vNodes, pNode, i )
        {
            if ( (int)pNode->Level != Level )
                continue;
            if ( Abc_ObjFaninNum(pNode) == 0 )
                continue;

            if ( Abc_NtkIsStrash(pNtk) )
                pSopString = const_cast<char*>(std::string("").c_str());
            else if ( Abc_NtkHasMapping(pNtk) && fGateNames )
                pSopString = Mio_GateReadName((Mio_Gate_t *)pNode->pData);
            else if ( Abc_NtkHasMapping(pNtk) )
                pSopString = Abc_NtkPrintSop(Mio_GateReadSop((Mio_Gate_t *)pNode->pData));
            else
                pSopString = Abc_NtkPrintSop((char *)pNode->pData);
            fprintf( pFile, "  Node%d [label = \"%s(%d)\\n%s\"", pNode->Id, Abc_ObjName(pNode), pNode->Id, pSopString );

            fprintf( pFile, ", shape = ellipse" );
            // cout << "Node " << pNode << ", fMarkA = " << pNode->fMarkA << ", fMarkB = " << pNode->fMarkB << ", fMarkC = " << pNode->fMarkC << endl;
            {
                if ( pNode->fMarkA )
                    fprintf( pFile, ", style = filled, color = hotpink" );
                else if ( pNode->fMarkB )
                    fprintf( pFile, ", style = filled, color = lightcoral" );
                else if ( pNode->fMarkC )
                    fprintf( pFile, ", style = filled, color = cyan3" );
            }
            fprintf( pFile, "];\n" );
        }
        fprintf( pFile, "}" );
        fprintf( pFile, "\n" );
        fprintf( pFile, "\n" );
    }

    // generate the PI nodes if any
    if ( LevelMin == 0 )
    {
        fprintf( pFile, "{\n" );
        fprintf( pFile, "  rank = same;\n" );
        // the labeling node of this level
        fprintf( pFile, "  Level%d;\n",  LevelMin );
        // generate the PI nodes
        Vec_PtrForEachEntry( AbcObj*, vNodes, pNode, i )
        {
            if ( !Abc_ObjIsCi(pNode) )
            {
                // check if the costant node is present
                if ( Abc_ObjFaninNum(pNode) == 0 && Abc_ObjFanoutNum(pNode) > 0 )
                {
                    fprintf( pFile, "  Node%d [label = \"Const%d(%d)\"", pNode->Id, Abc_NtkIsStrash(pNode->pNtk) || Abc_NodeIsConst1(pNode), pNode->Id );
                    fprintf( pFile, ", shape = ellipse" );
                    if ( pNode->fMarkB )
                        fprintf( pFile, ", style = filled" );
                    fprintf( pFile, ", color = coral, fillcolor = coral" );
                    fprintf( pFile, "];\n" );
                }
                continue;
            }
            {
                fprintf( pFile, "  Node%d [label = \"%s\n(%d)\"", 
                    pNode->Id, 
                    (Abc_ObjIsBo(pNode)? Abc_ObjName(Abc_ObjFanin0(pNode)):Abc_ObjName(pNode)),
                    pNode->Id );
            }
            fprintf( pFile, ", shape = box" );
            {
                if ( pNode->fMarkA )
                    fprintf( pFile, ", style = filled, color = hotpink" );
                else if ( pNode->fMarkB )
                    fprintf( pFile, ", style = filled, color = lightcoral" );
                else if ( pNode->fMarkC )
                    fprintf( pFile, ", style = filled, color = cyan3" );
            }
            // fprintf( pFile, ", color = coral, fillcolor = coral" );
            fprintf( pFile, "];\n" );
        }
        fprintf( pFile, "}" );
        fprintf( pFile, "\n" );
        fprintf( pFile, "\n" );
    }

    // generate invisible edges from the square down
    fprintf( pFile, "title1 -> title2 [style = invis];\n" );
    Vec_PtrForEachEntry( AbcObj*, vNodes, pNode, i )
    {
        if ( (int)pNode->Level != LevelMax )
            continue;
        fprintf( pFile, "title2 -> Node%d [style = invis];\n", pNode->Id );
    }
    // generate invisible edges among the COs
    Prev = -1;
    Vec_PtrForEachEntry( AbcObj*, vNodes, pNode, i )
    {
        if ( (int)pNode->Level != LevelMax )
            continue;
        if ( !Abc_ObjIsPo(pNode) )
            continue;
        if ( Prev >= 0 )
            fprintf( pFile, "Node%d -> Node%d [style = invis];\n", Prev, pNode->Id );
        Prev = pNode->Id;
    }

    // generate edges
    Vec_PtrForEachEntry( AbcObj*, vNodes, pNode, i )
    {
        if ( Abc_ObjIsLatch(pNode) )
            continue;
        Abc_ObjForEachFanin( pNode, pFanin, k )
        {
            if ( Abc_ObjIsLatch(pFanin) )
                continue;
            fCompl = 0;
            if ( Abc_NtkIsStrash(pNtk) )
                fCompl = Abc_ObjFaninC(pNode, k);
            // generate the edge from this node to the next
            fprintf( pFile, "Node%d",  pNode->Id );
            fprintf( pFile, " -> " );
            fprintf( pFile, "Node%d",  pFanin->Id );
            fprintf( pFile, " [style = %s", fCompl? "dotted" : "solid" );
//            fprintf( pFile, ", label = \"%c\"", 'a' + k );
            fprintf( pFile, "]" );
            fprintf( pFile, ";\n" );
        }
    }

    fprintf( pFile, "}" );
    fprintf( pFile, "\n" );
    fprintf( pFile, "\n" );
    fclose( pFile );

    // convert the network back into BDDs if this is how it was
    if ( fHasBdds )
        Abc_NtkSopToBdd(pNtk);
}


void NetMan::WriteDot(const string& fileName) const {
    fmt::print("write dot to {}\n", fileName);
    Vec_Ptr_t * vNodes;
    vNodes = Abc_NtkCollectObjects( pNtk );
    Net_WriteDotNtk( pNtk, vNodes, fileName.c_str(), 0, 0 );
    Vec_PtrFree( vNodes );
}


int NetMan::GetNodeMffcSize(AbcObj* pNode) const {
    assert(IsNode(pNode));
    Vec_Ptr_t * vCone = Vec_PtrAlloc( 100 );
    Abc_NodeDeref_rec(pNode);
    Abc_NodeMffcConeSupp(pNode, vCone, nullptr);
    Abc_NodeRef_rec( pNode );
    int ret = Vec_PtrSize(vCone);
    Vec_PtrFree(vCone);
    return ret;
}


AbcObjVect NetMan::GetNodeMffc(AbcObj* pNode) const {
    assert(IsNode(pNode));
    Vec_Ptr_t * vCone = Vec_PtrAlloc( 100 );
    Abc_NodeDeref_rec(pNode);
    Abc_NodeMffcConeSupp(pNode, vCone, nullptr);
    Abc_NodeRef_rec( pNode );
    AbcObjVect mffc;
    mffc.reserve(Vec_PtrSize(vCone));
    AbcObj* pObj = nullptr;
    int i = 0;
    Vec_PtrForEachEntry(AbcObj*, vCone, pObj, i)
        mffc.emplace_back(pObj);
    Vec_PtrFree(vCone);
    return mffc;
}


static inline int Vec_IntFindFrom(Vec_Int_t * p, int Entry, int start) {
    int i = 0;
    for ( i = start; i < p->nSize; i++ )
        if ( p->pArray[i] == Entry )
            return i;
    assert(0);
    return -1;
}


/**
 * @brief Transfer fanouts from pNodeFrom to pNodeTo
 * @brief Auxiliary function for TempRepl_v2
 * 
 * @param pFanouts    the targeted fanouts of pNodeFrom
 * @param pNodeFrom   the from node
 * @param pNodeTo     the to node
 */
static void Abc_ObjTransferFanout_v2(const AbcObjVect& pFanouts, AbcObj* pNodeFrom, AbcObj* pNodeTo)  {
    assert( !Abc_ObjIsComplement(pNodeFrom) );
    assert( !Abc_ObjIsComplement(pNodeTo) );
    assert( !Abc_ObjIsPo(pNodeFrom) && !Abc_ObjIsPo(pNodeTo) );
    assert( pNodeFrom->pNtk == pNodeTo->pNtk );
    assert( pNodeFrom != pNodeTo );
    assert( !Abc_ObjIsNode(pNodeFrom) || Abc_ObjFanoutNum(pNodeFrom) > 0 );
    // get the fanouts of the old node
    int nFanoutsOld = Abc_ObjFanoutNum(pNodeTo);
    // patch the fanin of each of them
    for (int i = 0; i < static_cast<int>(pFanouts.size()); ++i)
        Abc_ObjPatchFanin(pFanouts[i], pNodeFrom, pNodeTo);
    // assert( Abc_ObjFanoutNum(pNodeFrom) == 0 );
    assert( Abc_ObjFanoutNum(pNodeTo) == nFanoutsOld + static_cast<int>(pFanouts.size()) );
}


/**
 * @brief Temporarily replace pTS with pSS
 * 
 * @param pTS          the target node to be replaced
 * @param pSS          the substitution node
 * @param replTrace    the replacement trace used for recovery; {pTS, pSS, fanout0, iFanin0, fanout1, iFanin1, ...}; pTS is the iFanin-k-th fanin of pTS's fanout-k
 * @param fVerb       whether to print verbose information
 * @retval replTrace   the replacement trace used for recovery
 * @return void
 */
void NetMan::TempRepl_v2(int tsId, int ssId, IntVect& replTrace, bool fVerb) {
    // check
    auto pTS = GetObj(tsId), pSS = GetObj(ssId);
    if (fVerb)
        fmt::print("temporarily replace {} with {}\n", *pTS, *pSS);
    assert(tsId != ssId);
    assert(abc::Abc_ObjFanoutNum(pTS));
    // collect fanouts of pTS
    AbcObjVect pFanouts;
    pFanouts.reserve(abc::Abc_ObjFanoutNum(pTS));
    AbcObj* pFanout = nullptr;
    int i = 0;
    Abc_ObjForEachFanout(pTS, pFanout, i) {
        // special case: pFanout is pSS, the substitution node
        if (pFanout == pSS) {
            if (fVerb)
                fmt::print("skip fanout {} of {}\n", *pFanout, *pTS);
            continue;
        }
        pFanouts.emplace_back(pFanout);
    }
    // record (fanout, iFanin) pairs: for each pFanout of pTS, iFanin is the index of pTS in pFanout's fanin list
    replTrace.clear();
    replTrace.emplace_back(pTS->Id);
    replTrace.emplace_back(pSS->Id);
    set<IntPair> foIFaninPair;
    for (const auto& pFanout: pFanouts) {
        replTrace.emplace_back(pFanout->Id);
        int start = 0;
        int iFanin = Vec_IntFindFrom(&pFanout->vFanins, pTS->Id, start);
        // special case: pTS appears multiple times in pFanout's fanin list
        while (foIFaninPair.count(pair(pFanout->Id, iFanin)))
            iFanin = Vec_IntFindFrom(&pFanout->vFanins, pTS->Id, iFanin + 1);
        replTrace.emplace_back(iFanin);
        foIFaninPair.emplace(pair(pFanout->Id, iFanin));
    }
    // transfer fanouts
    Abc_ObjTransferFanout_v2(pFanouts, pTS, pSS);
}


// /**
//  * @brief Recover the replacement from the trace
//  * 
//  * @param replTrace  the replacement trace; the replacement trace used to recover the network; {pTS, pSS, fanout0, iFanin0, fanout1, iFanin1, ...}; pTS is the iFanin-k-th fanin of pTS's fanout-k
//  * @param fVerb     whether to print the recovery process
//  * @return void
//  */
// void NetMan::Recov(const IntVect& replTrace, bool fVerb) {
//     assert(replTrace.size() > 2);
//     assert((replTrace.size() & 1) == 0);
//     auto pTS = GetObj(replTrace[0]), pSS = GetObj(replTrace[1]);
//     if (fVerb) 
//         fmt::print("recover [{}, {}]\n", *pTS, *pSS);
//     for (int i = 1; i < replTrace.size() / 2; ++i) {
//         auto pFanout = GetObj(replTrace[i * 2]);
//         auto iFanin = replTrace[i * 2 + 1];
//         PatchFanin(pFanout, iFanin, pSS, pTS);
//     }
// }


/**
 * @brief Recover the replacement from the trace
 * 
 * @param replTrace  the replacement trace used to recover the network; {pTS, pSS, fanout0, iFanin0, fanout1, iFanin1, ..., -1, delObj0, delObj1, ...}; pTS is the iFanin-k-th fanin of pTS's fanout-k; delObj is the object to be deleted
 * @param fVerb     whether to print the recovery process
 * @return void
 */
void NetMan::Recov_v2(const IntVect& replTrace, bool fVerb) {
    // prepare
    assert(replTrace.size() > 2);
    auto itDelObj = std::ranges::find(replTrace, -1);
    int patchFaninEnd = 0;
    if (itDelObj != replTrace.end()) // have objects to be deleted
        patchFaninEnd = itDelObj - replTrace.begin();
    else // no objects to be deleted
        patchFaninEnd = replTrace.size();
    assert((patchFaninEnd & 1) == 0);
    // get pTS and pSS
    auto pTS = GetObj(replTrace[0]), pSS = GetObj(replTrace[1]);
    if (fVerb) 
        fmt::print("recover [pTS={}, pSS={}]: ", *pTS, *pSS);
    // patch fanins
    for (int i = 1; i < patchFaninEnd / 2; ++i) {
        auto pFanout = GetObj(replTrace[i * 2]);
        auto iFanin = replTrace[i * 2 + 1];
        PatchFanin(pFanout, iFanin, pSS, pTS);
        if (fVerb)
            fmt::print("patch [fanout={}, iFanin={}, pSS={}, pTS={}], ", *pFanout, iFanin, *pSS, *pTS);
    }
    // delete objects
    if (itDelObj == replTrace.end()) { // no objects to be deleted
        if (fVerb)
            fmt::print("\n");
        return;
    }
    assert(*itDelObj == -1);
    for (auto it = itDelObj + 1; it != replTrace.end(); ++it) {
        auto pObj = GetObj(*it);
        if (fVerb)
            fmt::print("delete {}, ", *pObj);
        DelObj(pObj);
    }
    if (fVerb)
        fmt::print("\n");
}


static inline int Vec_IntFindRev( Vec_Int_t * p, int Entry ) {
    int i;
    // for ( i = 0; i < p->nSize; i++ )
    for (i = p->nSize - 1; i >= 0; --i)
        if ( p->pArray[i] == Entry )
            return i;
    return -1;
}


static inline int Vec_IntRemoveRev( Vec_Int_t * p, int Entry ) {
    int i;
    // for ( i = 0; i < p->nSize; i++ )
    for (i = p->nSize - 1; i >= 0; --i)
        if ( p->pArray[i] == Entry )
            break;
    if ( i == p->nSize )
        return 0;
    assert( i < p->nSize );
    for ( i++; i < p->nSize; i++ )
        p->pArray[i-1] = p->pArray[i];
    p->nSize--;
    return 1;
}


static inline void Vec_IntPushMem( Mem_Step_t * pMemMan, Vec_Int_t * p, int Entry ) {
    if ( p->nSize == p->nCap )
    {
        int * pArray;
        int i;

        if ( p->nSize == 0 )
            p->nCap = 1;
        if ( pMemMan )
            pArray = (int *)Mem_StepEntryFetch( pMemMan, p->nCap * 8 );
        else
            pArray = ABC_ALLOC( int, p->nCap * 2 );
        if ( p->pArray )
        {
            for ( i = 0; i < p->nSize; i++ )
                pArray[i] = p->pArray[i];
            if ( pMemMan )
                Mem_StepEntryRecycle( pMemMan, (char *)p->pArray, p->nCap * 4 );
            else
                ABC_FREE( p->pArray );
        }
        p->nCap *= 2;
        p->pArray = pArray;
    }
    p->pArray[p->nSize++] = Entry;
}


void NetMan::PatchFanin( AbcObj* pObj, int iFanin, AbcObj* pFaninOld, AbcObj* pFaninNew ) {
    AbcObj* pFaninNewR = Abc_ObjRegular(pFaninNew);
    assert( !Abc_ObjIsComplement(pObj) );
    assert( !Abc_ObjIsComplement(pFaninOld) );
    assert( pFaninOld != pFaninNewR );
    assert( pObj->pNtk == pFaninOld->pNtk );
    assert( pObj->pNtk == pFaninNewR->pNtk );
    assert( abc::Abc_ObjFanin(pObj, iFanin) == pFaninOld );

    // remember the attributes of the old fanin
    Vec_IntWriteEntry( &pObj->vFanins, iFanin, pFaninNewR->Id );
    if ( Abc_ObjIsComplement(pFaninNew) )
        Abc_ObjXorFaninC( pObj, iFanin );

    // update the fanout of the fanin
    if ( !Vec_IntRemoveRev( &pFaninOld->vFanouts, pObj->Id ) ) {
        printf( "Node %s is not among", Abc_ObjName(pObj) );
        printf( " the fanouts of its old fanin %s...\n", Abc_ObjName(pFaninOld) );
    }
    Vec_IntPushMem( pObj->pNtk->pMmStep, &pFaninNewR->vFanouts, pObj->Id );
}


void NetMan::Trunc(int truncBit) {
    cout << "***** truncate " << truncBit << " bits" << endl;
    // truncation
    auto consts = CreateConstsIfNotExist();
    assert(truncBit <= GetPoNum());
    for (int poId = 0; poId < truncBit; ++poId) {
        auto pPo = GetPo(poId);
        assert(GetFaninNum(pPo) == 1);
        auto pDriv = GetFanin(pPo, 0);
        Abc_ObjPatchFanin(pPo, pDriv, GetObj(consts.first));
    }
    // clean up
    CleanUp();
}


// bool NetMan::ProcHalfAndFullAdd() {
//     // special processing for half/full adder
//     if (GetNetType() != NET_TYPE::GATE)
//         return false;
//     bool isUpd = false;
//     int idMaxPlus1 = GetIdMaxPlus1();
//     for (int nodeId = 0; nodeId < idMaxPlus1; ++nodeId) {
//         if (!IsNode(nodeId)) continue;
//         auto pNode = GetObj(nodeId);
//         if (GetGateName(pNode).find("HA1") != -1) {
//             if (GetTwinNode(pNode) == nullptr) {
//                 cout << "cannot find twin for "; PrintObj(pNode, true);
//                 auto sop = string(Mio_GateReadSop((Mio_Gate_t *)pNode->pData));
//                 auto pLib = (Mio_Library_t *)Abc_FrameReadLibGen();
//                 auto pNewNode = Abc_NtkCreateNode(GetNet());
//                 for (int faninId = 0; faninId < GetFaninNum(pNode); ++faninId)
//                     Abc_ObjAddFanin(pNewNode, GetFanin(pNode, faninId));
//                 AbcObj* pSub = nullptr;
//                 if (sop == "11 1\n") { // CO=A B
//                     auto pGate = Mio_LibraryReadGateByName(pLib, "CKAN2D1BWP7T30P140HVT", nullptr);
//                     assert(pGate != nullptr);
//                     pNewNode->pData = pGate;
//                     pSub = pNewNode;
//                 }
//                 else if (sop == "10 1\n01 1\n" || sop == "01 1\n10 1\n") { // S=A^B
//                     auto pGate = Mio_LibraryReadGateByName(pLib, "XOR2D0BWP7T30P140HVT", nullptr);
//                     assert(pGate != nullptr);
//                     pNewNode->pData = pGate;
//                     pSub = pNewNode;
//                 }
//                 else {
//                     cout << sop;
//                     assert(0);
//                 }
//                 assert(pSub != nullptr);
//                 cout << "replace " << pNode << " with new node " << pSub << endl;
//                 TransfFanout(pNode, pSub);
//                 isUpd = true;
//             }
//         }
//         else if (GetGateName(pNode).find("FA1") != -1) {
//             if (GetTwinNode(pNode) == nullptr) {
//                 cout << "cannot find twin for "; PrintObj(pNode, true);
//                 auto sop = string(Mio_GateReadSop((Mio_Gate_t *)pNode->pData));
//                 auto pLib = (Mio_Library_t *)Abc_FrameReadLibGen();
//                 auto pNewNode = Abc_NtkCreateNode(GetNet());
//                 for (int faninId = 0; faninId < GetFaninNum(pNode); ++faninId)
//                     Abc_ObjAddFanin(pNewNode, GetFanin(pNode, faninId));
//                 AbcObj* pSub = nullptr;
//                 if (sop == "1-1 1\n-11 1\n11- 1\n") { // CO=A B+B CI+A CI
//                     auto pGate = Mio_LibraryReadGateByName(pLib, "MAOI222D0BWP7T30P140HVT", nullptr);
//                     assert(pGate != nullptr);
//                     pNewNode->pData = pGate;
//                     pSub = Abc_NtkCreateNodeInv(GetNet(), pNewNode);
//                 }
//                 else if (sop == "100 1\n010 1\n111 1\n001 1\n") { // S=A^B^CI
//                     auto pGate = Mio_LibraryReadGateByName(pLib, "XOR3D1BWP7T30P140HVT", nullptr);
//                     assert(pGate != nullptr);
//                     pNewNode->pData = pGate;
//                     pSub = pNewNode;
//                 }
//                 else {
//                     cout << sop;
//                     assert(0);
//                 }
//                 assert(pSub != nullptr);
//                 cout << "replace " << pNode << " with new node " << pSub << endl;
//                 TransfFanout(pNode, pSub);
//                 isUpd = true;
//             }
//         }
//     }
//     return isUpd;
// }


// void NetMan::ProcHalfAndFullAddNew() {
//     // special processing for half/full adder
//     if (GetNetType() != NET_TYPE::GATE)
//         return;
//     unordered_set <int> vis;
//     IntVect targNodes;
//     for (int iNode = 0; iNode < GetIdMaxPlus1(); ++iNode) {
//         if (!IsNode(iNode))
//             continue;
//         if (vis.count(iNode))
//             continue;
//         vis.emplace(iNode);
//         auto pNode = GetObj(iNode);
//         auto pTwin = GetTwinNode(pNode);
//         if (pTwin == nullptr)
//             continue;
//         vis.emplace(pTwin->Id);
//         if (GetGateName(pNode).find("HA1") != -1)
//             targNodes.emplace_back(iNode);
//         else if (GetGateName(pNode).find("FA1") != -1)
//             targNodes.emplace_back(iNode);
//         else
//             assert(0);
//     }

//     for (int targId: targNodes) {
//         auto pNode = GetObj(targId);
//         auto pTwin = GetTwinNode(pNode);
//         assert(pTwin != nullptr);
//         if (GetGateName(pNode).find("HA1") != -1) {
//             // print
//             // PrintObj(pNode, true); 
//             auto sop0 = string(Mio_GateReadSop((Mio_Gate_t *)pNode->pData));
//             // cout << sop0;
//             // PrintObj(pTwin, true);
//             auto sop1 = string(Mio_GateReadSop((Mio_Gate_t *)pTwin->pData));
//             // cout << sop1;
//             // cout << endl;

//             // pNode sop0 S, pTwin sop1 CO
//             assert(sop0 == "10 1\n01 1\n" && sop1 == "11 1\n");
//             AbcObjVect fanins;
//             int nFanin = GetFaninNum(pNode);
//             assert(nFanin == GetFaninNum(pTwin) && nFanin == 2);
//             for (int iFanin = 0; iFanin < nFanin; ++iFanin) {
//                 assert(GetFanin(pNode, iFanin) == GetFanin(pTwin, iFanin));
//                 fanins.emplace_back(GetFanin(pNode, iFanin));
//             }

//             // create gates
//             auto pNodeCo = CreateGate(AbcObjVect ({fanins[0], fanins[1]}), "CKAN2D1BWP7T30P140HVT");
//             auto pNodeN6 = CreateGate(AbcObjVect ({fanins[0], fanins[1]}), "NR2D0BWP7T30P140HVT");
//             auto pNodeS = CreateGate(AbcObjVect ({pNodeCo, pNodeN6}), "NR2D0BWP7T30P140HVT");

//             // cout << "replace " << pTwin << " with new node " << pNodeCo << endl;
//             TransfFanout(pTwin, pNodeCo);
//             // cout << "replace " << pNode << " with new node " << pNodeS << endl;
//             TransfFanout(pNode, pNodeS);
//         }
//         else if (GetGateName(pNode).find("FA1") != -1) {
//             // print
//             // PrintObj(pNode, true); 
//             auto sop0 = string(Mio_GateReadSop((Mio_Gate_t *)pNode->pData));
//             // cout << sop0;
//             // PrintObj(pTwin, true);
//             auto sop1 = string(Mio_GateReadSop((Mio_Gate_t *)pTwin->pData));
//             // cout << sop1;
//             // cout << endl;

//             // pNode sop0 S, pTwin sop1 CO
//             assert(sop0 == "100 1\n010 1\n111 1\n001 1\n" && sop1 == "1-1 1\n-11 1\n11- 1\n");
//             AbcObjVect fanins;
//             int nFanin = GetFaninNum(pNode);
//             assert(nFanin == GetFaninNum(pTwin) && nFanin == 3);
//             for (int iFanin = 0; iFanin < nFanin; ++iFanin) {
//                 assert(GetFanin(pNode, iFanin) == GetFanin(pTwin, iFanin));
//                 fanins.emplace_back(GetFanin(pNode, iFanin));
//             }

//             // create gates
//             auto pNodeN6 = CreateGate(AbcObjVect ({fanins[1], fanins[2], fanins[1], fanins[2]}), "MOAI22D0BWP7T30P140HVT");
//             auto pNodeS = CreateGate(AbcObjVect ({fanins[0], pNodeN6, fanins[0], pNodeN6}), "MOAI22D0BWP7T30P140HVT");
//             auto pNodeCo = CreateGate(AbcObjVect ({fanins[1], fanins[2], fanins[0], pNodeN6}), "OA22D0BWP7T30P140HVT");

//             // cout << "replace " << pTwin << " with new node " << pNodeCo << endl;
//             TransfFanout(pTwin, pNodeCo);
//             // cout << "replace " << pNode << " with new node " << pNodeS << endl;
//             TransfFanout(pNode, pNodeS);
//         }
//     }
//     CleanUp();
// }


AbcObj* NetMan::CreateNode(const AbcObjVect& pFanins, const std::string& sop) {
    auto pNewNode = abc::Abc_NtkCreateNode(GetNet());
    for (const auto& pFanin: pFanins)
        Abc_ObjAddFanin(pNewNode, pFanin);
    assert(GetNetType() == NET_TYPE::SOP);
    pNewNode->pData = Abc_SopRegister((Mem_Flex_t *)GetNet()->pManFunc, sop.c_str());
    return pNewNode;
}


/**
 * @brief Create a new node with a SOP-style function
 * 
 * @param faninIds   fanin ids
 * @param sop        sum-of-product representing the function
 * @return int       id of the created node 
 */
int NetMan::CreateNode(const IntVect& faninIds, const std::string& sop) {
    auto pNewNode = abc::Abc_NtkCreateNode(GetNet());
    for (const auto& faninId: faninIds)
        Abc_ObjAddFanin(pNewNode, GetObj(faninId));
    assert(GetNetType() == NET_TYPE::SOP);
    pNewNode->pData = Abc_SopRegister((Mem_Flex_t *)GetNet()->pManFunc, sop.c_str());
    return pNewNode->Id;
}


/**
 * @brief Create AIG-style node in a SOP network
 * 
 * @param faninIds  fanin node ids
 * @param sop       the sum-of-product expression of the created node
 * @return int      the id of the created node
 */
int NetMan::CreateAIGStyleNodes(const IntVect& faninIds, const std::string& sop) {
    assert(GetNetType() == NET_TYPE::SOP);
    if (sop == "01 1\n10 1\n" || sop == "10 1\n01 1\n" || sop == "00 0\n11 0\n" || sop == "11 0\n00 0\n") { // xor
        auto and0 = CreateNode(faninIds, "01 1\n");
        auto and1 = CreateNode(faninIds, "10 1\n");
        cout << "create xor" << endl;
        return CreateNode({and0, and1}, "00 0\n");
    }
    else if (sop == "01 0\n10 0\n" || sop == "10 0\n01 0\n" || sop == "00 1\n11 1\n" || sop == "11 1\n00 1\n") { // xnor
        auto and0 = CreateNode(faninIds, "00 1\n");
        auto and1 = CreateNode(faninIds, "11 1\n");
        cout << "create xnor" << endl;
        return CreateNode({and0, and1}, "00 0\n");
    }
    else if (faninIds.size() == 3) {
        auto v0 = faninIds[0], v1 = faninIds[1], v2 = faninIds[2];
        auto pSop = const_cast<char*>(sop.c_str());
        bool isComplement = abc::Abc_SopIsComplement(pSop);
        int nCube = abc::Abc_SopGetCubeNum(pSop);
        int nVar = abc::Abc_SopGetVarNum(pSop);
        assert(nVar == 3);
        if (nCube == 1) {
            auto pCube = pSop;
            assert(pCube[0] != '-' && pCube[1] != '-' && pCube[2] != '-');
            auto and0 = CreateNode({v0, v1}, to_string(pCube[0] - '0') + to_string(pCube[1] - '0') + " 1\n");
            return CreateNode({and0, v2}, "1" + to_string(pCube[2] - '0') + (isComplement? " 0\n": " 1\n"));
        }
        else if (nCube == 2) {
            auto pCube0 = pSop;
            vector<int> cv0;
            for (int i = 0; i < nVar; ++i) {
                if (pCube0[i] != '-')
                    cv0.emplace_back(i);
            }
            auto pCube1 = pSop + nVar + 3;
            vector<int> cv1;
            for (int i = 0; i < nVar; ++i) {
                if (pCube1[i] != '-')
                    cv1.emplace_back(i);
            }
            // cout << "cv0.size() = " << cv0.size() << ", cv1.size() = " << cv1.size() << endl;
            if (cv0.size() == 1 && cv1.size() == 2) {
                string phaseVar0((pCube0[cv0[0]] == '0')? "1": "0");
                // cout << "phaseVar0 = " << phaseVar0 << endl;
                auto and0 = CreateNode({faninIds[cv1[0]], faninIds[cv1[1]]}, to_string(pCube1[cv1[0]] - '0') + to_string(pCube1[cv1[1]] - '0') + " 1\n");
                // cout << "sop = " << to_string(pCube1[cv1[0]] - '0') + to_string(pCube1[cv1[1]] - '0') + " 1\n" << endl;
                // cout << "sop = " << "0" + phaseVar0 + (isComplement? " 1\n": " 0\n") << endl;
                return CreateNode({and0, faninIds[cv0[0]]}, "0" + phaseVar0 + (isComplement? " 1\n": " 0\n"));
            }
            else if (cv0.size() == 2 && cv1.size() == 1) {
                string phaseVar1((pCube1[cv1[0]] == '0')? "1": "0");
                auto and0 = CreateNode({faninIds[cv0[0]], faninIds[cv0[1]]}, to_string(pCube0[cv0[0]] - '0') + to_string(pCube0[cv0[1]] - '0') + " 1\n");
                // cout << "sop = " << to_string(pCube0[cv0[0]] - '0') + to_string(pCube0[cv0[1]] - '0') + " 1\n" << endl;
                // cout << "sop = " << "0" + phaseVar1 + (isComplement? " 1\n": " 0\n") << endl;
                return CreateNode({and0, faninIds[cv1[0]]}, "0" + phaseVar1 + (isComplement? " 1\n": " 0\n"));
            }
            else
                assert(0);
        }
        else
            assert(0);
    }
    else {
        auto pNewNode = abc::Abc_NtkCreateNode(GetNet());
        for (const auto & faninId: faninIds)
            Abc_ObjAddFanin(pNewNode, GetObj(faninId));
        pNewNode->pData = Abc_SopRegister((Mem_Flex_t *)GetNet()->pManFunc, sop.c_str());
        return pNewNode->Id;
    }
    return 0;
}


AbcObj* NetMan::CreateGate(AbcObjVect&& fanins, const std::string& gateName) {
    auto pLib = (Mio_Library_t *)Abc_FrameReadLibGen();
    auto pGate = Mio_LibraryReadGateByName(pLib, const_cast <char *> (gateName.c_str()), nullptr);
    assert(pGate != nullptr);
    auto pNewNode = Abc_NtkCreateNode(GetNet());
    for (const auto& fanin: fanins)
        Abc_ObjAddFanin(pNewNode, fanin);
    pNewNode->pData = pGate;
    return pNewNode;
}


AbcObj* NetMan::DupObj(AbcObj* pObj, const char* pSuff) {
    AbcObj* pObjNew;
    // create the new object
    pObjNew = Abc_NtkCreateObj( pNtk, (Abc_ObjType_t)pObj->Type );
    // transfer names of the terminal objects
    Abc_ObjAssignName(pObjNew, Abc_ObjName(pObj), const_cast<char*>(pSuff));
    // copy functionality/names
    if ( Abc_ObjIsNode(pObj) ) // copy the function if functionality is compatible
    {
        if ( pNtk->ntkFunc == pObj->pNtk->ntkFunc ) 
        {
            if ( Abc_NtkIsStrash(pNtk) ) 
            {}
            else if ( Abc_NtkHasSop(pNtk) || Abc_NtkHasBlifMv(pNtk) )
                pObjNew->pData = Abc_SopRegister( (Mem_Flex_t *)pNtk->pManFunc, (char *)pObj->pData );
#ifdef ABC_USE_CUDD
            else if ( Abc_NtkHasBdd(pNtk) )
                pObjNew->pData = Cudd_bddTransfer((DdManager *)pObj->pNtk->pManFunc, (DdManager *)pNtk->pManFunc, (DdNode *)pObj->pData), Cudd_Ref((DdNode *)pObjNew->pData);
#endif
            else if ( Abc_NtkHasAig(pNtk) )
                pObjNew->pData = Hop_Transfer((Hop_Man_t *)pObj->pNtk->pManFunc, (Hop_Man_t *)pNtk->pManFunc, (Hop_Obj_t *)pObj->pData, Abc_ObjFaninNum(pObj));
            else if ( Abc_NtkHasMapping(pNtk) )
                pObjNew->pData = pObj->pData, pNtk->nBarBufs2 += !pObj->pData;
            else assert( 0 );
        }
    }
    else if ( Abc_ObjIsNet(pObj) ) // copy the name
    {
    }
    else if ( Abc_ObjIsLatch(pObj) ) // copy the reset value
        pObjNew->pData = pObj->pData;
    pObjNew->fPersist = pObj->fPersist;
    // transfer HAIG
//    pObjNew->pEquiv = pObj->pEquiv;
    // remember the new node in the old node
    pObj->pCopy = pObjNew;
    return pObjNew;
}


void NetMan::LimFanout() {
    assert(GetNetType() == NET_TYPE::SOP);
    auto nodes = CalcTopoOrdOfIds();
    for (int id: nodes) {
        int nFo = GetFanoutNum(id);
        if (nFo <= 2)
            continue;
        // PrintObj(id, true);
        list<Abc_Obj_t*> dealtPFos;
        for (int iFo = 1; iFo < nFo; ++iFo)
            dealtPFos.emplace_back(GetFanout(id, iFo));
        auto pDealt = GetObj(id);
        while (!dealtPFos.empty()) {
            auto pBuf = CreateBuf(pDealt);
            Rename(pBuf, (GetName(pDealt) + "_b").c_str());
            for (auto pFo: dealtPFos)
                Abc_ObjPatchFanin(pFo, pDealt, pBuf);
            dealtPFos.pop_front();
            pDealt = pBuf;
        }
    }
}


/**
 * @brief Replace the target node by the substitution node + inverter
 * 
 * @param targId  target node id
 * @param subId   substitution node id
 * @return void
 */
void NetMan::ReplaceByComplementedObj(int targId, int subId) {
    assert(GetNetType() == NET_TYPE::SOP);
    auto pTarg = GetObj(targId);
    auto pSub = GetObj(subId);
    assert(pTarg != nullptr && pSub != nullptr);
    // collect fanouts of pTarg
    AbcObjVect fanouts;
    for (int i = 0; i < GetFanoutNum(pTarg); ++i)
        fanouts.emplace_back(GetFanout(pTarg, i));
    // for each pTarg's fanout (marked as fo), collect the fanin index of pTarg; 
    // meanwhile, make sure that pTarg only appears once in fo's fanin list
    vector <int> foIFanin;
    for (const auto& fo: fanouts) {
        vector <int> iFanins;
        for (int i = 0; i < GetFaninNum(fo); ++i) {
            if (GetFanin(fo, i) == pTarg)
                iFanins.emplace_back(i);
        }
        assert(iFanins.size() == 1);
        foIFanin.emplace_back(iFanins[0]);
    }
    // transfer fanouts
    Abc_ObjTransferFanout( pTarg, pSub );
    // fix SOP functions for fanouts
    AbcObj* pInv = nullptr;
    for (int i = 0; i < static_cast<int>(fanouts.size()); ++i) {
        auto fo = fanouts[i];
        auto iFanin = foIFanin[i];
        if (IsObjPo(fo)) {
            if (pInv == nullptr) {
                pInv = Abc_NtkCreateNodeInv( GetNet(), pSub );
                // cout << "create inverter for " << pSub << ":" << pInv << endl;
                fmt::print("create inverter for {}: {}\n", *pSub, *pInv);
            } 
            // because pTarg's fanouts have been transferred to pSub, so pSub is the only fanin of fo
            assert(GetFaninNum(fo) == 1 && GetFanin(fo, 0) == pSub);
            Abc_ObjPatchFanin( fo, pSub, pInv );
        }
        else {
            auto pSop = static_cast <char *> (fo->pData);
            Abc_SopComplementVar( pSop, iFanin );
        }
    }
    // remove redundant nodes
    Abc_NtkDeleteObj_rec( pTarg, 1 );
}


/**
 * @brief Auxiliary function for constant propagation
 * @brief pFanin is a fanin of pNode; set pNode's fanin pFanin to be constant 0 or 1
 * 
 * @param pNode     the fanout node
 * @param pFanin    the fanin node
 * @param fConst0   whether to set the fanin to be constant 0
 */
static void SetConstInput(AbcObj* pNode, AbcObj* pFanin, int fConst0) {
    assert(pNode != nullptr && pFanin != nullptr);
    assert(pNode->pNtk == pFanin->pNtk);
    auto pNtk = pNode->pNtk;
    int iFanin = Vec_IntFind( &pNode->vFanins, pFanin->Id );
    if (iFanin == -1 ) {
        printf( "Node %s should be among", Abc_ObjName(pFanin) );
        printf( " the fanins of node %s...\n", Abc_ObjName(pNode) );
        assert(0);
        return;
    }

    // construct new sop
    string newSop = "";
    auto pOldSop = static_cast<char*>(pNode->pData);
    bool isOldSopComplemented = Abc_SopIsComplement(pOldSop);
    int nVars = Abc_SopGetVarNum(pOldSop);
    assert(iFanin < nVars);
    char * pCube = nullptr;
    Abc_SopForEachCube(pOldSop, nVars, pCube) {
        if ((fConst0 && pCube[iFanin] != '1') || (!fConst0 && pCube[iFanin] != '0')) // consider don't care
        {
            // append pCube to newSop except for the iFanin-th bit
            for (int i = 0; i < nVars; ++i) {
                if (i == iFanin)
                    continue;
                newSop += pCube[i];
            }
            newSop += isOldSopComplemented? " 0\n": " 1\n";
        }
    }
    if (newSop == "")
        newSop = isOldSopComplemented? " 1\n": " 0\n";
    // cout << "old fanins: ";
    // for (int i = 0; i < GetFaninNum(pNode); ++i)
    //     cout << GetFanin(pNode, i) << " ";
    // cout << endl;
    // cout << "old sop: " << pOldSop;

    // remove pNode's fanin: pFanin
    if (newSop == " 1\n" || newSop == " 0\n")
        Abc_ObjRemoveFanins(pNode);
    else
        Abc_ObjDeleteFanin( pNode, pFanin );
    // cout << "pFanin: " << pFanin << endl;
    // cout << "pFanin's fanouts: ";
    // for (int i = 0; i < GetFanoutNum(pFanin); ++i)
    //     cout << GetFanout(pFanin, i) << " ";
    // cout << endl;

    // update pNode's sop
    pNode->pData = Abc_SopRegister((Mem_Flex_t *)pNtk->pManFunc, newSop.c_str());
    // cout << "new fanins: ";
    // for (int i = 0; i < GetFaninNum(pNode); ++i)
    //     cout << GetFanin(pNode, i) << " ";
    // cout << endl;
    // cout << "new sop: " << newSop;
    // cout << endl;
}


/**
 * @brief Constant propagation from a starting node
 * @brief Refer to the function Abc_NtkSweep() in the ABC library
 * 
 * @param startId  the starting node id
 * @param fVerb   whether to print the information
 * @retval pNtk    the network after propagating the constant from the starting node
 * @return void
 */
void NetMan::PropConst(int startId, bool fKeepDanglNodes, bool fVerb) {
    // check & print
    assert(GetNetType() == NET_TYPE::SOP);
    assert(IsConst(startId));
    if (fVerb)
        fmt::print("propagate const {}\n", *GetObj(startId));
    // initialize
    auto vNodes = Vec_PtrAlloc( 100 );
    Vec_PtrPush( vNodes, GetObj(startId) );
    // sweep the nodes
    while ( Vec_PtrSize(vNodes) > 0 ) {
        // get any sweepable node
        auto pNode = (AbcObj*)Vec_PtrPop(vNodes);
        if ( !Abc_ObjIsNode(pNode) )
            continue;
        // get any non-CO fanout of this node
        auto pFanout = Abc_NodeFindNonCoFanout(pNode);
        if ( pFanout == nullptr ) 
            continue;
        assert( Abc_ObjIsNode(pFanout) );
        // transform the function of the fanout
        if ( Abc_ObjFaninNum(pNode) == 0 ) // constant node
            SetConstInput( pFanout, pNode, Abc_NodeIsConst0(pNode) );
        else { // buffer or inverter
            assert( Abc_ObjFaninNum(pNode) == 1 );
            auto pDriver = Abc_ObjFanin0(pNode);
            if ( Abc_NodeIsInv(pNode) )
                Abc_NodeComplementInput( pFanout, pNode );
            Abc_ObjPatchFanin( pFanout, pNode, pDriver );
        }
        // check if the fanout should be added
        if ( Abc_ObjFaninNum(pFanout) < 2 )
            Vec_PtrPush( vNodes, pFanout );
        // check if the node has other fanouts
        if ( Abc_ObjFanoutNum(pNode) > 0 )
            Vec_PtrPush( vNodes, pNode );
        else {
            if (!fKeepDanglNodes)
                Abc_NtkDeleteObj_rec( pNode, 1 );
        }
    }
    Vec_PtrFree( vNodes );
    // // clean up (SetConstInput might generate dangling nodes)
    // if (!fKeepDanglNodes)
    //     CleanUp();
}


/**
 * @brief Constant propagation for the whole network
 * 
 * @param fVerb   whether to print the information
 * @retval pNtk    the network after propagating the constant from the starting node
 * @return void
 */
void NetMan::PropConst(bool fVerb) {
    // merge constant nodes
    MergeConst(fVerb);

    // propagate constant
    auto constIds = GetConstIds(fVerb);
    // fmt::print("constIds: {} {}\n", constIds.first, constIds.second);
    if (constIds.first != -1)
        PropConst(constIds.first, false, fVerb);
    if (constIds.second != -1)
        PropConst(constIds.second, false, fVerb);

    // clean up
    CleanUp();
    
    // check: after constant propagation, constant nodes' fanouts should be all POs
    constIds = GetConstIds(fVerb);
    if (constIds.first != -1) {
        for (int i = 0; i < GetFanoutNum(constIds.first); ++i)
            assert(IsObjPo(GetFanoutId(constIds.first, i)));
    }
    if (constIds.second != -1) {
        for (int i = 0; i < GetFanoutNum(constIds.second); ++i)
            assert(IsObjPo(GetFanoutId(constIds.second, i)));
    }
}


/**
 * @brief Recursively dereference the node pNode and its MFFCs
 * 
 * @param pRoot    the source node of the dereference
 * @param pNode    the current node to be dereferenced
 * @param divSet   the set of nodes should be kept
 * @return int     the number of nodes dereferenced
 */
int NetMan::NodeDeref_rec(AbcObj* pRoot, AbcObj* pNode, unordered_set<AbcObj*>& divSet) const {
    AbcObj* pFanin;
    int i, Counter = 1;
    // if pNode is a PI, then it should be kept
    // if pNode is in divSet, then it should be kept
    // if pNode is a PO driver, and pNode is not the root (the source node of the dereference), then it should be kept
    if ( Abc_ObjIsCi(pNode) || divSet.count(pNode) || (pRoot != pNode && IsPoDriver(pNode)) )
        return 0;
    Abc_ObjForEachFanin( pNode, pFanin, i ) {
        assert( pFanin->vFanouts.nSize > 0 );
        if ( --pFanin->vFanouts.nSize == 0 )
            Counter += NodeDeref_rec(pRoot, pFanin, divSet);
    }
    return Counter;
}


/**
 * @brief Recursively reference the node pNode and its MFFCs
 * @brief The inverse operation of NodeDeref_rec
 * 
 * @param pRoot    the source node of the reference
 * @param pNode    the current node to be referenced
 * @param divSet   the set of nodes should be kept
 * @return int     the number of nodes referenced
 */
int NetMan::NodeRef_rec(AbcObj* pRoot, AbcObj* pNode, unordered_set<AbcObj*>& divSet) const {
    AbcObj* pFanin;
    int i, Counter = 1;
    if ( Abc_ObjIsCi(pNode) || divSet.count(pNode) || (pRoot != pNode && IsPoDriver(pNode)) )
        return 0;
    Abc_ObjForEachFanin( pNode, pFanin, i ) {
        if ( pFanin->vFanouts.nSize++ == 0 )
            Counter += NodeRef_rec(pRoot, pFanin, divSet);
    }
    return Counter;
}


/**
 * @brief Recursively dereference the node pNode and its MFFCs, and collect the nodes to be deleted
 * @brief delNodes is not cleared in this function
 * 
 * 
 * @param rootId    the source node of the dereference
 * @param nodeId    the current node to be dereferenced
 * @param delNodes  the vector of nodes to be deleted
 * @retval delNodes the vector of nodes to be deleted
 * @return int      the number of nodes dereferenced
 */
int NetMan::NodeDeref_rec_v2(int rootId, int nodeId, IntVect& delNodes) const {
    int count = 1;
    // if pNode is a PI, then it should be kept
    // if pNode is in divSet, then it should be kept
    // if pNode is a PO driver, and pNode is not the root (the source node of the dereference), then it should be kept
    if (IsObjPi(nodeId) || (rootId != nodeId && IsPoDriver(nodeId)))
        return 0;
    for (int i = 0; i < GetFaninNum(nodeId); ++i) {
        auto pFanin = GetFanin(nodeId, i);
        assert(pFanin->vFanouts.nSize > 0);
        if (--pFanin->vFanouts.nSize == 0)
            count += NodeDeref_rec_v2(rootId, pFanin->Id, delNodes);
    }
    delNodes.emplace_back(nodeId);
    return count;
}


/**
 * @brief Recursively reference the node pNode and its MFFCs
 * @brief The inverse operation of NodeDeref_rec_v2
 * 
 * @param pRoot    the source node of the reference
 * @param pNode    the current node to be referenced
 * @return int     the number of nodes referenced
 */
int NetMan::NodeRef_rec_v2(int rootId, int nodeId) const {
    int count = 1;
    if (IsObjPi(nodeId) || (rootId != nodeId && IsPoDriver(nodeId)))
        return 0;
    for (int i = 0; i < GetFaninNum(nodeId); ++i) {
        auto pFanin = GetFanin(nodeId, i);
        if (pFanin->vFanouts.nSize++ == 0)
            count += NodeRef_rec_v2(rootId, nodeId);
    }
    return count;
}


int NetMan::GetSizeGain(int rootId, const IntVect& divIds) const {
    // get divisor set
    AbcObjSet divSet; 
    for (int divId: divIds)
        divSet.emplace(GetObj(divId));

    // get size gain
    auto pRoot = GetObj(rootId);
    assert(IsNode(pRoot));
    int nSizeGain = NodeDeref_rec(pRoot, pRoot, divSet);
    NodeRef_rec(pRoot, pRoot, divSet);
    return nSizeGain;
}


int NetMan::GetSizeGain(const IntVect& targetIds, const IntVect& divIds) const {
    // get divisor set
    AbcObjSet divSet; 
    for (int divId: divIds)
        divSet.emplace(GetObj(divId));
    
    // get size gain
    int nSizeGain = 0;
    IntSet skipNodes;
    for (int targetId: targetIds) {
        auto pTarget = GetObj(targetId);
        assert(IsNode(pTarget));
        if (pTarget->vFanouts.nSize == 0) {
            skipNodes.emplace(targetId);
            continue;
        }
        nSizeGain += NodeDeref_rec(pTarget, pTarget, divSet);
    }

    int nSizeGain2 = 0;
    for (auto iter = targetIds.rbegin(); iter != targetIds.rend(); ++iter) {
        auto pTarget = GetObj(*iter);
        assert(IsNode(pTarget));
        if (skipNodes.count(*iter))
            continue;
        nSizeGain2 += NodeRef_rec(pTarget, pTarget, divSet);
    }
    assert(nSizeGain == nSizeGain2);
    
    return nSizeGain;
}




