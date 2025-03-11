/**
 * @file lac.cc
 * @author Chang Meng (chang.meng@epfl.ch)
 * @brief Local approximate change (LAC) management
 * 
 */
#include "lac.h"


using namespace std;
using namespace abc;


/**
 * @brief Generate constant LACs: Shin, Doochul, and Sandeep K. Gupta. "Approximate logic synthesis for error tolerant applications." 2010 Design, Automation & Test in Europe Conference & Exhibition (DATE 2010). IEEE, 2010.
 * 
 * @param net     the network
 * @retval pLacs  generated LACs
 * @return void
 */
void LACMan::GenConstLACs(const NetMan& net) {
    if (net.GetNetType() != NET_TYPE::SOP) {
        fmt::print(stderr, "Error: only support generating LACs on SOP network.\n");
        assert(0);
    }
    pLacs.clear();
    pLacs.reserve(net.GetNodeNum() << 1);
    for (int nodeId = 0; nodeId < net.GetIdMaxPlus1(); ++nodeId) {
        if (net.IsNode(nodeId) && !net.IsConst(nodeId)) {
            int sizeGain = net.GetSizeGain(nodeId, IntVect{});
            pLacs.emplace_back(make_shared<LAC>(nodeId, sizeGain, IntVect{}, string(" 1\n")));
            pLacs.emplace_back(make_shared<LAC>(nodeId, sizeGain, IntVect{}, string(" 0\n")));
        }
    }

    // print
    fmt::print("generated {} constant LACs\n", pLacs.size());
}


/**
 * @brief Generate SASIMI LACs: Venkataramani, Swagath, Kaushik Roy, and Anand Raghunathan. "Substitute-and-simplify: A unified design paradigm for approximate and quality configurable circuits." 2013 Design, Automation & Test in Europe Conference & Exhibition (DATE). IEEE, 2013.
 * 
 * @param net     the network
 * @retval pLacs  generated LACs
 * @return void
 */
void LACMan::GenSasimiLACs(const NetMan& net, int maxCandResub, bool inclConst) {
    // prepare
    if (net.GetNetType() != NET_TYPE::SOP) {
        fmt::print(stderr, "Error: only support generating LACs on SOP network.\n");
        assert(0);
    }
    pLacs.clear();

    // special case: constant replacement
    if (inclConst) {
    for (int nodeId = 0; nodeId < net.GetIdMaxPlus1(); ++nodeId) {
        if (net.IsNode(nodeId) && !net.IsConst(nodeId)) {
            int sizeGain = net.GetSizeGain(nodeId, IntVect{});
            pLacs.emplace_back(make_shared<LAC>(nodeId, sizeGain, IntVect{}, string(" 1\n")));
            pLacs.emplace_back(make_shared<LAC>(nodeId, sizeGain, IntVect{}, string(" 0\n")));
        }
    }
    }

    // general case: substitute the target signal with the substitution signal
    IntVect targIds;
    targIds.reserve(net.GetNodeNum());
    for (int nodeId = 0; nodeId < net.GetIdMaxPlus1(); ++nodeId) {
        if (net.IsNode(nodeId) && !net.IsConst(nodeId) && net.GetFaninNum(nodeId) > 1)
            targIds.emplace_back(nodeId);
    }
    if (targIds.empty()) {
        fmt::print("no target node for SASIMI LACs\n");
        return;
    }
    net.GetLev();
    const int MAX_LAC_PER_NODE = std::max(1, static_cast<int>(maxCandResub / targIds.size()));
    for (int targId: targIds) {
        int lacNum = 0;
        for (int subId = 0; subId < net.GetIdMaxPlus1(); ++subId) {
            if (!net.IsObj(subId) || net.IsObjPo(subId) || net.IsConst(subId) || targId == subId)
                continue;
            if (net.GetObjLev(subId) < net.GetObjLev(targId)) {
                int sizeGain = net.GetSizeGain(targId, IntVect{subId});
                pLacs.emplace_back(make_shared<LAC>(targId, sizeGain, IntVect{subId}, string("1 1\n")));
                pLacs.emplace_back(make_shared<LAC>(targId, sizeGain, IntVect{subId}, string("0 1\n")));
                lacNum += 2;
                if (lacNum >= MAX_LAC_PER_NODE)
                    break;
            }
        }
    }

    // print
    fmt::print("generated {} SASIMI LACs\n", pLacs.size());
}


/**
 * @brief Generate resubstitution-based LACs: Chang Meng, Alan Mishchenko, Weikang Qian, and Giovanni De Micheli, "Efficient Resubstitution-Based Approximate Logic Synthesis," in IEEE Transactions on Computer-Aided Design of Integrated Circuits and Systems (TCAD), 2024.
 * 
 * @param net              the network
 * @param seed             random seed
 * @param nFrame4ResubGen  number of simulation frames for approximate resubstitution generation
 * @param maxLevelDiff     maximum level difference when generating approximate resubstitutions
 * @param maxCandResub     maximum number of candidate approximate resubstitutions
 * @retval pLacs           generated LACs
 * @return void
 */
void LACMan::GenResubLACs(const NetMan& net, unsigned seed, int nFrame4ResubGen, int maxLevelDiff, int maxCandResub, bool inclConst) {
    fmt::print("generating resubstitution-based LACs\n");
    // initialize
    int simulationFrame = nFrame4ResubGen;
    int halfFrame = simulationFrame >> 1;
    int LAC_NUM_LIMIT = maxCandResub;
    assert(net.GetNetType() == NET_TYPE::SOP);
    pLacs.clear();
    net.GetLev();
    Abc_NtkStartReverseLevels(net.GetNet(), 0);

    // simulate
    Simulator smlt(net, seed, simulationFrame, DISTR_TYPE::UNIF);
    smlt.LogicSim();

    // collect target nodes
    IntVect targIds;
    targIds.reserve(net.GetNodeNum());
    for (int nodeId = 0; nodeId < net.GetIdMaxPlus1(); ++nodeId) {
        if (net.IsNode(nodeId) && !net.IsConst(nodeId) && net.GetFaninNum(nodeId) > 1)
            targIds.emplace_back(nodeId);
    }

    // collect divisors
    vector<IntVect> divs4Nodes;
    divs4Nodes.resize(net.GetIdMaxPlus1());
    boost::timer::progress_display pd(targIds.size());
    for (int targId: targIds) {
        auto pTarget = net.GetObj(targId);
        GetDivs(pTarget, abc::Abc_ObjRequiredLevel(pTarget) - 1, divs4Nodes[targId]);
        ++pd;
    }
    // generate 0 resubstitution
    if (inclConst) {
    for (int targId: targIds) {
        int sizeGain = net.IsTheOnlyPoDriver(targId)? net.GetSizeGain(targId, IntVect{}): net.GetSizeGain(targId, IntVect{}) + 1;
        if (static_cast<int>(smlt.GetDat(targId).count()) <= halfFrame)
            pLacs.emplace_back(make_shared<LAC>(targId, sizeGain, IntVect{}, string(" 0\n")));
        else
            pLacs.emplace_back(make_shared<LAC>(targId, sizeGain, IntVect{}, string(" 1\n")));
    }
    fmt::print("generated 0-resubs (constant LACs), total #lacs = {}\n", pLacs.size());
    }
    // generate 1 resubstitution
    for (int targId: targIds) {
        auto& divs = divs4Nodes[targId];
        for (int div: divs) {
            auto diff = (smlt.GetDat(div) ^ smlt.GetDat(targId)).count();
            if (diff == 0) {
                int sizeGain = net.GetSizeGain(targId, IntVect{div});
                pLacs.emplace_back(make_shared<LAC>(targId, sizeGain, IntVect{div}, string("1 1\n")));
                if (static_cast<int>(pLacs.size()) > LAC_NUM_LIMIT)
                    break;
            }
            else if (static_cast<int>(diff) == simulationFrame) {
                int sizeGain = net.GetSizeGain(targId, IntVect{div});
                pLacs.emplace_back(make_shared<LAC>(targId, sizeGain, IntVect{div}, string("0 1\n")));
                if (static_cast<int>(pLacs.size()) > LAC_NUM_LIMIT)
                    break;
            }
        }
    }
    fmt::print("generated 1-resubs, total #lacs = {}\n", pLacs.size());
    // generate 2 resubstitution: try replacing the i-th fanin with another divisor
    if (static_cast<int>(pLacs.size()) <= LAC_NUM_LIMIT) {
    bool _break = false;
    for (int targId: targIds) {
        if (_break)
            break;
        int nFanin = net.GetFaninNum(targId);
        assert(nFanin == 2);
        int fanin0 = net.GetFaninId(targId, 0), fanin1 = net.GetFaninId(targId, 1);
        auto& divs = divs4Nodes[targId];
        for (int i = 0; i < nFanin; ++i) {
            int remainedFanin = (i == 0)? fanin1: fanin0;
            int replacedFanin = net.GetFaninId(targId, i);
            auto faninIds = IntVect{remainedFanin, -1};
            for (int div: divs) {
                if (div == replacedFanin || div == remainedFanin)
                    continue;
                if (_break)
                    break;
                faninIds[1] = div;
                int sizeGain = net.GetSizeGain(targId, faninIds) - 1;
                if (sizeGain >= 1) {
                    for (int comb = 0; comb < 4; ++comb) {
                        int var0 = (comb >> 1) & 1, var1 = comb & 1;
                        auto dat0 = var0? smlt.GetDat(faninIds[0]): ~(smlt.GetDat(faninIds[0]));
                        auto dat1 = var1? smlt.GetDat(faninIds[1]): ~(smlt.GetDat(faninIds[1]));
                        auto res = dat0 & dat1;
                        auto diff = (res ^ smlt.GetDat(targId)).count();
                        if (diff == 0) {
                            pLacs.emplace_back(make_shared<LAC>(targId, sizeGain, faninIds, to_string(var0) + to_string(var1) + string(" 1\n")));
                            _break = (static_cast<int>(pLacs.size()) > LAC_NUM_LIMIT);
                        }
                        else if (static_cast<int>(diff) == simulationFrame) {
                            pLacs.emplace_back(make_shared<LAC>(targId, sizeGain, faninIds, to_string(var0) + to_string(var1) + string(" 0\n")));
                            _break = (static_cast<int>(pLacs.size()) > LAC_NUM_LIMIT);
                        }
                    }
                }
            }
        }
    }
    fmt::print("generated 2-resubs, total #lacs = {}\n", pLacs.size());
    }
    // for efficiency issues, do not generate LACs that replaces both fanins with new divisors
    // generate 2 resubstitution: try replacing both fanins with new divisors, using AND-based functions
    // if (static_cast<int>(pLacs.size()) <= LAC_NUM_LIMIT) {
    // bool _break = false;
    // for (int targId: targIds) {
    //     if (_break)
    //         break;
    //     assert(net.GetFaninNum(targId) == 2);
    //     int fanin0 = net.GetFaninId(targId, 0), fanin1 = net.GetFaninId(targId, 1);
    //     auto& divs = divs4Nodes[targId];
    //     IntVect posDiv4And, negDiv4And, posDiv4Or, negDiv4Or;
    //     posDiv4And.reserve(divs.size());
    //     negDiv4And.reserve(divs.size());
    //     posDiv4Or.reserve(divs.size());
    //     negDiv4Or.reserve(divs.size());
    //     for (int div: divs) {
    //         if (div == fanin0 || div == fanin1)
    //             continue;
    //         int levelDiff = abc::Abc_ObjReverseLevel(net.GetObj(div)) - abc::Abc_ObjReverseLevel(net.GetObj(targId));
    //         if (levelDiff > maxLevelDiff)
    //             continue;
    //         auto check = (smlt.GetDat(targId) & smlt.GetDat(div)) | (~smlt.GetDat(targId)); // targ = div & xxx, targ = 1 \Rightarrow div = 1
    //         if (check.all())
    //             posDiv4And.emplace_back(div);
    //         check = (smlt.GetDat(targId) & ~smlt.GetDat(div)) | (~smlt.GetDat(targId)); // targ = ~div & xxx, targ = 1 \Rightarrow div = 0
    //         if (check.all())
    //             negDiv4And.emplace_back(div);
    //         check = (~smlt.GetDat(targId) & ~smlt.GetDat(div)) | (smlt.GetDat(targId)); // targ = div | xxx, targ = 0 \Rightarrow div = 0
    //         if (check.all())
    //             posDiv4Or.emplace_back(div);
    //         check = (~smlt.GetDat(targId) & smlt.GetDat(div)) | (smlt.GetDat(targId)); // targ = ~div | xxx, targ = 0 \Rightarrow div = 1
    //         if (check.all())
    //             negDiv4Or.emplace_back(div);
    //     }
    //     // try targ = div0 & xxx
    //     for (int div0: posDiv4And) {
    //         // try targ = div0 & div1
    //         for (int div1: posDiv4And) {
    //             int sizeGain = net.GetSizeGain(targId, IntVect{div0, div1}) - 1;
    //             if (sizeGain >= 1) {
    //                 auto diff = (smlt.GetDat(div0) & smlt.GetDat(div1)) ^ smlt.GetDat(targId);
    //                 if (diff.none()) {
    //                     pLacs.emplace_back(make_shared<LAC>(targId, sizeGain, IntVect{div0, div1}, string("11 1\n")));
    //                     _break = (static_cast<int>(pLacs.size()) > LAC_NUM_LIMIT);
    //                 }
    //             }
    //         }
    //         // try targ = div0 & ~div1
    //         for (int div1: negDiv4And) {
    //             int sizeGain = net.GetSizeGain(targId, IntVect{div0, div1}) - 1;
    //             if (sizeGain >= 1) {
    //                 auto diff = (smlt.GetDat(div0) & ~smlt.GetDat(div1)) ^ smlt.GetDat(targId);
    //                 if (diff.none()) {
    //                     pLacs.emplace_back(make_shared<LAC>(targId, sizeGain, IntVect{div0, div1}, string("10 1\n")));
    //                     _break = (static_cast<int>(pLacs.size()) > LAC_NUM_LIMIT);
    //                 }
    //             }
    //         }
    //     }
    //     // try targ = ~div0 & xxx
    //     for (int div0: negDiv4And) {
    //         // try targ = ~div0 & div1
    //         for (int div1: posDiv4And) {
    //             int sizeGain = net.GetSizeGain(targId, IntVect{div0, div1}) - 1;
    //             if (sizeGain >= 1) {
    //                 auto diff = (~smlt.GetDat(div0) & smlt.GetDat(div1)) ^ smlt.GetDat(targId);
    //                 if (diff.none()) {
    //                     pLacs.emplace_back(make_shared<LAC>(targId, sizeGain, IntVect{div0, div1}, string("01 1\n")));
    //                     _break = (static_cast<int>(pLacs.size()) > LAC_NUM_LIMIT);
    //                 }
    //             }
    //         }
    //         // try targ = ~div0 & ~div1
    //         for (int div1: negDiv4And) {
    //             int sizeGain = net.GetSizeGain(targId, IntVect{div0, div1}) - 1;
    //             if (sizeGain >= 1) {
    //                 auto diff = (~smlt.GetDat(div0) & ~smlt.GetDat(div1)) ^ smlt.GetDat(targId);
    //                 if (diff.none()) {
    //                     pLacs.emplace_back(make_shared<LAC>(targId, sizeGain, IntVect{div0, div1}, string("00 1\n")));
    //                     _break = (static_cast<int>(pLacs.size()) > LAC_NUM_LIMIT);
    //                 }
    //             }
    //         }
    //     }
    //     // try targ = div0 | xxx
    //     for (int div0: posDiv4Or) {
    //         // try targ = div0 | div1
    //         for (int div1: posDiv4Or) {
    //             int sizeGain = net.GetSizeGain(targId, IntVect{div0, div1}) - 1;
    //             if (sizeGain >= 1) {
    //                 auto diff = (smlt.GetDat(div0) | smlt.GetDat(div1)) ^ smlt.GetDat(targId);
    //                 if (diff.none()) {
    //                     pLacs.emplace_back(make_shared<LAC>(targId, sizeGain, IntVect{div0, div1}, string("00 0\n")));
    //                     _break = (static_cast<int>(pLacs.size()) > LAC_NUM_LIMIT);
    //                 }
    //             }
    //         }
    //         // try targ = div0 | ~div1
    //         for (int div1: negDiv4Or) {
    //             int sizeGain = net.GetSizeGain(targId, IntVect{div0, div1}) - 1;
    //             if (sizeGain >= 1) {
    //                 auto diff = (smlt.GetDat(div0) | ~smlt.GetDat(div1)) ^ smlt.GetDat(targId);
    //                 if (diff.none()) {
    //                     pLacs.emplace_back(make_shared<LAC>(targId, sizeGain, IntVect{div0, div1}, string("01 0\n")));
    //                     _break = (static_cast<int>(pLacs.size()) > LAC_NUM_LIMIT);
    //                 }
    //             }
    //         }
    //     }
    //     // try targ = ~div0 | xxx
    //     for (int div0: negDiv4Or) {
    //         // try targ = ~div0 | div1
    //         for (int div1: posDiv4Or) {
    //             int sizeGain = net.GetSizeGain(targId, IntVect{div0, div1}) - 1;
    //             if (sizeGain >= 1) {
    //                 auto diff = (~smlt.GetDat(div0) | smlt.GetDat(div1)) ^ smlt.GetDat(targId);
    //                 if (diff.none()) {
    //                     pLacs.emplace_back(make_shared<LAC>(targId, sizeGain, IntVect{div0, div1}, string("10 0\n")));
    //                     _break = (static_cast<int>(pLacs.size()) > LAC_NUM_LIMIT);
    //                 }
    //             }
    //         }
    //         // try targ = ~div0 | ~div1
    //         for (int div1: negDiv4Or) {
    //             int sizeGain = net.GetSizeGain(targId, IntVect{div0, div1}) - 1;
    //             if (sizeGain >= 1) {
    //                 auto diff = (~smlt.GetDat(div0) | ~smlt.GetDat(div1)) ^ smlt.GetDat(targId);
    //                 if (diff.none()) {
    //                     pLacs.emplace_back(make_shared<LAC>(targId, sizeGain, IntVect{div0, div1}, string("11 0\n")));
    //                     _break = (static_cast<int>(pLacs.size()) > LAC_NUM_LIMIT);
    //                 }
    //             }
    //         }
    //     }
    // }
    // fmt::print("generated more 2-resubs, total #lacs = {}\n", pLacs.size());
    // }

    // clean up
    abc::Abc_NtkStopReverseLevels(net.GetNet());

}


/**
 * @brief Re-group LACs by node if not done
 * 
 * @retval node2Lacs    node2Lacs[nodeId] = LACs of the node
 * @return void
 */
void LACMan::RegroupLACsByNode(bool forceUpd) {
    // fmt::print("Regrouping LACs by nodes\n");
    if (!forceUpd && !node2Lacs.empty()) {
        fmt::print("LACs have been grouped by nodes\n");
        return;
    }
    node2Lacs.clear();
    for (const auto& pLac: pLacs) {
        int targId = pLac->GetTargId();
        if (node2Lacs.count(targId) == 0)
            node2Lacs[targId] = LACPtrVect{};
        node2Lacs[targId].emplace_back(pLac);
    }
}


extern "C" {
    Vec_Ptr_t * Abc_MfsWinMarkTfi(Abc_Obj_t * pNode);
    void Abc_MfsWinSweepLeafTfo_rec(Abc_Obj_t * pObj, int nLevelLimit);
}
/**
 * @brief Get divisors for a node
 * 
 * @param pNode       the target node
 * @param nLevDivMax  maximum level difference
 * @param divs        the collected divisors
 * @retval divs       the collected divisors
 * @return void
 */
void LACMan::GetDivs(Abc_Obj_t * pNode, int nLevDivMax, IntVect& divs) {
    const int nWinMax = 300;
    const int nFanoutsMax = 30;
    divs.clear();
    Vec_Ptr_t * vCone, * vDivs;
    Abc_Obj_t * pObj, * pFanout, * pFanin;
    int k, f, m;
    int nDivsPlus = 0, nTrueSupp;

    // mark the TFI with the current trav ID
    Abc_NtkIncrementTravId( pNode->pNtk );
    vCone = Abc_MfsWinMarkTfi( pNode );

    // count the number of PIs
    nTrueSupp = 0;
    Vec_PtrForEachEntry( Abc_Obj_t *, vCone, pObj, k )
        nTrueSupp += Abc_ObjIsCi(pObj);
//    printf( "%d(%d) ", Vec_PtrSize(p->vSupp), m );

    // mark with the current trav ID those nodes that should not be divisors:
    // (1) the node and its TFO
    // (2) the MFFC of the node
    // (3) the node's fanins (these are treated as a special case)
    Abc_NtkIncrementTravId( pNode->pNtk );
    Abc_MfsWinSweepLeafTfo_rec( pNode, nLevDivMax );
//    Abc_MfsWinVisitMffc( pNode );
    Abc_ObjForEachFanin( pNode, pObj, k )
        Abc_NodeSetTravIdCurrent( pObj );

    // at this point the nodes are marked with two trav IDs:
    // nodes to be collected as divisors are marked with previous trav ID
    // nodes to be avoided as divisors are marked with current trav ID

    // start collecting the divisors
    vDivs = Vec_PtrAlloc( nWinMax );
    Vec_PtrForEachEntry( Abc_Obj_t *, vCone, pObj, k )
    {
        if ( !Abc_NodeIsTravIdPrevious(pObj) )
            continue;
        if ( (int)pObj->Level > nLevDivMax )
            continue;
        Vec_PtrPush( vDivs, pObj );
        if ( Vec_PtrSize(vDivs) >= nWinMax )
            break;
    }
    Vec_PtrFree( vCone );

    // explore the fanouts of already collected divisors
    if ( Vec_PtrSize(vDivs) < nWinMax )
    Vec_PtrForEachEntry( Abc_Obj_t *, vDivs, pObj, k )
    {
        // consider fanouts of this node
        Abc_ObjForEachFanout( pObj, pFanout, f )
        {
            // stop if there are too many fanouts
            if ( nFanoutsMax && f > nFanoutsMax )
                break;
            // skip nodes that are already added
            if ( Abc_NodeIsTravIdPrevious(pFanout) )
                continue;
            // skip nodes in the TFO or in the MFFC of node
            if ( Abc_NodeIsTravIdCurrent(pFanout) )
                continue;
            // skip COs
            if ( !Abc_ObjIsNode(pFanout) )
                continue;
            // skip nodes with large level
            if ( (int)pFanout->Level > nLevDivMax )
                continue;
            // skip nodes whose fanins are not divisors  -- here we skip more than we need to skip!!! (revise later)  August 7, 2009
            Abc_ObjForEachFanin( pFanout, pFanin, m )
                if ( !Abc_NodeIsTravIdPrevious(pFanin) )
                    break;
            if ( m < Abc_ObjFaninNum(pFanout) )
                continue;
            // make sure this divisor in not among the nodes
//            Vec_PtrForEachEntry( Abc_Obj_t *, p->vNodes, pFanin, m )
//                assert( pFanout != pFanin );
            // add the node to the divisors
            Vec_PtrPush( vDivs, pFanout );
            // Vec_PtrPushUnique( p->vNodes, pFanout );
            Abc_NodeSetTravIdPrevious( pFanout );
            nDivsPlus++;
            if ( Vec_PtrSize(vDivs) >= nWinMax )
                break;
        }
        if ( Vec_PtrSize(vDivs) >= nWinMax )
            break;
    }

    // sort the divisors by level in the increasing order
    Vec_PtrSort( vDivs, (int (*)(const void *, const void *))Abc_NodeCompareLevelsIncrease );

    // add the fanins of the node
    Abc_ObjForEachFanin( pNode, pFanin, k )
        Vec_PtrPush( vDivs, pFanin );
    
    divs.reserve(Vec_PtrSize(vDivs));
    Vec_PtrForEachEntry(Abc_Obj_t *, vDivs, pObj, k)
        divs.emplace_back(pObj->Id);

    // clean up
    Vec_PtrFree(vDivs);
}


/**
 * @brief Print LACs
 * 
 * @return void
 */
void LACMan::PrintLACs(int firstK) const {
    if (firstK == -1 || firstK > static_cast<int>(pLacs.size()))
        firstK = pLacs.size();
    fmt::print("{}first {} LACs{}\n", HALF_DASH_LINE, firstK, HALF_DASH_LINE);
    for (int i = 0; i < firstK; ++i)
        fmt::print("{}\n", pLacs[i]->ToStr());
    fmt::print("{}\n", DASH_LINE);
}


/**
 * @brief Get the best LAC
 * 
 * @return the best LAC; nullptr if no LAC is found
 */
shared_ptr<LAC> LACMan::GetBestLac(double errUppBound) const {
    // special case: no LAC
    if (pLacs.empty())
        return nullptr;
    // special case: only one LAC
    if (pLacs.size() == 1)
        return pLacs[0];
    // general case
    shared_ptr<LAC> pBestLac = pLacs[0];
    for (int i = 1; i < static_cast<int>(pLacs.size()); ++i) {
        auto pLac = pLacs[i];
        double err = pLac->GetErr();
        if (DoubleGreat(err, errUppBound))
            continue;
        if (Lac0BetterThanLac1(pLac, pBestLac))
            pBestLac = pLac;
    }
    return pBestLac;
}


/**
 * @brief Sort and keep the top k LACs
 * 
 * @param  k      the number of LACs to keep; if k = -1, keep all LACs
 * @retval pLacs  the top k LACs
 * @return void
 */
void LACMan::SortAndKeepTopKLACs(int k) {
    if (k < -1) {
        fmt::print(stderr, "Error: k should be larger than or equal to -1.\n");
        assert(0);
        return;
    }
    // sort the LACs
    // std::sort(pLacs.begin(), pLacs.end(), LacBetterThanUsingErr2);
    std::sort(pLacs.begin(), pLacs.end(), Lac0BetterThanLac1);
    // if k is larger than the number of LACs, keep all LACs
    if (k >= static_cast<int>(pLacs.size()))
        return;
    // keep the top k LACs
    if (k != -1)
        pLacs.resize(k);
}


/**
 * @brief Remove LACs with the given representative strings
 * 
 * @param blackList     a black list of LACs
 * @retval pLacs        the LACs after removal
 * @return void
 */
void LACMan::RemLacsFromBlackList(std::unordered_set<std::string>& blackList) {
    if (blackList.empty())
        return;
    pLacs.erase(
        std::remove_if(pLacs.begin(), pLacs.end(),
            [&blackList](const LACPtr& pLac) {
                return blackList.count(pLac->ToStrShort());
            }
        ),
        pLacs.end()
    );
}

/**
 * @brief Remove LACs with error larger than the given upper bound
 * 
 * @param errUppBound     the upper bound of the error
 * @retval pLacs           the LACs after removal
 * @return void
 */
void LACMan::RemLargeErrLacs(double errUppBound) {
    pLacs.erase(
        std::remove_if(pLacs.begin(), pLacs.end(),
            [errUppBound](const LACPtr& pLac) {
                return DoubleGreat(pLac->GetErr(), errUppBound);
            }
        ),
        pLacs.end()
    );
}


/**
 * @brief Apply the LAC to the network net
 * 
 * @param net   the network
 * @param lac   the LAC to be applied
 * @retval net  the network after applying the LAC
 * @return void
 */
void ApplyLac(NetMan& net, const LAC& lac) {
    // prepare
    assert(net.GetNetType() == NET_TYPE::SOP);
    net.GetLev();
    auto targId = lac.GetTargId();
    auto faninIds = lac.GetDivIds();
    auto sop = lac.GetSop();
    auto _sop = sop;
    replace(_sop.begin(), _sop.end(), '\n', ';');
    fmt::print("replace {}(l={}) with old fanins [", *net.GetObj(targId), net.GetObjLev(targId));
    for (ll i = 0; i < net.GetFaninNum(targId); ++i)
        fmt::print("{}(l={}),", *net.GetFanin(targId, i), net.GetObjLev(net.GetFanin(targId, i)));
    fmt::print("] by divisors [");
    for (const auto & faninId: faninIds)
        fmt::print("{}(l={}),", *net.GetObj(faninId), net.GetObjLev(faninId));
    fmt::print("] using function [{}], ", _sop);
    fmt::print("max error = {}, ", lac.GetErr());
    fmt::print("estimated size gain = {}\n", lac.GetSizeGain());

    // perform replacement
    auto consts = net.CreateConstsIfNotExist();
    if (sop == " 0\n") {
        net.Replace(targId, consts.first);
        net.PropConst(consts.first, false, false);
    }
    else if (sop == " 1\n") {
        net.Replace(targId, consts.second);
        net.PropConst(consts.second, false, false);
    }
    else if (sop == "1 1\n") {
        assert(faninIds.size() == 1);
        net.Replace(targId, faninIds[0]);
    }
    else if (sop == "0 1\n") {
        assert(faninIds.size() == 1);
        net.ReplaceByComplementedObj(targId, faninIds[0]);
    }
    else {
        int newNodeId = net.CreateAIGStyleNodes(faninIds, sop);
        net.Replace(targId, newNodeId);
    }

    // remove redundant nodes
    net.CleanUp(true);
}


/**
 * @brief Temporarily apply the LACs to the network net
 * 
 * @param net         the network
 * @param lac         the LAC to be applied
 * @param replTrace   the replacement trace used to recover the network; {pTS, pSS, fanout0, iFanin0, fanout1, iFanin1, ..., -1, delObj0, delObj1, ...}; pTS is the iFanin-k-th fanin of pTS's fanout-k; delObj is the object to be deleted
 * @param fVerb       whether to print verbose information
 * @retval replTrace  the replacement trace used to recover the network
 * @return ssId       the substitution signal ID
 */
int TempApplyLac(NetMan& net, const LAC& lac, IntVect& replTrace, bool fVerb) {
    assert(net.GetNetType() == NET_TYPE::SOP);
    int ssId = -1;
    if (lac.IsConst()) {
        auto constIds = net.CreateConstsIfNotExist();
        ssId = lac.IsConst0()? constIds.first: constIds.second;
        net.TempRepl_v2(lac.GetTargId(), ssId, replTrace, fVerb);
    }
    else {
        ssId = net.CreateNode(lac.GetDivIds(), lac.GetSop());
        net.TempRepl_v2(lac.GetTargId(), ssId, replTrace, fVerb);
        replTrace.emplace_back(-1);
        replTrace.emplace_back(ssId);
    }
    assert(ssId != -1);
    return ssId;
}


/**
 * @brief Temporarily apply the LAC to the network net
 * @brief Use a controlling signal to control whether to apply the LAC or not
 * 
 * @param net         the network
 * @param lac         the LAC to be applied
 * @param ctrlId      the controlling signal ID
 * @param replTrace   the replacement trace used to recover the network; {pTS, pSS, fanout0, iFanin0, fanout1, iFanin1, ..., -1, delObj0, delObj1, ...}; pTS is the iFanin-k-th fanin of pTS's fanout-k; delObj is the object to be deleted
 * @param fVerb      whether to print verbose information
 * @retval replTrace  the replacement trace used to recover the network
 * @return void
 */
void TempApplyLacWithController(NetMan& net, const LAC& lac, int ctrlId, IntVect& replTrace, bool fVerb) {
    // fmt::print("Temprarily apply LAC: {}\n", lac);
    // create a new node for substitution
    int tsId = lac.GetTargId();
    int ssId = net.CreateNode(lac.GetDivIds(), lac.GetSop());
    // create a controlling signal
    auto pControl = net.CreatePi(fmt::format("ctrl_{}", ctrlId).c_str());
    // create a MUX
    // fmt::print("{}, {}, {}\n", *net.GetObj(tsId), *net.GetObj(ssId), *pControl);
    int muxId = net.CreateNode({tsId, ssId, pControl->Id}, "1-0 1\n-11 1\n");
    // replace the target 
    // fmt::print("{}, {}\n", *pTS, *net.GetObj(muxId));
    net.TempRepl_v2(tsId, muxId, replTrace, fVerb);
    // add nodes to be deleted
    replTrace.emplace_back(-1);
    replTrace.emplace_back(muxId);
    replTrace.emplace_back(pControl->Id);
    replTrace.emplace_back(ssId);
}


/**
 * @brief Recover the network net from the replacement traces
 * 
 * @param net         the network
 * @param replTraces  the replacement traces
 * @param fVerb      whether to print verbose information
 */
void RecovNet(NetMan& net, const Int2DVect& replTraces, bool fVerb) {
    for (auto it = replTraces.rbegin(); it != replTraces.rend(); ++it)
        net.Recov_v2(*it, fVerb);
}