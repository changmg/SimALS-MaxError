#include <gtest/gtest.h>
#include <assert.h>
#include <vector>
#include <cryptominisat.h>
#include "error.h"


using std::vector;
using namespace CMSat;


/**
 * @brief Demonstrate some basic assertions.
 * 
 */
TEST(HelloTest, BasicAssertions) {
    // Expect two strings not to be equal.
    EXPECT_STRNE("hello", "world");
    // Expect equality.
    EXPECT_EQ(7 * 6, 42);
}


/**
 * @brief Official example of cryptominisat
 * 
 */
TEST(CMSatTest, OfficialExample) {
    SATSolver solver;
    vector<Lit> clause;

    //Let's use 4 threads
    solver.set_num_threads(4);

    //We need 3 variables. They will be: 0,1,2
    //Variable numbers are always trivially increasing
    solver.new_vars(3);

    //add "1 0"
    clause.push_back(Lit(0, false));
    solver.add_clause(clause);

    //add "-2 0"
    clause.clear();
    clause.push_back(Lit(1, true));
    solver.add_clause(clause);

    //add "-1 2 3 0"
    clause.clear();
    clause.push_back(Lit(0, true));
    clause.push_back(Lit(1, false));
    clause.push_back(Lit(2, false));
    solver.add_clause(clause);

    //solve the problem
    lbool ret = solver.solve();
    EXPECT_EQ(ret, l_True);
    // std::cout
    // << "Solution is: "
    // << solver.get_model()[0]
    // << ", " << solver.get_model()[1]
    // << ", " << solver.get_model()[2]
    // << std::endl;

    //assumes 3 = FALSE, no solutions left
    vector<Lit> assumptions;
    assumptions.push_back(Lit(2, true));
    ret = solver.solve(&assumptions);
    EXPECT_EQ(ret, l_False);

    //without assumptions we still have a solution
    ret = solver.solve();
    EXPECT_EQ(ret, l_True);

    //add "-3 0"
    //No solutions left, UNSATISFIABLE returned
    clause.clear();
    clause.push_back(Lit(2, true));
    solver.add_clause(clause);
    ret = solver.solve();
    EXPECT_EQ(ret, l_False);
}


/**
 * @brief Test logic simulator using a multiplier
 * 
 */
TEST(ALSTest, LogicSimulator) {
    GlobStartAbc();

    NetMan net("./als/tests/benchmarks/am8.blif");
    int nFrame = 1 << 16;
    Simulator smlt(net, 0, nFrame, DISTR_TYPE::ENUM);
    smlt.LogicSim();

    for (int i = 0; i < nFrame; ++i) {
        int res = smlt.GetOutputFast(i);
        int op0 = static_cast<int>(smlt.GetInput(i, 0, 7));
        int op1 = static_cast<int>(smlt.GetInput(i, 8, 15));
        if (res != op0 * op1) {
            fmt::print("Error: wrong simulation result\n");
            assert(0);
        }
    }

    GlobStopAbc();
}


/**
 * @brief Test error manager for measuring maximum error
 * 
 */
TEST(ALSTest, ErrManTestC1355) {
    GlobStartAbc();
    NetMan net0("./als/tests/benchmarks/c1355.blif");

    {NetMan net1("./als/tests/benchmarks/c1355_r5_MAXHDxxx_s365_d16.blif");
    ErrMan errMan(net0, net1);
    auto maxErr = errMan.ComputeMaxErr(METR_TYPE::MAXHD);
    fmt::print("maxErr = {}\n", maxErr.str());
    EXPECT_EQ(maxErr, 3);}

    {NetMan net1("./als/tests/benchmarks/c1355_final.blif");
    ErrMan errMan(net0, net1);
    auto maxErr = errMan.ComputeMaxErr(METR_TYPE::MAXHD);
    fmt::print("maxErr = {}\n", maxErr.str());
    EXPECT_EQ(maxErr, 3);}

    GlobStopAbc();
}


/**
 * @brief Test error manager for measuring maximum error
 * 
 */
TEST(ALSTest, ErrManTestMac) {
    GlobStartAbc();
    NetMan net0("./als/tests/benchmarks/mac.aig");
    net0.Comm("logic; sop;");

    {NetMan net1("./als/tests/benchmarks/mac_wce0.blif");
    ErrMan errMan(net0, net1);
    ll maxErr = static_cast<ll>(errMan.ComputeMaxErr(METR_TYPE::MAXED));
    ll maxErrEnum = errMan.GetMaxErrDistUsingEnum();
    EXPECT_EQ(maxErr, maxErrEnum);}

    {NetMan net1("./als/tests/benchmarks/mac_wce8.blif");
    ErrMan errMan(net0, net1);
    ll maxErr = static_cast<ll>(errMan.ComputeMaxErr(METR_TYPE::MAXED));
    ll maxErrEnum = errMan.GetMaxErrDistUsingEnum();
    EXPECT_EQ(maxErr, maxErrEnum);}

    GlobStopAbc();
}


/**
 * @brief Test error manager for measuring maximum error
 * 
 */
TEST(ALSTest, ErrManTestAbsdiff) {
    GlobStartAbc();
    NetMan net0("./als/tests/benchmarks/absdiff.blif");

    {NetMan net1("./als/tests/benchmarks/absdiff_r1_MAXED0_s111_d14.blif");
    ErrMan errMan(net0, net1);
    ll maxErr = static_cast<ll>(errMan.ComputeMaxErr(METR_TYPE::MAXED));
    ll maxErrEnum = errMan.GetMaxErrDistUsingEnum();
    EXPECT_EQ(maxErr, maxErrEnum);}

    {NetMan net1("./als/tests/benchmarks/absdiff_r5_MAXED1_s88_d13.blif");
    ErrMan errMan(net0, net1);
    ll maxErr = static_cast<ll>(errMan.ComputeMaxErr(METR_TYPE::MAXED));
    ll maxErrEnum = errMan.GetMaxErrDistUsingEnum();
    EXPECT_EQ(maxErr, maxErrEnum);}

    {NetMan net1("./als/tests/benchmarks/absdiff_r6_MAXED3_s82_d13.blif");
    ErrMan errMan(net0, net1);
    ll maxErr = static_cast<ll>(errMan.ComputeMaxErr(METR_TYPE::MAXED));
    ll maxErrEnum = errMan.GetMaxErrDistUsingEnum();
    EXPECT_EQ(maxErr, maxErrEnum);}

    {NetMan net1("./als/tests/benchmarks/absdiff_r9_MAXED7_s65_d11.blif");
    ErrMan errMan(net0, net1);
    ll maxErr = static_cast<ll>(errMan.ComputeMaxErr(METR_TYPE::MAXED));
    ll maxErrEnum = errMan.GetMaxErrDistUsingEnum();
    EXPECT_EQ(maxErr, maxErrEnum);}

    {NetMan net1("./als/tests/benchmarks/absdiff_r11_MAXED15_s47_d11.blif");
    ErrMan errMan(net0, net1);
    ll maxErr = static_cast<ll>(errMan.ComputeMaxErr(METR_TYPE::MAXED));
    ll maxErrEnum = errMan.GetMaxErrDistUsingEnum();
    EXPECT_EQ(maxErr, maxErrEnum);}

    {NetMan net1("./als/tests/benchmarks/absdiff_r13_MAXED31_s33_d9.blif");
    ErrMan errMan(net0, net1);
    ll maxErr = static_cast<ll>(errMan.ComputeMaxErr(METR_TYPE::MAXED));
    ll maxErrEnum = errMan.GetMaxErrDistUsingEnum();
    EXPECT_EQ(maxErr, maxErrEnum);}

    {NetMan net1("./als/tests/benchmarks/absdiff_r16_MAXED63_s8_d3.blif");
    ErrMan errMan(net0, net1);
    ll maxErr = static_cast<ll>(errMan.ComputeMaxErr(METR_TYPE::MAXED));
    ll maxErrEnum = errMan.GetMaxErrDistUsingEnum();
    EXPECT_EQ(maxErr, maxErrEnum);}

    GlobStopAbc();
}


/**
 * @brief Test error manager for measuring maximum error
 * 
 */
TEST(ALSTest, ErrManTestAm8) {
    GlobStartAbc();
    NetMan net0("./als/tests/benchmarks/am8.blif");

    {NetMan net1("./als/tests/benchmarks/am8_r4_MAXED0_s467_d35.blif");
    ErrMan errMan(net0, net1);
    ll maxErr = static_cast<ll>(errMan.ComputeMaxErr(METR_TYPE::MAXED));
    ll maxErrEnum = errMan.GetMaxErrDistUsingEnum();
    EXPECT_EQ(maxErr, maxErrEnum);}

    {NetMan net1("./als/tests/benchmarks/am8_r46_MAXED1_s551_d38.blif");
    ErrMan errMan(net0, net1);
    ll maxErr = static_cast<ll>(errMan.ComputeMaxErr(METR_TYPE::MAXED));
    ll maxErrEnum = errMan.GetMaxErrDistUsingEnum();
    EXPECT_EQ(maxErr, maxErrEnum);}

    {NetMan net1("./als/tests/benchmarks/am8_r49_MAXED3_s540_d38.blif");
    ErrMan errMan(net0, net1);
    ll maxErr = static_cast<ll>(errMan.ComputeMaxErr(METR_TYPE::MAXED));
    ll maxErrEnum = errMan.GetMaxErrDistUsingEnum();
    EXPECT_EQ(maxErr, maxErrEnum);}

    {NetMan net1("./als/tests/benchmarks/am8_r53_MAXED5_s534_d38.blif");
    ErrMan errMan(net0, net1);
    ll maxErr = static_cast<ll>(errMan.ComputeMaxErr(METR_TYPE::MAXED));
    ll maxErrEnum = errMan.GetMaxErrDistUsingEnum();
    EXPECT_EQ(maxErr, maxErrEnum);}

    {NetMan net1("./als/tests/benchmarks/am8_r54_MAXED9_s526_d38.blif");
    ErrMan errMan(net0, net1);
    ll maxErr = static_cast<ll>(errMan.ComputeMaxErr(METR_TYPE::MAXED));
    ll maxErrEnum = errMan.GetMaxErrDistUsingEnum();
    EXPECT_EQ(maxErr, maxErrEnum);}

    {NetMan net1("./als/tests/benchmarks/am8_r55_MAXED12_s524_d38.blif");
    ErrMan errMan(net0, net1);
    ll maxErr = static_cast<ll>(errMan.ComputeMaxErr(METR_TYPE::MAXED));
    ll maxErrEnum = errMan.GetMaxErrDistUsingEnum();
    EXPECT_EQ(maxErr, maxErrEnum);}

    {NetMan net1("./als/tests/benchmarks/am8_r57_MAXED13_s517_d38.blif");
    ErrMan errMan(net0, net1);
    ll maxErr = static_cast<ll>(errMan.ComputeMaxErr(METR_TYPE::MAXED));
    ll maxErrEnum = errMan.GetMaxErrDistUsingEnum();
    EXPECT_EQ(maxErr, maxErrEnum);}

    {NetMan net1("./als/tests/benchmarks/am8_r58_MAXED21_s502_d38.blif");
    ErrMan errMan(net0, net1);
    ll maxErr = static_cast<ll>(errMan.ComputeMaxErr(METR_TYPE::MAXED));
    ll maxErrEnum = errMan.GetMaxErrDistUsingEnum();
    EXPECT_EQ(maxErr, maxErrEnum);}

    {NetMan net1("./als/tests/benchmarks/am8_r17_MAXED25_s421_d35.blif");
    ErrMan errMan(net0, net1);
    ll maxErr = static_cast<ll>(errMan.ComputeMaxErr(METR_TYPE::MAXED));
    ll maxErrEnum = errMan.GetMaxErrDistUsingEnum();
    EXPECT_EQ(maxErr, maxErrEnum);}

    {NetMan net1("./als/tests/benchmarks/am8_r64_MAXED44_s479_d38.blif");
    ErrMan errMan(net0, net1);
    ll maxErr = static_cast<ll>(errMan.ComputeMaxErr(METR_TYPE::MAXED));
    ll maxErrEnum = errMan.GetMaxErrDistUsingEnum();
    EXPECT_EQ(maxErr, maxErrEnum);}

    GlobStopAbc();
}

