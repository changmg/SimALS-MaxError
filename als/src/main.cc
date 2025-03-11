/**
 * @file main.cc
 * @author Chang Meng (chang.meng@epfl.ch)
 * @brief Main function
 * @date 2025-01-07
 * 
 */
#include "cmdline.hpp"
#include "als.h"


using cmdline::parser;
using std::string;


/**
 * @brief Command line parser
 * 
 * @param argc     number of arguments
 * @param argv     arguments
 * @return parser  parsed options
 */
parser CommPars(int argc, char * argv[]) {
    parser option;
    option.add<string>("accCirc", 'i', "path to accurate circuit", true);
    option.add<string>("standCellLib", 'l', "path to standard cell library", false, "./input/standard-cell/nangate_45nm_typ.lib");
    option.add<string>("outpPath", 'o', "path to approximate circuits", false, "./tmp/");
    option.add<string>("metrType", 'm', "error metric type: MAXED, MAXHD", false, "MAXED");
    option.add<unsigned>("seed", 's', "seed", false, 199608224);
    option.add<int>("nFrame", 'f', "#simulation patterns for maximum error estimation", false, 8192);
    option.add<int>("fUseMecals1_0", 'u', "use MECALS 1.0 (DATE'23 version)", false, 0);
    option.add<int>("fFastFlow", '\0', "use fast flow for EPFL large benchmarks", false, 0);
    option.add<double>("exactPBDPerc", 'p', "proportion of exact PBD (only used in MECALS 1.0)", false, 1.0);
    option.add<ll>("errUppBound", 'e', "upper bound of maximum error", false, 64);
    option.parse_check(argc, argv);
    return option;
}


/**
 * @brief Approximate logic synthesis
 * 
 * @param option  parsed options
 * @return void
 */
void ALS(const parser& option) {
    // print current date and time
    auto now = std::chrono::system_clock::now();
    std::time_t now_time = std::chrono::system_clock::to_time_t(now);
    fmt::print("Current date and time: {}\n", std::ctime(&now_time));

    // load configurations
    auto accCirc = option.get<string>("accCirc");
    auto standCellLib = option.get<string>("standCellLib");
    auto outpPath = option.get<string>("outpPath");
    // auto lacType = option.get<string>("lacType");
    auto metrType = option.get<string>("metrType");
    auto seed = option.get<unsigned> ("seed");
    auto nFrame = option.get<int>("nFrame");
    auto fUseMecals1_0 = option.get<int>("fUseMecals1_0");
    auto fFastFlow = option.get<int>("fFastFlow");
    auto exactPBDPerc = option.get<double>("exactPBDPerc");
    auto errUppBound = option.get<ll>("errUppBound");

    // extract circuit name
    if (!accCirc.ends_with(".blif") && !accCirc.ends_with(".aig")) {
        fmt::print(stderr, "Error: the accurate circuit should be in BLIF or AIG format.\n");
        assert(0);
    }
    if (!IsPathExist(accCirc)) {
        fmt::print(stderr, "Error: the accurate circuit file {} does not exist.\n", accCirc);
        assert(0);
    }
    std::filesystem::path accCircPath(accCirc);
    string accCircName = accCircPath.stem().string();
    fmt::print("accurate circuit: {}\n", accCirc);

    // fix & create output path
    FixPath(outpPath);
    CreateDir(outpPath);
    outpPath += accCircName + "_";

    // set ALS options
    ALSOpt alsOpt(metrType, seed, nFrame, fUseMecals1_0, exactPBDPerc, errUppBound, outpPath);
    alsOpt.ProcSeed();
    fmt::print("{}", alsOpt);

    // read standard cell library and accurate circuit
    AbcMan abc;
    abc.ReadStandCell(standCellLib);
    NetMan accNet(accCirc);

    // approximate logic synthesis
    ALSMan alsMan(accNet, alsOpt);
    if (alsOpt.fUseMecals1_0) 
        alsMan.Run_v1(); 
    else if (fFastFlow)
        alsMan.Run_FastFlow();
    else
        alsMan.Run_v2();
}


/**
 * @brief Main function
 * 
 * @param argc  number of arguments
 * @param argv  arguments
 * @return int  0 if success
 */
int main(int argc, char *argv[]) {
    parser option = CommPars(argc, argv);
    GlobStartAbc();
    ALS(option);
    GlobStopAbc();
    return 0;
}