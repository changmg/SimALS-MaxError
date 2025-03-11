/**
 * @file als.h
 * @author Chang Meng (chang.meng@epfl.ch)
 * @brief Approximate logic synthesis
 * 
 */
#pragma once


#include "header.h"
#include "my_abc.h"
#include "error.h"
#include "lac.h"
#include "my_util.hpp"


/**
 * @brief Options for approximate logic synthesis
 * 
 */
class ALSOpt {
public:
    // LAC_TYPE lacType;                 // LAC type
    METR_TYPE metrType;               // error metric type
    unsigned seed;                    // seed
    int nFrame;                       // number of simulation frames used for maximum error lower bound computation
    int fUseMecals1_0;                // flag of using MECALS 1.0 (DATE'23 version)
    int appResub_nFrame4ResubGen;     // Resub-based LAC: number of simulation frames for approximate resubstitution generation
    int appResub_maxLevelDiff;        // Resub-based LAC: maximum level difference between the target node and divisors when generating approximate resubstitutions
    int maxCandLacs;                  // Resub-based/SASIMI LAC: maximum number of candidate LACs
    double mecals1_exactPBDPerc;      // MECALS1.0: proportion of exact partial Boolean difference
    ll errUppBound;                   // upper bound of error
    std::string outpPath;             // output path

    explicit ALSOpt(const std::string& metr_type, unsigned _seed, int n_frame, int f_use_mecals1_0, double exact_pbd_perc, ll err_upp_bound, const std::string& outp_path): 
        // lacType(LAC_TYPE::CONSTANT),
        metrType(Str2MetrType(metr_type)),
        seed(_seed), 
        nFrame(n_frame),
        fUseMecals1_0(f_use_mecals1_0),
        // appResub_nFrame4ResubGen(64), 
        appResub_nFrame4ResubGen(4),
        appResub_maxLevelDiff(INT_MAX), 
        maxCandLacs(100000), 
        mecals1_exactPBDPerc(exact_pbd_perc),
        errUppBound(err_upp_bound),
        outpPath(outp_path)
    {
        if (errUppBound < 0) {
            fmt::print(stderr, "Error: errUppBound should be non-negative.\n");
            assert(0);
        }
    }

    ~ALSOpt() = default;

    inline void ProcSeed() {
        if (seed == 0) {
            boost::random::mt19937 rng(time(0));
            boost::uniform_int <> unif(INT_MIN, INT_MAX);
            seed = static_cast <unsigned> (unif(rng));
        }
    }
};
template <>
struct fmt::formatter<ALSOpt> {
    constexpr auto parse(fmt::format_parse_context& ctx) {
        return ctx.begin();
    }
    template <typename FormatContext>
    auto format(const ALSOpt& value, FormatContext& ctx) {
        std::string str = "---------- ALS options ----------\n";
        // str += fmt::format("lacType = {}\n", value.lacType);
        str += fmt::format("metrType = {}\n", value.metrType);
        str += fmt::format("seed = {}\n", value.seed);
        str += fmt::format("nFrame = {}\n", value.nFrame);
        str += fmt::format("fUseMecals1_0 = {}\n", value.fUseMecals1_0);
        str += fmt::format("appResub_nFrame4ResubGen = {}\n", value.appResub_nFrame4ResubGen);
        str += fmt::format("appResub_maxLevelDiff = {}\n", value.appResub_maxLevelDiff);
        str += fmt::format("maxCandLacs = {}\n", value.maxCandLacs);
        str += fmt::format("mecals1_exactPBDPerc = {}\n", value.mecals1_exactPBDPerc);
        str += fmt::format("errUppBound = {}\n", value.errUppBound);
        str += fmt::format("outpPath = {}\n", value.outpPath);
        str += fmt::format("--------------------\n");
        return fmt::format_to(ctx.out(), "{}", str);
    }
};


/**
 * @brief ALS manager
 */
class ALSMan {
private:
    const NetMan& accNet;                             // accurate network
    ALSOpt options;                                   // options
    std::shared_ptr<Simulator> pAccSmlt;              // simulator for accurate network
    std::shared_ptr<NetMan> pDevCompNet;              // network for maximum error computation (PI: accNet PI, appNetPI, ref_err; single PO = error > ref_err)
    std::shared_ptr<NetMan> pDevCompNetEmbErr;        // network for maximum error checking with embedded reference error (PI: accNet PI, appNetPI; single PO = error > ref_err; ref_err is a constant vector embedded in the network)

public:
    explicit ALSMan(NetMan& acc_net, ALSOpt& _options): accNet(acc_net), options(_options) {
        if (options.metrType == METR_TYPE::MAXHD) {
            if (options.errUppBound >= accNet.GetPoNum()) {
                fmt::print(stderr, "Error: the upper bound of the maximum Hamming distance should be less than the output width\n");
                assert(0);
            }
        }
        pAccSmlt = std::make_shared<Simulator>(accNet, options.seed, options.nFrame, DISTR_TYPE::UNIF);
        pAccSmlt->LogicSim();
        pDevCompNet = GenDevCompNet(options.metrType, accNet.GetPoNum());
        if (options.errUppBound > 0)
            pDevCompNetEmbErr = GenDevCompNetEmbedErrBound(pDevCompNet, accNet.GetPoNum(), options.errUppBound);
        else
            pDevCompNetEmbErr = nullptr;
    }
    ~ALSMan() = default;
    ALSMan(const ALSMan&) = delete;
    ALSMan(ALSMan&&) = delete;
    ALSMan& operator = (const ALSMan&) = delete;
    ALSMan& operator = (ALSMan&&) = delete;

    void Run_v1();
    void Run_v2();
    void Run_Trunc();
    void Run_FastFlow();
    void SimplifyWithSingleLac(LAC_TYPE lacType, NetMan& appNet, int& round, const std::chrono::time_point<std::chrono::high_resolution_clock>& startTime, bool inclConst, bool fSimplify);
    bool ApplyMultValidLacs(const LACMan& lacMan, NetMan& appNet, std::unordered_set<std::string>& lacBlackList, int& countExNum);
    bool ApplyMultValidLacs_NoSimPrune(const LACMan& lacMan, NetMan& appNet, std::unordered_set<std::string>& lacBlackList, int& countExNum);
    // bool ApplyMultValidLacs_IncSAT(const LACMan& lacMan, NetMan& appNet, std::unordered_set<std::string>& lacBlackList, int& countExNum);
    // bool ApplyMultValidLacsUsingBaseErr(const LACMan& lacMan, NetMan& appNet, std::unordered_set<std::string>& lacBlackList, int& countExNum);
};