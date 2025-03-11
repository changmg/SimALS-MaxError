/**
 * @file sat_wrapper.h
 * @author Chang Meng (chang.meng@epfl.ch)
 * @brief Wrapper for the cryptominisat SAT solver
 * @date 2024-12-29
 * 
 */

#pragma once


#include <fmt/core.h>
#include "cryptominisat.h"


/**
 * @brief Self-defined fmt formatters
 */
template <>
struct fmt::formatter<CMSat::lbool> {
    constexpr auto parse(fmt::format_parse_context& ctx) {
        return ctx.begin();
    }
    template <typename FormatContext>
    auto format(const CMSat::lbool& value, FormatContext& ctx) {
        std::string str;
        if (value == CMSat::l_False)
            str = "False";
        else if (value == CMSat::l_True)
            str = "True";
        else
            str = "Undef";
        return fmt::format_to(ctx.out(), "{}", str);
    }
};
template <>
struct fmt::formatter<CMSat::Lit> {
    constexpr auto parse(fmt::format_parse_context& ctx) {
        return ctx.begin();
    }
    template <typename FormatContext>
    auto format(const CMSat::Lit& value, FormatContext& ctx) {
        return fmt::format_to(ctx.out(), "{}{}", value.sign()? "-": "", value.var());
    }
};
template <>
struct fmt::formatter<std::vector<CMSat::Lit>> {
    constexpr auto parse(fmt::format_parse_context& ctx) {
        return ctx.begin();
    }
    template <typename FormatContext>
    auto format(const std::vector<CMSat::Lit>& value, FormatContext& ctx) {
        std::string str;
        for (const auto& lit: value)
            str += fmt::format("{} ", lit);
        return fmt::format_to(ctx.out(), "{}", str);
    }
};


/**
 * @brief Solve the SAT problem defined by pSolver
 * 
 * @return CMSat::lbool  the result of the SAT problem
 */
static inline CMSat::lbool SolveSat_(CMSat::SATSolver& solver, bool printTime = false) {
    auto startTime = std::chrono::high_resolution_clock::now();
    auto res = solver.solve();
    if (printTime) 
        PrintRuntime(startTime, "SAT instance solving");
    return res;
}


/**
 * @brief Solve the SAT problem defined by pSolver
 * @brief If the result is SAT, return the counter-example
 * 
 * @param ids              the IDs of the variables to extract the counter-example
 * @param counterExample   the counter-example on the variables (if the result is SAT)
 * @param printTime        whether to print the runtime
 * @retval counterExample  the counter-example on the variables (if the result is SAT)
 * @return CMSat::lbool    the result of the SAT problem
 */
static inline CMSat::lbool SolveSatAndGetCountEx(CMSat::SATSolver& solver, const std::vector<CMSat::Lit>& assumpts, const std::vector<int>& cnfVarIds, std::vector<int>& counterExample, bool printTime = false) {
    auto startTime = std::chrono::high_resolution_clock::now();
    CMSat::lbool res;
    if (assumpts.empty())
        res = solver.solve();
    else
        res = solver.solve(&assumpts);
    if (printTime) 
        PrintRuntime(startTime, "SAT instance solving");
    if (res == CMSat::l_True) {
        counterExample.clear();
        for (auto cnfVarId: cnfVarIds)
            counterExample.emplace_back(solver.get_model()[cnfVarId] == CMSat::l_True? 1: 0);
    }
    return res;
}


/**
 * @brief Print the clauses in the SAT solver
 * 
 * @param solver  the SAT solver
 * @return void
 */
static inline void PrintSolverClauses(CMSat::SATSolver& solver) {
    fmt::print("Clauses in the SAT solver:\n");
    solver.start_getting_constraints(false);
    std::vector<CMSat::Lit> clause; 
    bool is_xor = false, rhs = false, ret = true;
    while (ret) {
        ret = solver.get_next_constraint(clause, is_xor, rhs);
        if (!ret) 
            break;
        assert(!is_xor);
        fmt::print("Clause: {}\n", clause);
    }
    solver.end_getting_constraints();
}