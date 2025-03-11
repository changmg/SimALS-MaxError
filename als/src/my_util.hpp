/**
 * @file my_util.h
 * @author Chang Meng (chang.meng@epfl.ch)
 * @brief Util functions
 * @date 2025-01-21
 * 
 */
#pragma once


#include "header.h"


// constant values
const double EPSILON = 1e-8;
const double DELAY_TOL = 1e-3;
const double AREA_TOL = 1e-3;
const std::string HALF_DASH_LINE = "--------------------";
const std::string DASH_LINE = "----------------------------------------";


/**
 * @brief Check if a path exists
 * 
 * @param _path  path to be checked
 * @return bool  true if the path exists
 */
static inline bool IsPathExist(const std::string & _path) {
    std::filesystem::path path(_path);
    return std::filesystem::exists(path); 
}


/**
 * @brief Create a directory if it does not exist
 * 
 * @param _path   path to the directory
 * @return void
 */
static inline void CreateDir(const std::string& _path) {
    std::filesystem::path path(_path);
    if (!std::filesystem::exists(path)) {
        std::filesystem::create_directories(path);
        fmt::print("Create directory {}\n", _path);
    }
}


/**
 * @brief Fix the path by adding a slash at the end
 * 
 * @param _path   path to be fixed
 * @retval _path  fixed path
 * @return void
 */
static inline void FixPath(std::string& _path) {
    if (_path.back() != '/')
        _path += "/";
}


// util of float-point numbers
static inline bool DoubleEqual     (const double a, const double b, double epsilon = EPSILON) {return fabs(a - b) < epsilon;}
static inline bool DoubleGreat     (const double a, const double b, double epsilon = EPSILON) {return a - b >= epsilon;}
static inline bool DoubleGreatEqual(const double a, const double b, double epsilon = EPSILON) {return a - b > -epsilon;}
static inline bool DoubleLess      (const double a, const double b, double epsilon = EPSILON) {return a - b <= -epsilon;}
static inline bool DoubleLessEqual (const double a, const double b, double epsilon = EPSILON) {return a - b < epsilon;}


// util of bit operation
template <typename T>
concept INT_NUMB = (
    std::is_same_v<T, uint32_t> ||
    std::is_same_v<T, ull> ||
    std::is_same_v<T, BigInt>
);


template <INT_NUMB T>
inline void SetBit(T & x, const int bit) {
    x |= (T(1) << bit);
}


template <INT_NUMB T>
inline void ResetBit(T & x, const int bit) {
    x &= ~(T(1) << bit);
}


template <INT_NUMB T>
inline void ChangeBit(T & x, const int bit, const bool val) {
    x |= (T(val) << bit);
}


template <INT_NUMB T>
inline bool GetBit(const T x, const int bit) {
    return static_cast<bool>((x >> bit) & T(1));
}


template <> 
inline void boost::to_block_range(const boost::dynamic_bitset<>& b, std::tuple<int, ull&> param)
{
    int iBlock = std::get<0>(param);
    std::get<1>(param) = b.m_bits[iBlock];
    return;
}


static inline ull GetBlockFromDynBitset(const boost::dynamic_bitset<>& b, int iBlock) {
    ull res = 0;
    boost::to_block_range(b, std::make_tuple(iBlock, std::ref(res)));
    return res;
}


/**
 * @brief Format the dynamic bit vector
 */
template <> 
struct fmt::formatter<BitVect> {
    constexpr auto parse(fmt::format_parse_context& ctx) {
        return ctx.begin();
    }
    template <typename FormatContext>
    auto format(const BitVect& value, FormatContext& ctx) {
        std::string str;
        // print from MSB to LSB
        for (int i = static_cast<int>(value.size()) - 1; i >= 0; --i)
            str += fmt::format("{}", static_cast<int>(value[i]));
        return fmt::format_to(ctx.out(), "{}", str);
    }
};


/**
 * @brief Execute system command
 * 
 * @param cmd       the command to be executed
 * @return int      1 if the command is executed successfully; 0 otherwise
 */
static inline int ExecSystComm(const std::string& cmd) {
    pid_t status;
    fmt::print("Execute system command: {}\n", cmd);
    status = system(cmd.c_str());
    if (status == -1) {
        fmt::print(stderr, "System error!\n");
        assert(0);
        return 0;
    }
    else {
        if (WIFEXITED(status)) {
            if (!WEXITSTATUS(status)) {
                fmt::print("Run shell script successfully.\n");
                return 1;
            }
            else {
                fmt::print("Run shell script fail, script exit code: {}\n", WEXITSTATUS(status));
                return 0;
            }
        }
        else {
            fmt::print("Exit status = {}\n", WEXITSTATUS(status));
            return 0;
        }
    }
}


/**
 * @brief Print the runtime
 * 
 * @param startTime   the start time
 * @param info        the information to be printed
 */
static inline void PrintRuntime(const std::chrono::time_point<std::chrono::high_resolution_clock> startTime, const std::string& info = "current") {
    auto currTime = std::chrono::high_resolution_clock::now();
    fmt::print("{} runtime = {}ms\n", info, std::chrono::duration_cast<std::chrono::milliseconds>(currTime - startTime).count());
}
