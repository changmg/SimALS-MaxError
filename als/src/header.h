/**
 * @file header.h
 * @author Chang Meng (chang.meng@epfl.ch)
 * @brief General header file
 * @date 2025-01-25
 * 
 */

#pragma once


#include <algorithm>
#include <boost/dynamic_bitset.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/timer/progress_display.hpp>
#include <boost/timer/timer.hpp>
#include <boost/random.hpp>
#include <chrono>
#include <concepts>
#include <filesystem>
#include <fmt/core.h>
#include <fmt/os.h>
#include <iostream>
#include <limits.h>
#include <list>
#include <memory>
#include <mutex>
#include <omp.h>
#include <queue>
#include <random>
#include <regex>
#include <set>
#include <sstream>
#include <stdint.h>
#include <thread>
#include <tuple>
#include <unordered_set>
#include <unordered_map>


using ll = int64_t;
using ull = uint64_t;
using BigFlt = boost::multiprecision::cpp_dec_float_100;
using BigInt = boost::multiprecision::int512_t;
using IntPair = std::pair<int, int>;
using IntTuple = std::tuple<int, int, int>;
using IntVect = std::vector<int>;
using IntSet = std::unordered_set<int>;
using Int2DVect = std::vector<IntVect>;
using LLVect = std::vector<ll>;
using LL2DVect = std::vector<LLVect>;
using DblVect = std::vector<double>;
using BigIntVect = std::vector<BigInt>;
using BigInt2DVect = std::vector<BigIntVect>;
using BitVect = boost::dynamic_bitset<ull>;


const ll LL_MAX = std::numeric_limits<ll>::max();
const BigInt MAX_BIGINT = std::numeric_limits<BigInt>::max();

