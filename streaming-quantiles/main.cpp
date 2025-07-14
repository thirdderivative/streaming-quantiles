#include <algorithm>
#include <bit>
#include <cassert>
#include <cinttypes>
#include <cmath>
#include <fmt/core.h>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <random>
#include <vector>

#include <cereal/archives/binary.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>

#include "relative_error_quantiles_sketch.h"

#define Type std::string

std::string GenerateKey() {
  static std::random_device rd;
  static std::mt19937_64 gen(rd());
  static std::uniform_int_distribution<uint64_t> distrib(
      0, std::numeric_limits<uint64_t>::max());

  uint64_t r1 = distrib(gen);
  uint64_t r2 = distrib(gen);
  uint64_t r3 = distrib(gen);
  uint64_t r4 = distrib(gen);
  uint64_t r5 = distrib(gen);
  return fmt::format("{:016x}:{:016x}:{:016x}:{:016x}:{:016x}", r1, r2, r3, r4,
                     r5);
}

uint32_t GenerateInt() {
  static std::random_device rd;
  static std::mt19937_64 gen(rd());
  static std::uniform_int_distribution<uint64_t> distrib(
      0, std::numeric_limits<uint64_t>::max());
  return distrib(gen);
}

int main(int argc, char **argv) {
  RelativeErrorQuantilesSketchOptions options = {
      // Rough estimate of the number of elements in input set.
      .n = 1'000'000'000,
      // Must be an even integer
      .k = 16384};
  assert(options.k % 2 == 0);

  RelativeErrorQuantilesSketch<Type> sketch(options);

  fmt::print("Attempting to insert {} keys\n", options.n);
  std::vector<Type> keys;
  for (uint64_t i = 0; i < options.n; ++i) {
    Type line = GenerateKey();
    // keys.push_back(line);
    sketch.Insert(line, 0);
  }
  fmt::print("Inserted {} keys\n", options.n);
  sketch.Close();

  const int numQuantiles = 1000;
  auto quantiles = sketch.Quantiles(numQuantiles);
  double error = 0.0;
  for (const auto &q : quantiles) {
    // TODO(tjackson): Modify EstimateRank() to get the weight for a quantile
    // then measure error as (weight - quantile fraction).
    const double e = q.cumulative_weight -
                     ((static_cast<double>(q.quantile) / numQuantiles) *
                      sketch.TotalWeight());
    fmt::print("Quantile {} at item {} with cumulative weight {} total weight "
               "{} error {}\n",
               q.quantile, q.item, q.cumulative_weight, sketch.TotalWeight(),
               e);
    error += std::pow(e, 2);
  }
  double rmse = std::sqrt(error / numQuantiles);
  fmt::print("RMSE: {}\n", rmse);
  return 0;
}
