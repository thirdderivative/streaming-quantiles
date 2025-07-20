#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <fmt/core.h>
#include <functional>
#include <map>
#include <vector>

#include "compactor.h"

struct RelativeErrorQuantilesSketchOptions {
  uint64_t n;
  uint64_t k;
};

template <typename T> class RelativeErrorQuantilesSketch {
public:
  explicit RelativeErrorQuantilesSketch(
      const RelativeErrorQuantilesSketchOptions &options)
      : options_(options), H_(0), compactors_(std::vector<Compactor<T>>()),
        total_weight_(0) {
    fmt::print("Creating Relative error quantiles sketch with parameters k {} "
               "n {}...\n",
               options.k, options.n);
    compactors_.push_back(Compactor<T>(options_.k, options_.n, 0));
  }

  void Insert(const T &element, const uint64_t h) {
    if (H_ < h) {
      fmt::print("Sketch H {} < intended h {}, creating new compactor\n", H_,
                 h);
      H_ = h;
      compactors_.push_back(Compactor<T>(options_.k, options_.n, h));
    }

    std::vector<T> output_stream = std::move(compactors_[h].Insert(element));
    std::for_each(output_stream.begin(), output_stream.end(),
                  [this, h](const T &t) -> void { Insert(t, h + 1); });
  }

  struct WeightedElement {
    T item;
    double weight;
  };

  void Close() {
    assert(H_ + 1 == compactors_.size());

    // Weight of an item is 2^h, where h is the position the compactor has in
    // the overall hierarchy.
    uint64_t h = 0;
    for (auto &compactor : compactors_) {
      const double weight = std::pow(2, h);
      compactor.buffer.shrink_to_fit();
      std::transform(compactor.buffer.begin(), compactor.buffer.end(),
                     std::back_inserter(weighted_elements_),
                     [weight](const T &t) -> WeightedElement {
                       assert(!t.empty());
                       return WeightedElement{.item = t, .weight = weight};
                     });
      ++h;
    };
    std::sort(weighted_elements_.begin(), weighted_elements_.end(),
              [](const WeightedElement &w1, const WeightedElement &w2) -> bool {
                return w1.item < w2.item;
              });
    total_weight_ =
        std::accumulate(weighted_elements_.begin(), weighted_elements_.end(),
                        static_cast<double>(0.0),
                        [](double d, const WeightedElement &element) -> double {
                          return d + element.weight;
                        });
    /*for (const auto &element : weighted_elements_) {
      fmt::print("WeightedElement[ .item = {}, .weight = {}]\n", element.item,
                 element.weight);
    }*/
  }

  [[nodiscard]] double EstimateRank(const T &item) const {
    auto it = std::lower_bound(
        weighted_elements_.begin(), weighted_elements_.end(), item,
        [](const WeightedElement &element, const T &t) -> bool {
          return element.item < t;
        });
    auto index = std::distance(weighted_elements_.begin(), it);
    fmt::print("Found {} elements smaller than item {} out of {}\n", index,
               item, weighted_elements_.size());
    double item_weight = std::accumulate(
        weighted_elements_.begin(), it, static_cast<double>(0.0),
        [](double d, const WeightedElement &element) -> double {
          return d + element.weight;
        });
    fmt::print("Item weight {} total weight {}\n", item_weight, total_weight_);
    return item_weight;
  }

  struct Quantile {
    int quantile;
    T item;
    double cumulative_weight;
  };

  [[nodiscard]] std::vector<Quantile> Quantiles(int n) {
    std::vector<Quantile> quantiles;

    int current_quantile = 1;
    double current_total_weight = 0.0;
    int index = 0;
    for (const auto &element : weighted_elements_) {
      current_total_weight += element.weight;
      if (current_total_weight / total_weight_ >=
          static_cast<double>(current_quantile) / n) {
        Quantile quantile = {.quantile = current_quantile,
                             .item = element.item,
                             .cumulative_weight = current_total_weight};
        quantiles.push_back(quantile);
        fmt::print(
            "Found {} out of {} quantiles at item {} [index {} current total "
            "weight {}, total weight {}]\n",
            current_quantile, n, element.item, index, current_total_weight,
            total_weight_);
        ++current_quantile;
      }
      ++index;
    }
    return quantiles;
  }

  [[nodiscard]] uint64_t Depth() const { return H_; }

  [[nodiscard]] double TotalWeight() const { return total_weight_; }

  void Print() const {
    fmt::print("Sketch n {} k {} H {}\n", options_.n, options_.k, H_);
    std::for_each(
        compactors_.begin(), compactors_.end(),
        [](const Compactor<T> &compactor) -> void { compactor.Print(); });
  }

private:
  const RelativeErrorQuantilesSketchOptions options_;
  uint64_t H_;
  std::vector<Compactor<T>> compactors_;
  std::vector<WeightedElement> weighted_elements_;
  double total_weight_;
};
