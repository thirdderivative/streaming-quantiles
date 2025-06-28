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

#define Type std::string

bool RandomBoolean() {
  static std::mt19937 generator;
  static std::random_device rd;
  generator.seed(rd());
  std::bernoulli_distribution distribution(0.5);
  return distribution(generator);
}

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

template <typename T> struct Compactor {
  uint64_t n;               // Number of sections for compaction
  uint64_t k;               // Section size
  uint64_t max_buffer_size; // Max buffer size
  uint64_t C;               // Compaction schedule
  uint64_t h;               // Position in the compaction hierarchy
  std::vector<T> buffer;    // Items

  Compactor(uint64_t k, uint64_t n, uint64_t h) : n(n), k(k), C(0), h(h) {
    const uint64_t m =
        std::ceil(std::log(static_cast<double>(n) / static_cast<double>(k)) /
                  std::log(2));
    max_buffer_size = 2 * k * m;
    fmt::print("{}: Creating compactor with buffer size {}\n", h,
               max_buffer_size);
  }

  std::vector<T> Insert(const T &element) {
    std::vector<T> output;
    if (buffer.size() == max_buffer_size) {
      int sections_to_compact = std::countr_one<uint64_t>(C) + 1;
      const uint64_t elements_to_compact = sections_to_compact * k;

      // Put the largest elements to compact at the back of the buffer by
      // figuring out where to pivot the buffer. See example at
      // https://en.cppreference.com/w/cpp/algorithm/partial_sort.html
      const uint64_t S = max_buffer_size - elements_to_compact + 1;
      std::partial_sort(buffer.rbegin(), buffer.rbegin() + S, buffer.rend(),
                        std::greater{});

      // Take even or odd indexes for compacted sections to add to the output
      // for pushing to the next compactor, then drop the compacted sections.
      bool even = RandomBoolean();
      uint64_t i = S - 1;
      if (!even && i % 2 == 0) {
        ++i;
      }
      while (i < max_buffer_size) {
        output.push_back(buffer[i]);
        i += 2;
      }

      // Clear elements after element S
      const uint64_t target_capacity = max_buffer_size - elements_to_compact;
      buffer.resize(target_capacity);

      // Resize the vector down so that we don't have extra memory growth.
      // Post "compaction" we should have max_buffer_size/2 elements in the
      // buffer and it should start growing again.
      buffer.shrink_to_fit();
      assert(buffer.size() == target_capacity);
      assert(buffer.capacity() == target_capacity);

      // Update the compactor schedule so we can "randomly" choose new
      // sections next time.
      ++C;
    }

    // Now that we're done compaction (if necessary) we can add the element to
    // the buffer.
    buffer.push_back(element);
    return output;
  };

  void Print() const {
    fmt::print("{}: Compactor n {} k {} max_buffer_size {} C {}\n", h, n, k,
               max_buffer_size, C);
    fmt::print("  Buffer:\n");
    std::for_each(buffer.begin(), buffer.end(),
                  [](const T &item) -> void { fmt::print("    {}\n", item); });
  }
};

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

  template <class Archive> void serialize(Archive &ar) {
    ar(H_, compactors_, options_);
  }

private:
  const RelativeErrorQuantilesSketchOptions options_;
  uint64_t H_;
  std::vector<Compactor<T>> compactors_;
  std::vector<WeightedElement> weighted_elements_;
  double total_weight_;
};

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
