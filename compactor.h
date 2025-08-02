#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <fmt/core.h>
#include <functional>
#include <random>
#include <vector>

bool RandomBoolean() {
  static std::mt19937 generator;
  static std::random_device rd;
  generator.seed(rd());
  std::bernoulli_distribution distribution(0.5);
  return distribution(generator);
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
      int sections_to_compact = 0;
      uint64_t c = C;
      while ((c & 1) == 1) {
        sections_to_compact++;
        c >>= 1;
      }
      sections_to_compact++;
      const uint64_t elements_to_compact = sections_to_compact * k;

      // Put the largest elements to compact at the back of the buffer by
      // figuring out where to pivot the buffer. See example at
      // https://en.cppreference.com/w/cpp/algorithm/partial_sort.html
      const uint64_t S = max_buffer_size - elements_to_compact;
      std::partial_sort(buffer.rbegin(), buffer.rbegin() + S, buffer.rend(),
                        std::greater{});

      // Take even or odd indexes for compacted sections to add to the output
      // for pushing to the next compactor, then drop the compacted sections.
      bool even = RandomBoolean();
      uint64_t i = S;
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
