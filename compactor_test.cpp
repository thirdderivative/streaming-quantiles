#include "compactor.h"
#include <gtest/gtest.h>

TEST(CompactorTest, BufferSize) {
  Compactor<int> compactor(16, 1024, 0);
  ASSERT_EQ(compactor.max_buffer_size, 192);
}

TEST(CompactorTest, Insert) {
  Compactor<int> compactor(16, 1024, 0);
  compactor.Insert(1);
  compactor.Insert(2);
  compactor.Insert(3);
  ASSERT_EQ(compactor.buffer.size(), 3);
  ASSERT_EQ(compactor.buffer[0], 1);
  ASSERT_EQ(compactor.buffer[1], 2);
  ASSERT_EQ(compactor.buffer[2], 3);
}

TEST(CompactorTest, Compaction) {
  Compactor<int> compactor(2, 8, 0);
  // max_buffer_size should be 2 * 2 * ceil(log2(8/2)) = 4 * ceil(log2(4)) = 4 * 2 = 8
  ASSERT_EQ(compactor.max_buffer_size, 8);

  for (int i = 0; i < 8; ++i) {
    compactor.Insert(i);
  }
  ASSERT_EQ(compactor.buffer.size(), 8);

  compactor.C = 1;
  // The next insert should trigger a compaction.
  std::vector<int> output = compactor.Insert(8);
  std::cout << "Buffer size after compaction: " << compactor.buffer.size() << std::endl;
  std::cout << "Output size after compaction: " << output.size() << std::endl;
  ASSERT_EQ(compactor.buffer.size(), 5); // 8 - 2*2 + 1 = 5
  ASSERT_EQ(output.size(), 2); // 4/2 = 2
}

TEST(CompactorTest, CompactionInitialState) {
  Compactor<int> compactor(2, 8, 0);
  ASSERT_EQ(compactor.max_buffer_size, 8);
  ASSERT_EQ(compactor.buffer.size(), 0);
  ASSERT_EQ(compactor.C, 0);
}
