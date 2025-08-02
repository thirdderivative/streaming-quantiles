#include "relative_error_quantiles_sketch.h"
#include <gtest/gtest.h>

TEST(RelativeErrorQuantilesSketchTest, InitialState) {
  RelativeErrorQuantilesSketchOptions options;
  options.k = 2;
  options.n = 8;
  RelativeErrorQuantilesSketch<int> sketch(options);
  ASSERT_EQ(sketch.Depth(), 0);
  ASSERT_EQ(sketch.TotalWeight(), 0);
}

TEST(RelativeErrorQuantilesSketchTest, Insert) {
  RelativeErrorQuantilesSketchOptions options;
  options.k = 2;
  options.n = 8;
  RelativeErrorQuantilesSketch<std::string> sketch(options);

  sketch.Insert("a", 0);
  sketch.Insert("b", 0);
  sketch.Insert("c", 0);

  ASSERT_EQ(sketch.Depth(), 0);
  // After 3 inserts, no compaction should have happened, so total weight is still 0 until Close() is called.
  ASSERT_EQ(sketch.TotalWeight(), 0);

  // This should trigger a new compactor to be created at level 1
  sketch.Insert("d", 1);
  ASSERT_EQ(sketch.Depth(), 1);
}

TEST(RelativeErrorQuantilesSketchTest, RankAndQuantiles) {
  RelativeErrorQuantilesSketchOptions options;
  options.k = 4;
  options.n = 100;
  RelativeErrorQuantilesSketch<int> sketch(options);

  for (int i = 1; i <= 100; ++i) {
    sketch.Insert(i, 0);
  }

  sketch.Close();

  // With k=4 and n=100, the buffer size of the first compactor is 4 * 2 * ceil(log2(100/4)) = 8 * ceil(log2(25)) = 8 * 5 = 40.
  // We inserted 100 items, so there will be compactions.
  // The total weight should be 100.
  ASSERT_EQ(sketch.TotalWeight(), 100);

  // The rank of 50 should be close to 50.
  // The exact value depends on the random choices during compaction, but it should be within a certain error bound.
  // For this test, we'll just check if it's in a reasonable range.
  // Since we insert numbers from 1 to 100, the rank of 50 should be the sum of weights of numbers less than 50.
  // Let's check the rank of 51, which should be approximately 50.
  double rank51 = sketch.EstimateRank(51);
  ASSERT_NEAR(rank51, 50, 15);

  // Check the 2-quantiles (median).
  std::vector<RelativeErrorQuantilesSketch<int>::Quantile> quantiles =
      sketch.Quantiles(2);
  ASSERT_EQ(quantiles.size(), 1);
  ASSERT_EQ(quantiles[0].quantile, 1);
  // The median should be around 50.
  ASSERT_NEAR(quantiles[0].item, 50, 15);
}

TEST(RelativeErrorQuantilesSketchTest, CompactionAndGrowth) {
  RelativeErrorQuantilesSketchOptions options;
  options.k = 2;
  options.n = 8;
  RelativeErrorQuantilesSketch<int> sketch(options);

  for (int i = 0; i < 1000; ++i) {
    sketch.Insert(i, 0);
  }

  ASSERT_GT(sketch.Depth(), 1);
}
