/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2022 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Leona Gottwald and Michael    */
/*    Feldmeier                                                          */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_SOLUTION_POOL_H_
#define HIGHS_SOLUTION_POOL_H_

#include <vector>

#include "mip/HighsMipSolver.h"
#include "util/HighsRbTree.h"

class HighsSolutionPool;

namespace highs {
template <>
struct RbTreeTraits<HighsSolutionPool> {
  using KeyType = std::pair<double, HighsInt>;
  using LinkType = HighsInt;
};

}  // namespace highs

class HighsSolutionPool : public highs::CacheMinRbTree<HighsSolutionPool> {
 public:
  struct Solution {
    double reduced_objective_value;
    double original_objective_value;
    std::vector<double> reduced_values;
    std::vector<double> original_values;
    highs::RbTreeLinks<HighsInt> links;
  };

 private:
  HighsMipSolver& mipsolver;
  std::vector<Solution> solutions;
  HighsInt root = -1;
  HighsInt best = -1;

 public:
  highs::RbTreeLinks<HighsInt>& getRbTreeLinks(HighsInt node) {
    return solutions[node].links;
  }

  const highs::RbTreeLinks<HighsInt>& getRbTreeLinks(HighsInt node) const {
    return solutions[node].links;
  }

  std::pair<double, HighsInt> getKey(HighsInt node) const {
    return std::make_pair(solutions[node].reduced_objective_value, node);
  }

  HighsSolutionPool(HighsMipSolver& mipsolver)
      : highs::CacheMinRbTree<HighsSolutionPool>(root, best),
        mipsolver(mipsolver) {
    solutions.reserve(mipsolver.options_mip_->mip_num_solutions);
  }

  void addSolution(double reduced_objective_value,
                   const std::vector<double>& reduced_values);

  void untransform();

  void retransform();
};

#endif
