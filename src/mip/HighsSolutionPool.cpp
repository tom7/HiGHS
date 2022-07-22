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

#include "mip/HighsSolutionPool.h"
#include "mip/HighsMipSolverData.h"

void HighsSolutionPool::addSolution(double reduced_objective_value,
                 const std::vector<double>& reduced_values) {
  HighsInt insertPos;

  if (solutions.size() == mipsolver.options_mip_->mip_num_solutions) {
    insertPos = this->last();
    if (solutions[insertPos].reduced_objective_value < reduced_objective_value)
      return;

    this->unlink(insertPos);

    solutions[insertPos].original_values.clear();
  } else {
    insertPos = solutions.size();
    solutions.emplace_back();
  }

  solutions.back().reduced_values = reduced_values;
  solutions.back().reduced_objective_value = reduced_objective_value;
}

void HighsSolutionPool::untransform() {
  for (Solution& sol : solutions) {
    if (sol.original_values.size() == mipsolver.orig_model_->num_col_) continue;

    HighsSolution solution;
    solution.col_value = sol.reduced_values;
    mipsolver.mipdata_->postSolveStack.undoPrimal(*mipsolver.options_mip_, solution);
    sol.original_values = std::move(solution.col_value);
  }
}

void HighsSolutionPool::retransform() {

}
