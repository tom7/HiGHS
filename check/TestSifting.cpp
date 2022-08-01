#include <algorithm>

#include "Highs.h"
#include "catch.hpp"

using std::min;

const bool dev_run = true;
const double inf = kHighsInf;

// No commas in test case name.
TEST_CASE("test-sifting", "[highs_sifting]") {
  const bool test_sifting = true;
  if (test_sifting) {
    HighsLp lp;
    Highs highs;
    highs.setOptionValue("output_flag", dev_run);
    highs.setOptionValue("presolve", kHighsOffString);
    highs.setOptionValue("sifting_strategy", kSiftingStrategyOff);
    const double profile = 10;
    const HighsInt num_row = 1;
    const HighsInt num_col = num_row * profile;
    const HighsInt min_nz_per_col = std::min((HighsInt)3, num_row);
    const double density = std::max(0.2, (1.0 * min_nz_per_col) / num_row);
    printf("Density = %g\n", density);
    const HighsInt row_count = num_col * density;
    
    HighsRandom random;
    lp.col_lower_.assign(num_col, 0);
    lp.col_upper_.assign(num_col, inf);
    for (HighsInt iCol = 0; iCol < num_col; iCol++) {
      double cost = 1 + random.fraction();
      if (num_row == 1) cost = 0.1 * (num_col - iCol);
      lp.col_cost_.push_back(cost);
    }
    std::vector<bool> check_index;
    check_index.assign(num_col, false);
    printf("row_count = %d\n", (int)row_count);
    lp.a_matrix_.format_ = MatrixFormat::kRowwise;
    lp.a_matrix_.start_.clear();
    for (HighsInt iRow = 0; iRow < num_row; iRow++) {
      double row_sum = 0;
      HighsInt el = lp.a_matrix_.index_.size();
      lp.a_matrix_.start_.push_back(el);
      HighsInt num_nz = 0;
      for (;;) {
	HighsInt iCol = random.integer(num_col);
	if (check_index[iCol]) continue;
	lp.a_matrix_.index_.push_back(iCol);
	check_index[iCol] = true;
	double value = random.fraction();
	if (num_row == 1) value = 1.0;
	lp.a_matrix_.value_.push_back(value);
	row_sum += value;
	num_nz++;
	if (num_nz >= row_count) break;
      }
      lp.row_lower_.push_back(row_count * row_sum);
      lp.row_upper_.push_back(inf);
      for (HighsInt el = lp.a_matrix_.start_[iRow];
	   el < lp.a_matrix_.start_[iRow] + row_count; el++)
	check_index[lp.a_matrix_.index_[el]] = false;
      for (HighsInt iCol = 0; iCol < num_col; iCol++) assert(!check_index[iCol]);
    }
    lp.a_matrix_.start_.push_back(lp.a_matrix_.index_.size());
    
    lp.num_col_ = lp.col_lower_.size();
    lp.num_row_ = lp.row_lower_.size();
    assert(lp.num_col_ == num_col);
    assert(lp.num_row_ == num_row);
    
    highs.passModel(lp);
    highs.run();
    const double optimal_objective = highs.getInfo().objective_function_value;
    const HighsSolution& solution = highs.getSolution();
    HighsInt num_nonzero_col_value = 0;
    for (HighsInt iCol = 0; iCol < num_col; iCol++) {
      double value = solution.col_value[iCol];
      if (value > 1e-4) num_nonzero_col_value++;
    }
    HighsInt num_nonzero_slack = 0;
    for (HighsInt iRow = 0; iRow < num_row; iRow++) {
      double value = solution.row_value[iRow] - lp.row_lower_[iRow];
      if (value > 1e-4) num_nonzero_slack++;
    }
    printf("num_nonzero_col_value = %d/%d; num_nonzero_slack = %d/%d\n",
	   (int)num_nonzero_col_value, (int)num_col, (int)num_nonzero_slack,
	   (int)num_row);
    highs.setOptionValue("simplex_strategy", kSimplexStrategyChoose);
    highs.setOptionValue("sifting_strategy", kSiftingStrategyOn);
    highs.setBasis();
    highs.setSolution();
    highs.run();
  }
}

