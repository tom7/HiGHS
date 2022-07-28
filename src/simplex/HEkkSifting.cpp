
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
/**@file simplex/HEkkSifting.cpp
 * @brief
 */
#include "lp_data/HighsLpSolverObject.h"
#include "simplex/HEkk.h"
#include "util/HighsSort.h"
#include "util/HighsRandom.h"
//#include "lp_data/HighsInfo.h"
#include "io/HMPSIO.h"
HighsStatus HEkk::sifting() {
  HighsStatus return_status = HighsStatus::kOk;
  std::vector<HighsInt> sifted_list;
  std::vector<bool> in_sifted_list;
  in_sifted_list.assign(lp_.num_col_, false);
  // Need to start from a primal feasible solution
  assert(info_.num_primal_infeasibilities >= 0);
  HighsInt sifted_list_max_count = lp_.num_row_;
  HEkk sifted_ekk_instance;
  HighsLp sifted_lp;
  HighsBasis sifted_basis;
  HighsSolution sifted_solution;
  HighsInfo sifted_highs_info;
  HighsOptions sifted_options = *options_;
  HighsLpSolverObject sifted_solver_object(sifted_lp, sifted_basis, sifted_solution,
    					   sifted_highs_info, sifted_ekk_instance,
    					   sifted_options, *timer_);
  // Prevent recursive sifting!
  sifted_options.sifting_strategy = kSiftingStrategyOff;
  sifted_options.log_dev_level = 3;
  //  sifted_options.highs_analysis_level = 4;

  HighsInt sifting_iter = 0;
  for (;;) {
    sifting_iter++;
    addToSiftedList(lp_.num_row_, sifted_solver_object, 
		    sifted_list, in_sifted_list);
    assert(okSiftedList(sifted_list, in_sifted_list));
    const bool write_lp = false;
    if (write_lp) {
      HighsModel model;
      model.lp_ = sifted_solver_object.lp_;
      writeModelAsMps(*options_, "sifted.mps", model);
    }
    sifted_ekk_instance.moveLp(sifted_solver_object);
    return_status = sifted_ekk_instance.solve();
    assert(return_status == HighsStatus::kOk);
    sifted_lp.moveBackLpAndUnapplyScaling(sifted_ekk_instance.lp_);

    updateIncumbentData(sifted_solver_object, sifted_list);
    sifted_ekk_instance.clear();
    if (sifting_iter>2) break;
  }
 
  assert(1 == 0);
  return return_status;
}

void HEkk::addToSiftedList(const HighsInt max_add_to_sifted_list,
			   HighsLpSolverObject& sifted_solver_object, 
			   std::vector<HighsInt>& sifted_list,
			   std::vector<bool>& in_sifted_list) {
  HighsLp& sifted_lp = sifted_solver_object.lp_;
  HighsBasis& sifted_basis = sifted_solver_object.basis_;
  const bool first_sifted_list = sifted_list.size() == 0;
  const bool primal_feasible =
    info_.num_primal_infeasibilities == 0;
  if (!primal_feasible) assert(first_sifted_list);
  if (first_sifted_list) {
    assert(sifted_lp.col_cost_.size() == 0);
    assert(sifted_lp.col_lower_.size() == 0);
    assert(sifted_lp.col_upper_.size() == 0);
    assert(sifted_lp.row_lower_.size() == 0);
    assert(sifted_lp.row_upper_.size() == 0);
  }
  if (primal_feasible) {
    for (HighsInt iX = 0; iX < sifted_list.size(); iX++) {
      HighsInt iCol = sifted_list[iX];
      double dual = info_.workCost_[iCol];
      for (HighsInt iEl = lp_.a_matrix_.start_[iCol]; iEl < lp_.a_matrix_.start_[iCol+1]; iEl++) {
	HighsInt iRow = lp_.a_matrix_.index_[iEl];
	dual += lp_.a_matrix_.value_[iEl] * info_.workDual_[lp_.num_col_ + iRow];
      }
      assert(std::fabs(dual-info_.workDual_[iCol]) < 1e-4);
    }
  }

  HighsInt num_add_to_sifted_list = 0;
  std::vector<HighsInt> heap_index;
  std::vector<double> heap_value;
  heap_index.push_back(0);
  heap_value.push_back(0);  
  for (HighsInt iCol = 0; iCol < lp_.num_col_; iCol++) {
    if (in_sifted_list[iCol]) continue;
    if (basis_.nonbasicFlag_[iCol] == 0) {
      // Basic, so in sifted list
      assert(first_sifted_list);
      num_add_to_sifted_list++;
      sifted_list.push_back(iCol);
      in_sifted_list[iCol] = true;
      sifted_lp.col_cost_.push_back(info_.workCost_[iCol]);
      sifted_lp.col_lower_.push_back(info_.workLower_[iCol]);
      sifted_lp.col_upper_.push_back(info_.workUpper_[iCol]);
      for (HighsInt iEl = lp_.a_matrix_.start_[iCol]; iEl < lp_.a_matrix_.start_[iCol+1]; iEl++) {
	sifted_lp.a_matrix_.index_.push_back(lp_.a_matrix_.index_[iEl]);
	sifted_lp.a_matrix_.value_.push_back(lp_.a_matrix_.value_[iEl]);
      }
      sifted_lp.a_matrix_.start_.push_back(sifted_lp.a_matrix_.index_.size());
      continue;
    }
    if (!primal_feasible) continue;
    // Nonbasic, and not in sifted list, so possible new entry
    const double check_dual = info_.workDual_[iCol];
    // Compute the dual for this column
    double dual = info_.workCost_[iCol];
    for (HighsInt iEl = lp_.a_matrix_.start_[iCol]; iEl < lp_.a_matrix_.start_[iCol+1]; iEl++) {
      HighsInt iRow = lp_.a_matrix_.index_[iEl];
      dual += lp_.a_matrix_.value_[iEl] * info_.workDual_[lp_.num_col_ + iRow];
    }
    if (first_sifted_list) assert(std::fabs(dual-check_dual) < 1e-4);
    // Determine the dual infeasibility for this column
    const double lower = info_.workLower_[iCol];
    const double upper = info_.workUpper_[iCol];
    double dual_infeasibility = 0;
    if (lower <= -kHighsInf && upper >= kHighsInf) {
      // Free: any nonzero dual value is infeasible
      dual_infeasibility = fabs(dual);
    } else {
      // Not free: any dual infeasibility is given by the dual value
      // signed by nonbasicMove
      dual_infeasibility = -basis_.nonbasicMove_[iCol] * dual;
    }
    if (dual_infeasibility > options_->dual_feasibility_tolerance) {
      heap_index.push_back(iCol);
      heap_value.push_back(dual_infeasibility);
    }
  }
  if (!primal_feasible) {
    // Construct an LP containing any basic columns and a random
    // collection of nonbasic columns
    HighsRandom random;
    for (;;) {
      HighsInt iCol = random.integer(lp_.num_col_);
      if (in_sifted_list[iCol]) continue;
      num_add_to_sifted_list++;
      sifted_list.push_back(iCol);
      in_sifted_list[iCol] = true;
      sifted_lp.col_cost_.push_back(info_.workCost_[iCol]);
      sifted_lp.col_lower_.push_back(info_.workLower_[iCol]);
      sifted_lp.col_upper_.push_back(info_.workUpper_[iCol]);
      for (HighsInt iEl = lp_.a_matrix_.start_[iCol]; iEl < lp_.a_matrix_.start_[iCol+1]; iEl++) {
	sifted_lp.a_matrix_.index_.push_back(lp_.a_matrix_.index_[iEl]);
	sifted_lp.a_matrix_.value_.push_back(lp_.a_matrix_.value_[iEl]);
      }
      sifted_lp.a_matrix_.start_.push_back(sifted_lp.a_matrix_.index_.size());
      if (num_add_to_sifted_list == max_add_to_sifted_list) break;
    }
  } else {
    HighsInt heap_num_en = heap_index.size() - 1;
    assert(heap_num_en >= 0);
    HighsInt sifted_list_count = sifted_list.size();
    if (heap_num_en == 0) {
      // Optimal!
      printf("Optimal!\n");
      assert(111==999);
    }
    // There are dual infeasibilities
    maxheapsort(&heap_value[0], &heap_index[0], heap_num_en);
    for (HighsInt iEl = heap_num_en; iEl >0; iEl--) {
      HighsInt iCol = heap_index[iEl];
      num_add_to_sifted_list++;
      sifted_list.push_back(iCol);
      in_sifted_list[iCol] = true;
      sifted_lp.col_cost_.push_back(info_.workCost_[iCol]);
      sifted_lp.col_lower_.push_back(info_.workLower_[iCol]);
      sifted_lp.col_upper_.push_back(info_.workUpper_[iCol]);
      for (HighsInt iEl = lp_.a_matrix_.start_[iCol]; iEl < lp_.a_matrix_.start_[iCol+1]; iEl++) {
	sifted_lp.a_matrix_.index_.push_back(lp_.a_matrix_.index_[iEl]);
	sifted_lp.a_matrix_.value_.push_back(lp_.a_matrix_.value_[iEl]);
      }
      sifted_lp.a_matrix_.start_.push_back(sifted_lp.a_matrix_.index_.size());
      if (num_add_to_sifted_list == max_add_to_sifted_list) break;
    }
  }
  std::vector<double> unsifted_row_activity;
  unsifted_row_activity.assign(lp_.num_row_, 0);
  for (HighsInt iCol = 0; iCol < lp_.num_col_; iCol++) {
    if (in_sifted_list[iCol]) continue;
    double value = info_.workValue_[iCol];
    for (HighsInt iEl = lp_.a_matrix_.start_[iCol]; iEl < lp_.a_matrix_.start_[iCol+1]; iEl++) {
      HighsInt iRow = lp_.a_matrix_.index_[iEl];
      unsifted_row_activity[iRow] += value * lp_.a_matrix_.value_[iEl];
    }
  }
  sifted_lp.row_lower_.resize(lp_.num_row_);
  sifted_lp.row_upper_.resize(lp_.num_row_);
  for (HighsInt iRow = 0; iRow < lp_.num_row_; iRow++) {
    HighsInt iVar = lp_.num_col_ + iRow;
    sifted_lp.row_lower_[iRow] = -info_.workUpper_[iVar] - unsifted_row_activity[iRow];
    sifted_lp.row_upper_[iRow] = -info_.workLower_[iVar] - unsifted_row_activity[iRow];
    /*
    printf("Row %3d has bounds [%11.4g, %11.4g] from  [%11.4g, %11.4g] and shift %11.4g\n",
	   (int)iRow, sifted_lp.row_lower_[iRow], sifted_lp.row_upper_[iRow],
	   info_.workLower_[iVar], info_.workUpper_[iVar], unsifted_row_activity[iRow]);
    */
  }
  sifted_lp.num_col_ = sifted_lp.col_lower_.size();
  sifted_lp.num_row_ = sifted_lp.row_lower_.size();
  assert(sifted_lp.num_col_ == sifted_list.size());
  assert(sifted_lp.num_row_ == lp_.num_row_);
  sifted_lp.setMatrixDimensions();
  
}


bool HEkk::okSiftedList(const std::vector<HighsInt>& sifted_list,
			const std::vector<bool>& in_sifted_list) {
  std::vector<bool> local_in_sifted_list = in_sifted_list;
  for (HighsInt iX = 0; iX < sifted_list.size(); iX++) {
    if (!local_in_sifted_list[sifted_list[iX]]) {
      printf("local_in_sifted_list[sifted_list[%d]] is false\n", (int)iX);
      return false;
    }
    local_in_sifted_list[sifted_list[iX]] = false;
  }
  for (HighsInt iCol = 0; iCol < sifted_list.size(); iCol++) {
    if (local_in_sifted_list[iCol]) {
      printf("local_in_sifted_list[%d] is true\n", (int)iCol);
      return false;
    }
  }
}

void HEkk::updateIncumbentData(const HighsLpSolverObject& sifted_solver_object, 
			       const std::vector<HighsInt>& sifted_list) {
  HighsLp& sifted_lp = sifted_solver_object.lp_;
  HighsBasis& sifted_basis = sifted_solver_object.basis_;
  HEkk& sifted_ekk_instance = sifted_solver_object.ekk_instance_;
  for (HighsInt iX = 0; iX < sifted_list.size(); iX++) {
    HighsInt iCol = sifted_list[iX];
    if (sifted_ekk_instance.basis_.nonbasicFlag_[iX]) {
      info_.workValue_[iCol] = sifted_ekk_instance.info_.workValue_[iX];
      info_.workDual_[iCol] = sifted_ekk_instance.info_.workDual_[iX];
      printf("Nonbasic: iX %2d: iCol = %3d: Value = %11.4g; Dual = %11.4g\n",
	   (int)iX, (int)iCol, info_.workValue_[iCol], info_.workDual_[iCol]);
    }
  }
  for (HighsInt iRow = 0; iRow < lp_.num_row_; iRow++) {
    HighsInt iX = sifted_ekk_instance.basis_.basicIndex_[iRow];
    if (iX >= sifted_lp.num_col_) continue;
     HighsInt iCol = sifted_list[iX];
      info_.workValue_[iCol] = sifted_ekk_instance.info_.baseValue_[iRow];
      info_.workDual_[iCol] = 0;
      printf("Basic:    iX %2d: iCol = %3d: Value = %11.4g\n",
	   (int)iX, (int)iCol, info_.workValue_[iCol]);
  }
  for (HighsInt iRow = 0; iRow < lp_.num_row_; iRow++) {
    info_.workValue_[lp_.num_col_+iRow] =
      sifted_ekk_instance.info_.workValue_[sifted_lp.num_col_+iRow];
    info_.workDual_[lp_.num_col_+iRow] =
      sifted_ekk_instance.info_.workDual_[sifted_lp.num_col_+iRow];
  }
  info_.num_primal_infeasibilities = sifted_ekk_instance.info_.num_primal_infeasibilities;
}
