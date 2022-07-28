
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
  
  if (info_.num_primal_infeasibilities > 0) {
    // Construct an LP containing any basic columns and a random
    // collection of nonbasic columns
    for (HighsInt iCol = 0; iCol < lp_.num_col_; iCol++) {
      if (in_sifted_list[iCol]) continue;
      if (basis_.nonbasicFlag_[iCol] == 0) {
	sifted_list.push_back(iCol);
	in_sifted_list[iCol] = true;
      }
    }
    HighsRandom random;
    for (;;) {
      if (sifted_list.size() == sifted_list_max_count) break;
      HighsInt iCol = random.integer(lp_.num_col_);
      if (in_sifted_list[iCol]) continue;
      sifted_list.push_back(iCol);
      in_sifted_list[iCol] = true;
    }
    
    getSiftedSolverObject(sifted_solver_object, in_sifted_list);
    sifted_ekk_instance.moveLp(sifted_solver_object);

     return_status = sifted_ekk_instance.solve();
 
  }

  assert(1 == 0);
  getNewSiftedList(sifted_list_max_count, sifted_list, in_sifted_list);
  return return_status;
}

void HEkk::getNewSiftedList(const HighsInt new_sifted_list_max_count,
			    std::vector<HighsInt>& new_sifted_list,
			    const std::vector<bool>& in_sifted_list) {
  assert(new_sifted_list.size() == 0);
  HighsInt new_sifted_list_count = 0;
  std::vector<HighsInt> heap_index;
  std::vector<double> heap_value;
  heap_index.push_back(0);
  heap_value.push_back(0);
  for (HighsInt iCol = 0; iCol < lp_.num_col_; iCol++) {
    if (in_sifted_list[iCol]) continue;
    if (basis_.nonbasicFlag_[iCol] == 0) {
      // Basic, so in sifted list
      new_sifted_list.push_back(iCol);
      continue;
    }
    const double dual = info_.workDual_[iCol];
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
  new_sifted_list_count = new_sifted_list.size();
  HighsInt heap_num_en = heap_index.size() - 1;
  assert(heap_num_en >= 0);
  if (heap_num_en > 0) {
    // There are dual infeasibilities
    maxheapsort(&heap_value[0], &heap_index[0], heap_num_en);
    for (HighsInt iEl = 1; iEl <= heap_num_en; iEl++) {
      if (new_sifted_list_count == new_sifted_list_max_count) break;
      HighsInt iCol = heap_index[iEl];
      new_sifted_list.push_back(iCol);
    }
  }
}

void HEkk::getSiftedSolverObject(HighsLpSolverObject& sifted_solver_object, 
				 const std::vector<bool>& in_sifted_list) {
  HighsLp& sifted_lp = sifted_solver_object.lp_;
  HighsBasis& sifted_basis = sifted_solver_object.basis_;
  std::vector<double> unsifted_row_activity;
  unsifted_row_activity.assign(lp_.num_row_, 0);
  sifted_lp.a_matrix_.start_.clear();
  for (HighsInt iCol = 0; iCol < lp_.num_col_; iCol++) {
    if (in_sifted_list[iCol]) {
      sifted_lp.a_matrix_.start_.push_back(sifted_lp.a_matrix_.index_.size());
      sifted_lp.col_cost_.push_back(info_.workCost_[iCol]);
      sifted_lp.col_lower_.push_back(info_.workLower_[iCol]);
      sifted_lp.col_upper_.push_back(info_.workUpper_[iCol]);
      for (HighsInt iEl = lp_.a_matrix_.start_[iCol]; iEl < lp_.a_matrix_.start_[iCol+1]; iEl++) {
	sifted_lp.a_matrix_.index_.push_back(lp_.a_matrix_.index_[iEl]);
	sifted_lp.a_matrix_.value_.push_back(lp_.a_matrix_.value_[iEl]);
      }
    } else {
      assert(basis_.nonbasicFlag_[iCol]);
      double value = info_.workValue_[iCol];
      for (HighsInt iEl = lp_.a_matrix_.start_[iCol]; iEl < lp_.a_matrix_.start_[iCol+1]; iEl++) {
	HighsInt iRow = lp_.a_matrix_.index_[iEl];
	unsifted_row_activity[iRow] += value * lp_.a_matrix_.value_[iEl];
      }
    }
  }
  sifted_lp.a_matrix_.start_.push_back(sifted_lp.a_matrix_.index_.size());
  for (HighsInt iRow = 0; iRow < lp_.num_row_; iRow++) {
    sifted_lp.row_lower_.push_back(-info_.workUpper_[iRow]-unsifted_row_activity[iRow]);
    sifted_lp.row_upper_.push_back(-info_.workLower_[iRow]-unsifted_row_activity[iRow]);
    printf("Row %3d has bounds [%11.4g, %11.4g] from  [%11.4g, %11.4g] and shift %11.4g\n",
	   (int)iRow, sifted_lp.row_lower_[iRow], sifted_lp.row_upper_[iRow],
	   info_.workLower_[iRow], info_.workUpper_[iRow], unsifted_row_activity[iRow]);
  }
}
