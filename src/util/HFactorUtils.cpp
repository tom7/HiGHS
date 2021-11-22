/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Qi Huangfu, Leona Gottwald    */
/*    and Michael Feldmeier                                              */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HFactorUtils.cpp
 * @brief Types of solution classes
 */
#include "util/HFactor.h"

void HFactor::invalidAMatrixAction() {
  this->a_matrix_valid = false;
  refactor_info_.clear();
}

void HFactor::reportLu(const HighsInt l_u_or_both, const bool full) const {
  if (l_u_or_both < kReportLuJustL || l_u_or_both > kReportLuBoth) return;
  if (l_u_or_both & 1) {
    printf("L");
    if (full) printf(" - full");
    printf(":\n");

    if (full) reportIntVector("l_pivot_lookup", l_pivot_lookup);
    if (full) reportIntVector("l_pivot_index", l_pivot_index);
    reportIntVector("l_start", l_start);
    reportIntVector("l_index", l_index);
    reportDoubleVector("l_value", l_value);
    if (full) {
      reportIntVector("lr_start", lr_start);
      reportIntVector("lr_index", lr_index);
      reportDoubleVector("lr_value", lr_value);
    }
  }
  if (l_u_or_both & 2) {
    printf("U");
    if (full) printf(" - full");
    printf(":\n");
    if (full) reportIntVector("u_pivot_lookup", u_pivot_lookup);
    reportIntVector("u_pivot_index", u_pivot_index);
    reportDoubleVector("u_pivot_value", u_pivot_value);
    reportIntVector("u_start", u_start);
    if (full) reportIntVector("u_last_p", u_last_p);
    reportIntVector("u_index", u_index);
    reportDoubleVector("u_value", u_value);
    if (full) {
      reportIntVector("ur_start", ur_start);
      reportIntVector("ur_lastp", ur_lastp);
      reportIntVector("ur_space", ur_space);
      for (HighsInt iRow = 0; iRow < ur_start.size(); iRow++) {
        const HighsInt start = ur_start[iRow];
        const HighsInt end = ur_lastp[iRow];
        if (start >= end) continue;
        printf("UR    Row %2d: ", (int)iRow);
        for (HighsInt iEl = start; iEl < end; iEl++)
          printf("%11d ", (int)ur_index[iEl]);
        printf("\n              ");
        for (HighsInt iEl = start; iEl < end; iEl++)
          printf("%11.4g ", ur_value[iEl]);
        printf("\n");
      }
      //      reportIntVector("ur_index", ur_index);
      //      reportDoubleVector("ur_value", ur_value);
    }
  }
  if (l_u_or_both == 3 && full) {
    reportDoubleVector("pf_pivot_value", pf_pivot_value);
    reportIntVector("pf_pivot_index", pf_pivot_index);
    reportIntVector("pf_start", pf_start);
    reportIntVector("pf_index", pf_index);
    reportDoubleVector("pf_value", pf_value);
  }
}

void HFactor::reportIntVector(const std::string name,
                              const vector<HighsInt> entry) const {
  const HighsInt num_en = entry.size();
  printf("%-12s: siz %4d; cap %4d: ", name.c_str(), (int)num_en,
         (int)entry.capacity());
  for (HighsInt iEn = 0; iEn < num_en; iEn++) {
    if (iEn > 0 && iEn % 10 == 0)
      printf("\n                                  ");
    printf("%11d ", (int)entry[iEn]);
  }
  printf("\n");
}
void HFactor::reportDoubleVector(const std::string name,
                                 const vector<double> entry) const {
  const HighsInt num_en = entry.size();
  printf("%-12s: siz %4d; cap %4d: ", name.c_str(), (int)num_en,
         (int)entry.capacity());
  for (HighsInt iEn = 0; iEn < num_en; iEn++) {
    if (iEn > 0 && iEn % 10 == 0)
      printf("\n                                  ");
    printf("%11.4g ", entry[iEn]);
  }
  printf("\n");
}

void HFactor::analyseActiveKernel(const std::string message, const bool report) {
  if (report)
    highsLogDev(log_options_, HighsLogType::kInfo, "\n%s\n", message.c_str());
  initialiseValueDistribution("Kernel value", "", 1e-20, 1e20, 10.0, analyse_kernel_value_);
  initialiseValueDistribution("Kernel row counts", "", 1, 16384, 2, analyse_kernel_row_count_);
  initialiseValueDistribution("Kernel col counts", "", 1, 16384, 2, analyse_kernel_col_count_);
  for (HighsInt count = 0; count <= num_row; count++) {
    for (HighsInt j = col_link_first[count]; j != -1; j = col_link_next[j]) {
      updateValueDistribution(count, analyse_kernel_col_count_);
      HighsInt start = mc_start[j];
      HighsInt end = start + mc_count_a[j];
      for (HighsInt k = start; k < end; k++)
        updateValueDistribution(mc_value[k], analyse_kernel_value_);
    }
  }
  for (HighsInt count = 0; count <= num_basic; count++) {
    for (HighsInt j = row_link_first[count]; j != -1; j = row_link_next[j])
      updateValueDistribution(count, analyse_kernel_row_count_);
  }
  if (report) {
    logValueDistribution(log_options_, analyse_kernel_value_);
    logValueDistribution(log_options_, analyse_kernel_row_count_);
    logValueDistribution(log_options_, analyse_kernel_col_count_);
  }
}

void HFactor::reportKernelValueChange(const std::string message,
                                      const HighsInt iRow, const HighsInt iCol,
                                      double& track_value) {
  if (iRow < 0 || iCol < 0) return;
  double latest_value = 0;
  HighsInt start = mc_start[iCol];
  HighsInt end = start + mc_count_a[iCol];
  for (HighsInt k = start; k < end; k++) {
    if (mc_index[k] == iRow) {
      latest_value = mc_value[k];
      break;
    }
  }
  if (track_value != latest_value)
    printf("\n%s: Change of %11.4g in B(%6d, %6d) to %11.4g\n", message.c_str(),
           latest_value - track_value, (int)iRow, (int)iCol, latest_value);
  track_value = latest_value;
}

void HFactor::reportMcColumn(const HighsInt num_pivot,
                             const HighsInt iCol) const {
  if (iCol >= num_basic) return;
  const HighsInt iCol_start = mc_start[iCol];
  const HighsInt iCol_count_a = mc_count_a[iCol];
  HighsInt next_col = -1;
  HighsInt next_start = kHighsIInf;
  const HighsInt num_mc_col = mc_start.size();
  assert(num_mc_col == num_basic);
  for (HighsInt i = 0; i < num_basic; i++) {
    if (mc_start[i] > iCol_start && mc_start[i] < next_start) {
      next_col = i;
      next_start = mc_start[i];
    }
  }
  printf(
      "McColumn %d: (Pivot %6d) Var %6d; Start %8d; Count(A %3d; N %3d); Space "
      "%6d: Next col is %6d; Start %8d",
      (int)iCol, (int)num_pivot, (int)mc_var[iCol], (int)mc_start[iCol],
      (int)mc_count_a[iCol], (int)mc_count_n[iCol], (int)mc_space[iCol],
      (int)next_col, (int)next_start);
  HighsInt en = 0;
  for (HighsInt iEl = mc_start[iCol]; iEl < mc_start[iCol] + mc_count_a[iCol];
       iEl++) {
    if (en % 5 == 0) printf("\n");
    en++;
    printf("[%6d %11.4g] ", (int)mc_index[iEl], mc_value[iEl]);
  }
  printf("\n");
  fflush(stdout);
  assert(iCol_start + iCol_count_a <= next_start);
}

void HFactor::AnalyseBuild::clear() {
  this->num_row = 0;
  this->num_col = 0;
  this->num_basic = 0;
  this->basic_num_nz = 0;
  this->num_simple_pivot = 0;
  this->num_kernel_pivot = 0;
  this->kernel_initial_num_nz = 0;
  this->kernel_final_num_nz = 0;
  this->invert_num_nz = 0;
  this->sum_merit = 0;
}

HighsInt HFactor::getKernelNumNz() const {
  HighsInt kernel_num_nz = 0;
  for (HighsInt iCol = 0; iCol < num_basic; iCol++) 
    kernel_num_nz += mc_count_a[iCol] + mc_count_n[iCol];
  return kernel_num_nz;
}

void HFactor::debugReportAnalyseBuild(const HighsLogOptions& local_log_options,
				      const std::string message) {
  const HighsInt report_rank_deficiency =
    analyse_build_record_.num_basic < analyse_build_record_.num_row ?
			       rank_deficiency - (analyse_build_record_.num_row-analyse_build_record_.num_basic) : rank_deficiency;
  highsLogDev(local_log_options, HighsLogType::kInfo,
	      "\n%s has %d rows, %d columns, %d nonzeros and rank deficiency %d\n",
	      message.c_str(),
	      (int)analyse_build_record_.num_row,
	      (int)analyse_build_record_.num_basic,
	      (int)analyse_build_record_.basic_num_nz,
	      (int)report_rank_deficiency);
  highsLogDev(local_log_options, HighsLogType::kInfo,
	      "   Number of simple pivots = %d\n", (int)analyse_build_record_.num_simple_pivot);
  if (analyse_build_record_.num_kernel_pivot) {
    assert(analyse_build_record_.kernel_initial_num_nz>0);
    const HighsInt kernel_num_row = analyse_build_record_.num_row - analyse_build_record_.num_simple_pivot;
    const HighsInt kernel_num_col = analyse_build_record_.num_basic - analyse_build_record_.num_simple_pivot;
    highsLogDev(local_log_options, HighsLogType::kInfo,
		"   Kernel has %d rows, %d columns, %d nonzeros\n", (int)kernel_num_row, (int)kernel_num_col,
		analyse_build_record_.kernel_initial_num_nz);
    highsLogDev(local_log_options, HighsLogType::kInfo,
		"   Number of kernel pivots = %d\n", (int)analyse_build_record_.num_kernel_pivot);
    double kernel_fill_factor = (1.0 * analyse_build_record_.kernel_final_num_nz) / analyse_build_record_.kernel_initial_num_nz;
    highsLogDev(local_log_options, HighsLogType::kInfo,
		"   Sum of pivot merit = %g (average %g)\n", analyse_build_record_.sum_merit,
		analyse_build_record_.sum_merit / analyse_build_record_.num_kernel_pivot);
    highsLogDev(local_log_options, HighsLogType::kInfo,
		"   Kernel final number of nonzeros = %d: fill factor of %g\n",
		(int)analyse_build_record_.kernel_final_num_nz, kernel_fill_factor);
  } 
  logValueDistribution(local_log_options, analyse_pivot_col_count_);
  logValueDistribution(local_log_options, analyse_pivot_row_count_);
  logValueDistribution(local_log_options, analyse_pivot_merit_);
  logValueDistribution(local_log_options, analyse_pivot_value_);
}
