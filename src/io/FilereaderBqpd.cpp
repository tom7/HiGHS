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
/**@file io/FilereaderBqpd.cpp
 * @brief
 */

#include "io/FilereaderBqpd.h"

#include "lp_data/HighsLpUtils.h"


FilereaderRetcode FilereaderBqpd::readModelFromFile(const HighsOptions& options,
                                                  const std::string filename,
                                                  HighsModel& model) {
  return FilereaderRetcode::kParserError;
}

HighsStatus FilereaderBqpd::writeModelToFile(const HighsOptions& options,
                                           const std::string filename,
                                           const HighsModel& model) {
  HighsLp rowwise = model.lp_;
  setFormat(rowwise, MatrixFormat::kRowwise);

  FILE* file = fopen(filename.c_str(), "w");

  // number of variables and number of constraints
  fprintf(file, "%d %d\n", model.lp_.num_col_, model.lp_.num_row_);

  // number of nonzeroes in c and A (write full vector)
  fprintf(file, "%d\n", model.lp_.num_col_ + model.lp_.a_start_[model.lp_.num_col_]);

  // number of nonzeroes in upper triangle of Q
  if (model.hessian_.format_ == HessianFormat::kTriangular)
    fprintf(file, "%d\n", model.hessian_.q_start_[model.lp_.num_col_]);
  else 
    return HighsStatus::kError;

  // estimate of x -> 0
  for (HighsInt i=0; i<model.lp_.num_col_; i++)
    fprintf(file, "%lf ", 0.0);
  fprintf(file, "\n");

  // number of nonzeroes of each [c, constraint]
  fprintf(file, "%d ", model.lp_.num_col_);
  for (HighsInt i=0; i<model.lp_.num_row_; i++)
    fprintf(file, "%d ", rowwise.a_start_[i+1] - rowwise.a_start_[i]); // TODO
  fprintf(file, "\n");

  // c
  for (HighsInt i=0; i<model.lp_.num_col_; i++) {
    fprintf(file, "%d %lf ", i+1, model.lp_.col_cost_[i]);
  }
  fprintf(file, "\n");

  // all rows, each value preceded by its index
  for (HighsInt row=0; row<model.lp_.num_row_; row++) {
    for (HighsInt i=rowwise.a_start_[row]; i<rowwise.a_start_[row+1]; i++)
      fprintf(file, "%d %lf ", rowwise.a_index_[i]+1, rowwise.a_value_[i]);
    fprintf(file, "\n");
  }

  // variable lower bounds, constraint lower bounds
  for (HighsInt i=0; i<model.lp_.num_col_; i++) {
    fprintf(file, "%lf ", model.lp_.col_lower_[i]);
  }
  for (HighsInt i=0; i<model.lp_.num_row_; i++) {
    fprintf(file, "%lf ", model.lp_.row_lower_[i]);
  }
  fprintf(file, "\n");

  // variable upper bounds, constraint upper bounds
  for (HighsInt i=0; i<model.lp_.num_col_; i++) {
    fprintf(file, "%lf ", model.lp_.col_upper_[i]);
  }
  for (HighsInt i=0; i<model.lp_.num_row_; i++) {
    fprintf(file, "%lf ", model.lp_.row_upper_[i]);
  }
  fprintf(file, "\n");

  // upper triangle of Q, each value preceded by its row and column index
  if (model.hessian_.format_ == HessianFormat::kTriangular) {
    for (HighsInt col=0; col<model.hessian_.dim_; col++) {
      for (HighsInt i=model.hessian_.q_start_[col]; i<model.hessian_.q_start_[col+1]; i++) {
        fprintf(file, "%d %d %lf ", col+1, model.hessian_.q_index_[i]+1, model.hessian_.q_value_[i]);
      }
    }
    fprintf(file, "\n");
  } else {
    return HighsStatus::kError;
  }

  // mode of operation (0: solve from scratch)
  fprintf(file, "0\n");

  fclose(file);

  
  return HighsStatus::kOk;
}
