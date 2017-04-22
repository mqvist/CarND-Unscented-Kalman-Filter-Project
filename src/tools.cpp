#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  // Check preconditions
  assert(estimations.size());
  assert(estimations.size() == ground_truth.size());
  assert(estimations[0].size() == 4);

  VectorXd rmse = VectorXd::Zero(4);
  // Check for invalid input values
  // Calculate the sum of squared residuals
  for (int i = 0; i < estimations.size(); ++i) {
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array() * residual.array();
    rmse += residual;
  }
  // Calculate mean
  rmse /= estimations.size();
  // Calculate square root of the mean and return
  rmse = rmse.array().sqrt();
  // Check postconditions
  assert(rmse.size() == 4);
  return rmse;
}
