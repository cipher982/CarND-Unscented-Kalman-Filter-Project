#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

  /// create and initialize with 0s
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;


  // Thanks to dariemnt for ideas on good/consistent error messages to use here / future!
  if(estimations.size() == 0){
    cout << "ERROR - CalculateRMSE () - The estimations vector is empty" << endl;
    return rmse;
  }

  if(ground_truth.size() == 0){
    cout << "ERROR - CalculateRMSE () - The ground-truth vector is empty" << endl;
    return rmse;
  }

  if(estimations.size() != ground_truth.size()){
    cout << "ERROR - CalculateRMSE () - The ground-truth and estimations vectors must have the same size." << endl;
    return rmse;
  }

  for (unsigned int i=0; i < estimations.size(); ++i) {
    VectorXd diff = estimations[i] - ground_truth[i]; // residual diffs
    diff = diff.array()*diff.array();
    rmse += diff; // increment
  }

  rmse = rmse / estimations.size(); // take the mean of values
  rmse = rmse.array().sqrt();
  return rmse;
}

