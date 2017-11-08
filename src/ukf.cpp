#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  /// if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  /// if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  /// initial state vector
  x_ = VectorXd(5);

  /// initial state covariance matrix
  P_ = MatrixXd(5, 5);

  /// Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  /// Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

  /// Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  /// Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  /// Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  /// Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  /// Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /// Weights of sigma points
  weights_ = VectorXd(2 * n_aug_ + 1); // 2n + 1 sigma points
  weights_.fill(0.5 / (n_aug_ + lambda_));
  weights_(0) = lambda_ / (lambda_ + n_aug_);

  /// State dimension
  n_x_ = 5;

  ///  Augmented state dimension
  n_aug_ = 7;

  ///  Sigma point spreading parameter
  lambda_ = 2;

  /// number of sigma points
  n_sig_ = 2 * n_aug_ + 1; // 2n + 1

  R_laser = MatrixXd(2,2);
  R_laser << std_laspx_ * std_laspx_, 0, 0, std_laspy_ * std_laspy_; // diagonals

  R_radar = MatrixXd(3,3);
  R_radar << std_radr_ * std_radr_, 0, 0, // rho
          0, std_radphi_*std_radphi_, 0, // phi
          0, 0,std_radrd_*std_radrd_; // rho_dot


}

/**
 *
 * @param phi The input angle from radar measurement
 * @return
 */
static double NormalizePhiAngle (double phi) {
  while (phi >= M_PI) phi -= 2. * M_PI;
  while (phi < M_PI) phi += 2. * M_PI;

  return phi;
}


UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  if ((meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_ == true) ||
      (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_ == true)) {

    if (is_initialized == false);
    {

      /*****************************************************************************
      *  Initialization
      ****************************************************************************/

      /// first measurement
      // x_ << 1, 1, 1, 1, 0;

      /// added for readability, place into x_ later near end
      double px;
      double py;
      double vx;
      double vy;


      /// initialize covariance matrix
      P_ << 1, 0, 0, 0, 0,
              0, 1, 0, 0, 0,
              0, 0, 1, 0, 0,
              0, 0, 0, 1, 0,
              0, 0, 0, 0, 1;

      /// Initialize the timestamp
      time_us_ = meas_package.timestamp_;

      /// laser measurements
      if (meas_package.sensor_type_ == MeasurementPackage::LASER) {

        /// initialize state
        px = meas_package.raw_measurements_(0);
        py = meas_package.raw_measurements_(1);

        /// convert back to x_
        x_ << px, py, 0, 0, 0; // laser has no velocity information

      } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {

        /// initialize state
        float rho = meas_package.raw_measurements_(0);
        float phi = meas_package.raw_measurements_(1);
        float rho_dot = meas_package.raw_measurements_(2);

        /// convert from polar to cartesian
        px = rho * cos(phi);
        py = rho * sin(phi);

        vx = rho_dot * cos(phi);
        vy = rho_dot * sin(phi);

        /// convert back to x_
        x_ << px, py, sqrt(vx * vx + vy * vy), 0, 0; // radar has velocity information
      }

      cout << "Initialized x_: " << x_ << endl;


      /// finished initialization;
      is_initialized_ = true;

      return;
    }
  }

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_ == true) {
    UpdateRadar(meas_package);
  }
  if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_ == true) {
    UpdateLidar(meas_package);
  }

}



/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
  void UKF::Prediction(double delta_t) {

  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  /// estimate object location??


  /// SIGMA POINTS ///

  ///  create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  /// create augmented state covariance matrix
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_); // augmented dimensions
  P_aug.fill(0.0); // fill with floating 0
  P_aug.topLeftCorner(n_aug_, n_aug_);
  P_aug(5, 5) = std_a_ * std_a_; // proc noise std dev long accel
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;

  /// create sigma points
  MatrixXd X_sig_aug = GenerateSigmaPoints(x_aug, P_aug, lambda_
  n_sig_);

  /// predict sigma points
  Xsig_pred_ = PredictSigmaPoints(Xsig_aug, delta_t, n_x_, n_sig_, std_yawdd_);

  /// predict state mean and covariance ///

  /// state mean
  x_ = Xsig_pred_ * weights_;

  /// initialize and fill prediction state covariance matrix
  P_.fill(0.0);

  for (int i = 0; i < n_sig_; i++) // iterate sigma points
  {
    /// calculate the state difference each column
    VectorXd x_diff = Xsig_pred_.col(i) - x_; // create vector to use below

    /// normalize the angle
    x_diff(3) = NormalizePhiAngle(x_diff(3)); // index 3 = angle to normalize

    /// predicted covariance formula from lecture
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose(); // sum up diffs.t
  }
}

/*****************************************************************************
 *  LASER
 ****************************************************************************/

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
  void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  /// set measurement prediction Matrix
  int n_z = 2; // px and py
  MatrixXd Zsig = Xsig_pred_.block(0, 0, n_z, n_sig_);

  /// initialize mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  /// initialize measurement covariance Matrix
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {

    /// calculate residuals
    VectorXd z_diff = Zsig.col(i) - z_pred;

    /// covariance formula again
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  /// add measurement noise covariance Matrix (Tc)
  S = S + R_laser;

  /// update the state now!! ///

  /// laser measurement
  VectorXd z = meas_package.raw_measurements_; // pull in data from measurements

  /// create matrix for cross validation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  Tc.fill(0.0); // populate blank template
  for (long i = 0; i < n_sig_; i++) { // 2n + 1 sigma points

    /// residuals
    VectorXd z_diff = Zsig.col(i) - z_pred;

    /// state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    /// add in the noise!!
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();

  }

  /// Kalman gain
  MatrixXd K = Tc * S.inverse();

  /// residual diff
  VectorXd z_diff = z - z_pred;

  /// update state mean and covariance matrix
  x_ = x_ + (K * z_diff); // x_ is mean
  P_ = P_ - (K * S * K.transpose());

  /// NIS laser
  NIS_laser = z_diff.transpose() * S.inverse() * z_diff;

}

/*****************************************************************************
 *  RADAR
 ****************************************************************************/

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
  void UKF::UpdateRadar(MeasurementPackage meas_package) {


  /// Radar dimensions
  int n_z = 3; // rho, phi, rho_dot

  /// measurement prediction matrix
  MatrixXd Zsig = MatrixXd(n_z, n_sig_);

  /// transform sigma points into measurement space
  for (int i = 0; i < n_sig_; i++) { // 2n+1 sigma points

    /// extract to named variables
    double px  = Xsig_pred_(0, i);
    double py  = Xsig_pred_(1, i);
    double v   = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    double vx = cos(yaw) * v;
    double vy = sin(yaw) * v;

    /// measurement model
    Zsig(0, i) = sqrt(px*px + py*py); // rho
    Zsig(1, i) = atan2(py, px); // phi
    Zsig(2, i) = (px * vx + py * vy) / sqrt(px * px + py * py);  // rho_dot



  }




  }
