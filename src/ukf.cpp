#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;


/// Apparently I can define numbers like this??
#define EPS  0.001;


/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {

  /// toggle for whether the state has been initialized yet
  is_initialized = false;

  /// if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  /// if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  /// State dimension
  n_x_ = 5;

  /// initial state vector
  x_ = VectorXd(n_x_);

  /// initial state covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  P_ <<   1, 0, 0, 0, 0, // Identity Matrix
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;

  /// Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1;

  /// Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;

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

  /// State vector dimension
  n_x_ = 5;

  ///  Augmented state dimension
  n_aug_ = n_x_ + 2;

  ///  Sigma point spreading parameter
  lambda_ = 3 - n_x_;

  /// number of sigma points
  n_sig_ = 2 * n_aug_ + 1; // 2n + 1

  /// Weights of sigma points
  weights_ = VectorXd(n_sig_); // 2n + 1 sigma points
  weights_.fill(0.5 / (n_aug_ + lambda_));
  weights_(0) = lambda_ / (lambda_ + n_aug_);


  R_laser_ = MatrixXd(2,2);
  R_laser_ << std_laspx_ * std_laspx_, 0, 0, std_laspy_ * std_laspy_; // diagonals

  R_radar_ = MatrixXd(3,3);
  R_radar_ << std_radr_ * std_radr_, 0, 0, // rho
          0, std_radphi_*std_radphi_, 0, // phi
          0, 0,std_radrd_*std_radrd_; // rho_dot


}

/**
 *
 * @param phi The input angle from radar measurement
 * @return
 */
/// hanging my computer
/*
static long NormalizePhiAngle (long phi) {
  while (phi >= M_PI) phi -= 2. * M_PI;
  while (phi < M_PI) phi += 2. * M_PI;

  return phi;
}
*/

/// new method from forum mentor driveWell
static double NormalizePhiAngle(double phi) {
  if (phi > M_PI) {
    double temp = fmod((phi - M_PI), (2 * M_PI)); // -= 2. * M_PI;
    phi = temp - M_PI;
  } // phi normalization
  if (phi < -M_PI) {
    double temp = fmod((phi + M_PI) ,(2 * M_PI));
    phi = temp + M_PI;
  }
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



    if (is_initialized == false)
    {

      /*****************************************************************************
      *  Initialization
      ****************************************************************************/

      /// first measurement
      // x_ << 1, 1, 1, 1, 0;

      /// print this out for debugging:
      cout << "UKF: " << endl;

      x_ << 0, 0, 0, 0, 0;
      cout << "did x_" << x_;

      /// added for readability, place into x_ later near end
      double px;
      double py;
      double vx;
      double vy;


      /// initialize covariance matrix
      P_ <<   1, 0, 0, 0, 0,
              0, 1, 0, 0, 0,
              0, 0, 1, 0, 0,
              0, 0, 0, 1, 0,
              0, 0, 0, 0, 1;

      /// laser measurements
      if (meas_package.sensor_type_ == MeasurementPackage::LASER) {

        cout << "initialize laser" << endl;

        /// initialize state
        px = meas_package.raw_measurements_(0);
        py = meas_package.raw_measurements_(1);
        vx = 0;
        vy = 0;

        cout << "got measurements! PX = " << px << "PY = " << py << endl;

        /// convert back to x_
        x_ << px, py, 0, 0, 0; // laser has no velocity information
        cout << "converted measurements! " << endl;


      } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {

        cout << "initialize radar" << endl;


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
        x_ << px, py, sqrt((vx * vx) + (vy * vy)), 0, 0; // radar has velocity information
      }

      cout << "Initialized x_: " << x_ << endl;


      /// grab timestamp in ultra-low seconds (us)
      time_us_ = meas_package.timestamp_ ; // Without this it stops after 1 step?? oops!

      /// finished initialization;
      is_initialized = true;
      cout << "initialized is now true" << endl;

      return;
    }
  }

  /// Initialize the timestamp
  cout << "Initialize the timestamp" << endl;
  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;
  Prediction(dt);


  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_ == true) {
    cout << "Update Radar line 223" << endl;
    UpdateRadar(meas_package);
  }
  if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_ == true) {
    cout << "Update laser line 227" << endl;
    UpdateLidar(meas_package);
  }

}



/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
  void UKF::Prediction(double delta_t) {

  /*****************************************************************************
  *  PREDICTION
  ****************************************************************************/

  /// estimate object location??

  cout << "starting prediction" << endl;

  /// SIGMA POINTS ///

  ///  create augmented mean vector
  cout << "n_aug_ is " << n_aug_ << endl;
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;
  cout << "Created augmented mean vector" << endl;

  /// create augmented state covariance matrix
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_); // augmented dimensions
  cout << "1" << endl;
  P_aug.fill(0.0); // fill with floating 0
  cout << "2" << endl;
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  cout << "3" << endl;
  P_aug(5, 5) = std_a_ * std_a_; // proc noise std dev long accel
  cout << "4" << endl;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;
  cout << "5" << endl;
  cout << "Created augmented state vector" << endl;

  /// create sigma points
  MatrixXd X_sig_aug = GenerateSigmaPoints(x_aug, P_aug, lambda_, n_sig_);
  cout << "Created sigma points" << endl;

  /// predict sigma points
  Xsig_pred = PredictSigmaPoints(X_sig_aug, delta_t, n_x_, n_sig_,std_a_, std_yawdd_);
  cout << "Predicted sigma points into Xsig_pred: "<< endl << Xsig_pred << endl;

  /// predict state mean and covariance ///

  cout << "now try to mult Xsig_pred and weights_ " << endl;
  cout << Xsig_pred * weights_ << endl;

  /// state mean
  x_ = Xsig_pred * weights_;
  cout << "Got state mean x_ " << x_ << endl;

  /// initialize and fill prediction state covariance matrix
  P_.fill(0.0);

  cout << "About to iterage sigma points in Prediction" << endl;
  for (int i = 0; i < n_sig_; i++) // iterate sigma points
  {
    /// calculate the state difference each column
    cout << "Xsig_pred is :" << Xsig_pred << endl;
    VectorXd x_diff = Xsig_pred.col(i) - x_; // create vector to use below

    cout << "Before: " << x_diff(3) << endl;
    /// normalize the angle
    x_diff(3) = NormalizePhiAngle(x_diff(3)); // index 3 = angle to normalize
    cout << "After: " << x_diff(3) << endl;

    /// predicted covariance formula from lecture
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose(); // sum up diffs.t

    cout << "iterated sigma prediction point " << i << endl;

  }

  cout << "Calculated x_diff and incremented P_" << endl;
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

  cout << "Start laser update" << endl;

  /// set measurement prediction Matrix
  int n_z = 2; // px and py
  MatrixXd Zsig = Xsig_pred.block(0, 0, n_z, n_sig_);
  cout << "Set the measurement pred matrix" << endl;

  /// initialize mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }
  cout << "initialized the mean prediction measurement" << endl;

  /// initialize measurement covariance Matrix
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {

    /// calculate residuals
    VectorXd z_diff = Zsig.col(i) - z_pred;

    /// covariance formula again
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  cout << "Set the measurement covariance Matrix (laser)" << endl;

  /// add measurement noise covariance Matrix (Tc)
  S = S + R_laser_;

  /// update the state now!! ///

  /// laser measurement
  VectorXd z = meas_package.raw_measurements_; // pull in data from measurements

  /// create matrix for cross validation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  cout << "grabbed laser measurement and cross val Mat Tc (laser)" << endl;

  Tc.fill(0.0); // populate blank template
  for (long i = 0; i < n_sig_; i++) { // 2n + 1 sigma points

    /// residuals
    VectorXd z_diff = Zsig.col(i) - z_pred;

    /// state difference
    VectorXd x_diff = Xsig_pred.col(i) - x_;

    /// add in the noise!!
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  cout << "calculated residual/diff/noise in the loop" << endl;

  /// Kalman gain
  MatrixXd K = Tc * S.inverse();

  /// residual diff
  VectorXd z_diff = z - z_pred;

  /// update state mean and covariance matrix
  x_ = x_ + (K * z_diff); // x_ is mean
  P_ = P_ - (K * S * K.transpose());

  /// NIS laser
  NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;

  cout << "Finished update (laser)" << endl;

}

/*****************************************************************************
 *  RADAR
 ****************************************************************************/

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
  void UKF::UpdateRadar(MeasurementPackage meas_package) {

  ///////////////////////////////////
  /// *** PREDICT MEASUREMENT *** ///
  ///////////////////////////////////

  cout << "Start radar update" << endl;

  /// Radar dimensions
  int n_z = 3; // rho, phi, rho_dot

  /// measurement prediction matrix
  MatrixXd Z_sig = MatrixXd(n_z, n_sig_);
  cout << "Finished X_sig Matrix" << endl;

  /// transform sigma points into measurement space
  for (int i = 0; i < n_sig_; i++) { // 2n+1 sigma points

    cout << "going to start the named variables now" << endl;
    /// extract to named variables
    cout << "going to try PX: " << Xsig_pred << endl;
    double px  = Xsig_pred(0, i);
    cout << "got px" << endl;
    double py  = Xsig_pred(1, i);
    double v   = Xsig_pred(2, i);
    double yaw = Xsig_pred(3, i);

    double vx = cos(yaw) * v;
    double vy = sin(yaw) * v;

    cout << "got the named variables in" << endl;
    /// measurement model
    Z_sig(0, i) = sqrt(px * px + py * py); // rho
    Z_sig(1, i) = atan2(py, px); // phi
    Z_sig(2, i) = ((px * vx) + (py * vy)) / sqrt((px * px) + (py * py));  // rho_dot

  }

  cout << "Finished sigma to measurement space" << endl;

  /// mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {
    z_pred = z_pred + (weights_(i) * Z_sig.col(i));
  }

  /// measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {  // 2n + 1 sigma points

    /// residual
    VectorXd z_diff = Z_sig.col(i) - z_pred;

    /// normalize angle in case of poor values
    z_diff(1) = NormalizePhiAngle(z_diff(1)); // 2nd value = phi

    /// add up the increments
    S = S + (weights_(i) * z_diff * z_diff.transpose());

  }


  /// add measurement noise covariance matrix
  S = S + R_radar_;

  ///////////////////////////////////
  /// ****** UPDATE STATE  ****** ///
  ///////////////////////////////////

  /// Incoming radar measurement
  VectorXd z = meas_package.raw_measurements_;

  /// create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);

  for (int i = 0; i < n_sig_; i++) {  // 2n + 1 sigma points

    /// residual
    VectorXd z_diff = Z_sig.col(i) - z_pred;

    /// angle normalization
    z_diff(1) = NormalizePhiAngle(z_diff(1));

    /// state difference
    VectorXd x_diff = Xsig_pred.col(i) - x_;

    /// normalize angle in case of poor values
    x_diff(3) = NormalizePhiAngle(x_diff(3));

    /// add in the noise!
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();

  }

  /// Kalman gain
  MatrixXd K = Tc * S.inverse();

  /// residual
  VectorXd z_diff = z - z_pred;

  /// normalize angle in case of poor values

  z_diff(1) = NormalizePhiAngle(z_diff(1));

  /// update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

  /// NIS Radar
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;

}


/**
 * Predits sigma points.
 * @param Xsig : Sigma points to predict.
 * @param delta_t : Time between k and k+1 in s
 * @param n_x : State dimension.
 * @param n_sig : Sigma points dimension.
 * @param nu_am : Process noise standard deviation longitudinal acceleration in m/s^2
 * @param nu_yawdd : Process noise standard deviation yaw acceleration in rad/s^2
 */
MatrixXd UKF::PredictSigmaPoints(MatrixXd Xsig, double delta_t, int n_x, int n_sig, double nu_am, double nu_yawdd) {
  MatrixXd Xsig_pred = MatrixXd(n_x, n_sig);

  /// predict sigma points
  cout << "Beginning to predict sigma points" << endl;
  for (int i = 0; i< n_sig; i++)
  {
    /// extract values for better readability
    double p_x      = Xsig(0,i);
    double p_y      = Xsig(1,i);
    double v        = Xsig(2,i);
    double yaw      = Xsig(3,i);
    double yawd     = Xsig(4,i);
    double nu_a     = Xsig(5,i);
    double nu_yawdd = Xsig(6,i);

    /// predicted state values
    double px_p, py_p;


    /// avoid division by zero
    if (fabs(yawd) > 0.001) {
      px_p = p_x + v / yawd * ( sin (yaw + yawd * delta_t) - sin(yaw));
      py_p = p_y + v / yawd * ( cos (yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
      px_p = p_x + v*delta_t*cos(yaw);
      py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p    = v;
    double yaw_p  = yaw + yawd*delta_t;
    double yawd_p = yawd;

    /// add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p  = v_p  + nu_a*delta_t;

    yaw_p  = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    /// write predicted sigma point into right column
    Xsig_pred(0,i) = px_p;
    Xsig_pred(1,i) = py_p;
    Xsig_pred(2,i) = v_p;
    Xsig_pred(3,i) = yaw_p;
    Xsig_pred(4,i) = yawd_p;
  }

  cout << "Made the sigma points loop" << endl;


  return Xsig_pred;
}

/**
 *   Generate sigma points:
 *  @param x : State vector.
 *  @param P : Covariance matrix.
 *  @param lambda: Sigma points spreading parameter.
 *  @param n_sig: Sigma points dimension.
 */
MatrixXd UKF::GenerateSigmaPoints(VectorXd x, MatrixXd P, double lambda, int n_sig) {

  cout << "Begin to generate sigma points" << endl;

  /// create sigma point matrix
  MatrixXd Xsig = MatrixXd( x.size(), n_sig );

  /// calculate square root of P
  MatrixXd A = P.llt().matrixL();

  Xsig.col(0) = x;

  double lambda_plue_n_x_sqrt = sqrt(lambda + x.size());
  for (int i = 0; i < x.size(); i++){
    Xsig.col( i + 1 ) = x + lambda_plue_n_x_sqrt * A.col(i);
    Xsig.col( i + 1 + x.size() ) = x - lambda_plue_n_x_sqrt * A.col(i);
  }
  return Xsig;
}



