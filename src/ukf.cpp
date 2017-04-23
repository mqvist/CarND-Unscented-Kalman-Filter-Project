#include "ukf.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  is_initialized_ = false;
  previous_timestamp_ = 0;
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_.setIdentity();

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  n_x_ = 5;
  n_aug_ = 2 * n_x_ + 1;
  lambda_ = 3 - n_aug_;
  weights_ = VectorXd(2 * n_aug_ + 1);
  // set weights
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < weights_.size(); ++i) {
    weights_(i) = 1.0 / (2.0 * (lambda_ + n_aug_));
  }
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (!is_initialized_) {
    VectorXd m = meas_package.raw_measurements_;
    switch (meas_package.sensor_type_) {
      case MeasurementPackage::RADAR:
        cout << "Init with radar measurement: " << m << endl;
        x_ << m[0] * cos(m[1]), m[0] * sin(m[1]), 0, 0, 0;
        break;
      case MeasurementPackage::LASER:
        cout << "Init with Lidar measurement: " << m << endl;
        x_ << m[0], m[1], 0, 0, 0;
        break;
      default:
        assert(0);
    }
    cout << "Initial state: " << x_ << endl;
    previous_timestamp_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }
  double delta_t = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;
  Prediction(delta_t);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  } else {
    UpdateLidar(meas_package);
  }

  previous_timestamp_ = meas_package.timestamp_;
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
  vector, x_. Predict sigma points, the state, and the state covariance
  matrix.
  */
  cout << "Prediction, delta_t = " << delta_t << endl;

  // Augmentation
  // create augmented mean state
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.setZero();
  x_aug.head(n_x_) = x_;
  // create augmented covariance matrix
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.setZero();
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;
  // create square root matrix
  MatrixXd A = P_aug.llt().matrixL();
  MatrixXd B = sqrt(lambda_ + n_aug_) * A;
  // create augmented sigma points
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i < n_aug_; ++i) {
    Xsig_aug.col(i + 1) << x_aug + B.col(i);
    Xsig_aug.col(i + n_aug_ + 1) << x_aug - B.col(i);
  }

  double delta_t_sq = delta_t * delta_t;

  // Prediction
  for (int i = 0; i < Xsig_pred_.cols(); ++i) {
    double px = Xsig_aug(0, i);
    double py = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yaw_d = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    if (yaw_d != 0.0) {
      Xsig_pred_(0, i) = px +
                         v / yaw_d * (sin(yaw + yaw_d * delta_t) - sin(yaw)) +
                         0.5 * delta_t_sq * cos(yaw) * nu_a;
      Xsig_pred_(1, i) = py +
                         v / yaw_d * (-cos(yaw + yaw_d * delta_t) + cos(yaw)) +
                         0.5 * delta_t_sq * sin(yaw) * nu_a;
      Xsig_pred_(2, i) = v + delta_t * nu_a;
      Xsig_pred_(3, i) = yaw + yaw_d * delta_t + 0.5 * delta_t_sq * nu_yawdd;
      Xsig_pred_(4, i) = yaw_d + delta_t * nu_yawdd;
    } else {
      Xsig_pred_(0, i) =
          px + v * cos(yaw) * delta_t + 0.5 * delta_t_sq * cos(yaw) * nu_a;
      Xsig_pred_(1, i) =
          py + v * sin(yaw) * delta_t + 0.5 * delta_t_sq * sin(yaw) * nu_a;
      Xsig_pred_(2, i) = v + delta_t * nu_a;
      Xsig_pred_(3, i) = yaw + 0.5 * delta_t_sq * nu_yawdd;
      Xsig_pred_(4, i) = delta_t * nu_yawdd;
    }
  }

  // predict state mean
  x_.setZero();
  for (int j = 0; j < Xsig_pred_.cols(); ++j) {
    x_ += weights_(j) * Xsig_pred_.col(j);
  }
  // predict state covariance matrix
  P_.setZero();
  for (int i = 0; i < Xsig_pred_.cols(); ++i) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    cout << "x_diff " << x_diff << endl;
    while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;
    P_ += weights_(i) * x_diff * x_diff.transpose();
  }
}

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
  assert(meas_package.raw_measurements_.size() == 2);
  VectorXd z = meas_package.raw_measurements_;

  // Predict Lidar measurement
  MatrixXd Zsig = MatrixXd(2, 2 * n_aug_ + 1);
  VectorXd z_pred = VectorXd(2);
  z_pred.setZero();

  // transform sigma points into measurement space
  for (int i = 0; i < Xsig_pred_.cols(); ++i) {
    double px = Xsig_pred_.col(i)(0);
    double py = Xsig_pred_.col(i)(1);

    Zsig.col(i)(0) = px;
    Zsig.col(i)(1) = py;
    z_pred += weights_(i) * Zsig.col(i);
  }
  // calculate measurement covariance matrix S
  MatrixXd S = MatrixXd(2, 2);
  S.setZero();
  for (int i = 0; i < Zsig.cols(); ++i) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S += weights_(i) * z_diff * z_diff.transpose();
  }

  S(0, 0) += std_laspx_ * std_laspx_;
  S(1, 1) += std_laspy_ * std_laspy_;
  MatrixXd Tc = MatrixXd(n_x_, 2);
  Tc.setZero();
  // calculate cross correlation matrix
  for (int i = 0; i < Xsig_pred_.cols(); ++i) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;
    VectorXd z_diff = Zsig.col(i) - z_pred;
    Tc += weights_(i) * x_diff * z_diff.transpose();
  }
  // calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // update state mean and covariance matrix
  VectorXd z_diff = z - z_pred;
  x_ += K * z_diff;
  P_ -= K * S * K.transpose();
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  assert(meas_package.raw_measurements_.size() == 3);
  VectorXd z = meas_package.raw_measurements_;
  // Predict radar measurment
  MatrixXd Zsig = MatrixXd(3, 2 * n_aug_ + 1);
  VectorXd z_pred = VectorXd(3);
  z_pred.setZero();

  // transform sigma points into measurement space
  for (int i = 0; i < Xsig_pred_.cols(); ++i) {
    double px = Xsig_pred_.col(i)(0);
    double py = Xsig_pred_.col(i)(1);
    double v = Xsig_pred_.col(i)(2);
    double yaw = Xsig_pred_.col(i)(3);
    double yaw_d = Xsig_pred_.col(i)(4);

    Zsig.col(i)(0) = sqrt(px * px + py * py);
    Zsig.col(i)(1) = atan2(py, px);
    Zsig.col(i)(2) = (px * cos(yaw) * v + py * sin(yaw) * v) / Zsig.col(i)(0);
    z_pred += weights_(i) * Zsig.col(i);
  }
  // calculate measurement covariance matrix S
  MatrixXd S = MatrixXd(3, 3);
  S.setZero();
  for (int i = 0; i < Zsig.cols(); ++i) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;
    S += weights_(i) * z_diff * z_diff.transpose();
  }

  S(0, 0) += std_radr_ * std_radr_;
  S(1, 1) += std_radphi_ * std_radphi_;
  S(2, 2) += std_radrd_ * std_radrd_;
  MatrixXd Tc = MatrixXd(n_x_, 3);
  Tc.setZero();
  // calculate cross correlation matrix
  for (int i = 0; i < Xsig_pred_.cols(); ++i) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // angle normalization
    while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;
    Tc += weights_(i) * x_diff * z_diff.transpose();
  }
  // calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // update state mean and covariance matrix
  VectorXd z_diff = z - z_pred;
  while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
  while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

  x_ += K * z_diff;
  P_ -= K * S * K.transpose();
}
