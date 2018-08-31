#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
    x_ = x_in;
    P_ = P_in;
    F_ = F_in;
    H_ = H_in;
    R_ = R_in;
    Q_ = Q_in;
}

void KalmanFilter::Predict() {
    x_ = F_*x_;
    MatrixXd Ft = F_.transpose();
    
    P_ = F_* P_*Ft + Q_;
    
}

void KalmanFilter::Update(const VectorXd &z) {
    VectorXd y = z - H_ * x_;
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd K =  P_ * Ht * Si;
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    
    P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    float x = x_(0);
    float y = x_(1);
    float v_x = x_(2);
    float v_y = x_(3);
    
    float rho = sqrt(x*x + y*y);
    if (fabs(rho) < .00001) {
        x += .001;
        y += .001;
        rho = sqrt(x*x + y*y);
    }
    float theta = atan2f(y, x);
    float ro_dot = (x* v_x + y*v_y)/rho;
    VectorXd z_pred = VectorXd(3);
    z_pred << rho, theta, ro_dot;
    VectorXd y_pred = z - z_pred;
    if (y_pred(1) > M_PI ) {
        y_pred(1) = y_pred(1) - 2 * M_PI;
    } else if (y_pred(1) < -1* M_PI) {
        y_pred(1) = y_pred(1) + 2 * M_PI;
    }
    
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_*P_*Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd K =  P_ * Ht * Si;
    x_ = x_ + (K * y_pred);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K*H_)* P_;
    
}
