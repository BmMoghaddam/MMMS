function [JJ_R ,JJ_R_inv]=Right_Jacob_SE(tau)
% This function takes the vector tau and calculates the Right Jacobian of
% SE(3)

ad_tau=adjoint1_EKF(tau);
Ad_T=expm(ad_tau);

JJ_L =Left_Jacob_SE(tau);

JJ_R=inv(Ad_T)*JJ_L;

JJ_R_inv=inv(JJ_R);

