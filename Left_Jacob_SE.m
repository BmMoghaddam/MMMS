function [JJ_L ,JJ_L_inv]=Left_Jacob_SE(tau)
% This function takes the vector tau and calculates the Right Jacobian of
% SE(3)

Phi=tau(1:3,1);
phi=norm(Phi);
ad_tau=adjoint1_EKF(tau);
 
if phi~=0
   JJ_L=eye(6)+(4-phi*sin(phi)-4*cos(phi))/(2*phi^2)*ad_tau+...
     (4*phi-5*sin(phi)+phi*cos(phi))/(2*phi^3)*ad_tau^2+...
     (2-phi*sin(phi)-2*cos(phi))/(2*phi^4)*ad_tau^3+...
     (2*phi-3*sin(phi)+phi*cos(phi))/(2*phi^5)*ad_tau^4;
else
   JJ_L=eye(6)+1/2*ad_tau;
end

JJ_L_inv=inv(JJ_L);