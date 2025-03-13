function [J_R ,J_R_inv]=Right_Jacob_SO(Phi)
% This function takes the vector Phi and calculates the Right Jacobian of
% SO(3)

phi=norm(Phi);

Phi_wedge=skew_EKF(Phi);

if phi~=0

J_R     = eye(3)-(1-cos(phi))/phi^2*Phi_wedge+(phi-sin(phi))/phi^3*Phi_wedge^2;
J_R_inv = eye(3)+1/2*Phi_wedge+(1/phi^2-(1+cos(phi))/(2*phi*sin(phi)))*Phi_wedge^2;
else
    J_R     = eye(3);
    J_R_inv = eye(3);
end
