function [J_L, J_L_inv]=Left_Jacob_SO(Phi)
% This function takes the vector Phi and calculates the Left Jacobian of
% SO(3)

[m,n]=size(Phi);

if m==3 & n==1

phi=norm(Phi);

Phi_wedge=skew_EKF(Phi);

if phi~=0

J_L     = eye(3)+(1-cos(phi))/phi^2*Phi_wedge+(phi-sin(phi))/phi^3*Phi_wedge^2;
J_L_inv = eye(3)-1/2*Phi_wedge+(1/phi^2-(1+cos(phi))/(2*phi*sin(phi)))*Phi_wedge^2;
else
    J_L     = eye(3);
    J_L_inv = eye(3);
end

else
    error('Error: the input must be a 3 by 1 vector')

end