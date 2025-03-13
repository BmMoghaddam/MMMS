function W=wedge_EKF(w)
% This function gets the vector w and applies wedge operator of group
% SE(3) on it to get matrix W which is an element of Lie algebra

phi=w(1:3,1);
rho1=w(4:6,1);


Phi=skew_EKF(phi);

W=[Phi         rho1  
   zeros(1,3)  0  ];

 
     