function ad_tau=adjoint1_EKF(tau)
% This function takes tau, an element of se_4(3), and returns its adjoint

% matrix representation 
Phi  = tau(1:3);     % rotation vector
rho1 = tau(4:6);     % translation vector 



ad_tau=[skew_EKF(Phi)    zeros(3)    
        skew_EKF(rho1)   skew_EKF(Phi)  ];

