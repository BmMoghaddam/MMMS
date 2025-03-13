function Ad_T=Adjoint_EKF(X)
% This function takes X, an element of SE_4(3), and returns its Adjoint
% matrix representation 
R  = X(1:3,1:3); % rotation matrix
t1 = X(1:3,4);   % translation vector 1 (velocity)


Ad_T=[R            zeros(3)   
      skew_EKF(t1)*R   R      ];

end

