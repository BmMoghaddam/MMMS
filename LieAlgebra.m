function zeta=LieAlgebra(X)
% this function gets an element of Lie group SE_k(n) and returns its
% corresponding Lie algebra element
% X is an element of SE_k(n)
% zeta is an element of Lie algebra


    R  = X(1:3,1:3);
    t1 = X(1:3,4);
    
    Phi=logm(R); 
    phi=vee_EKF(Phi);

    J_L =Left_Jacob_SO(phi);

    rho=inv(J_L)*t1;

    zeta=[phi; rho ];




