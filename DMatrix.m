function D=DMatrix(X,Dt,Mass_matrix,mu_vector)
% X: state
% Dt: time increment

% the 6 generative of lie algebra
g1 = [1;0;0;0;0;0];
g2 = [0;1;0;0;0;0];
g3 = [0;0;1;0;0;0];
g4 = [0;0;0;1;0;0];
g5 = [0;0;0;0;1;0];
g6 = [0;0;0;0;0;1];


d1 = inv(Mass_matrix)*Adjoint_EKF(inv(X))'*(-adjoint1_EKF(g1))'*mu_vector;
d2 = inv(Mass_matrix)*Adjoint_EKF(inv(X))'*(-adjoint1_EKF(g2))'*mu_vector;
d3 = inv(Mass_matrix)*Adjoint_EKF(inv(X))'*(-adjoint1_EKF(g3))'*mu_vector;
d4 = inv(Mass_matrix)*Adjoint_EKF(inv(X))'*(-adjoint1_EKF(g4))'*mu_vector;
d5 = inv(Mass_matrix)*Adjoint_EKF(inv(X))'*(-adjoint1_EKF(g5))'*mu_vector;
d6 = inv(Mass_matrix)*Adjoint_EKF(inv(X))'*(-adjoint1_EKF(g6))'*mu_vector;

D=Dt*[d1 d2 d3 d4 d5 d6];




