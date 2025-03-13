
%clear all;
%close all;
%clc;

%
% covariance of process noise (symmetric and positive definite) 
Q=0.00001.*eye(6);
% covariance of measurement noise (symmetric and Positive definite)
R=0.00001.*eye(6);

% initial value of state
X0=[eye(3) [0;0;0] 
    zeros(1,3) 1 ];   

% initial value of covariance (symmetric and positive definite)
P0=0.00001.*eye(6); % usually a big number which means we are not confident at the begining

% assigning initial conditions to k=1
X_pos(:,:,1) = X0; % the posterior estimate of state 
P_pos(:,:,1) = P0; % the posterior estimate of covariance 

% mass matrix
Mass_matrix = 10*eye(6);

% mu vector
mu_vector = [1;0;0;0;0;0];

% time increment
 Dt=0.5;

% total time steps; time=k*Dt  where k=1:total_time_steps
total_time_steps=100;



%% Generating noise signal
% m: measurment noise; in your code you set it with a random normal variable
% it is zero-mean and its covariance is R

s  = RandStream.create('mrg32k3a','NumStreams',1);

% measurement noise signal  
m  = randn(s,total_time_steps,6)*chol(R); 
    
% % % Ground truth; comment this if you have it
 X_gt(:,:,1)=X0; % ground truth state
for k=2:total_time_steps
   Omega_gt(k-1,:) = Dt*inv(Mass_matrix)*Adjoint_EKF(inv(X_gt(:,:,k-1)))*mu_vector; % generalized velocity
   X_gt(:,:,k)     = X_gt(:,:,k-1)*expm(wedge_EKF(Omega_gt(k-1,:)'));
end

% % simulated measurments in case you do not have measurements; comment this if you have it 
for k=1:total_time_steps
    Measured_data(:,:,k) = h_fcn(X_gt(:,:,k))*expm(wedge_EKF(m(k,:)'));
end

%% EKF iterations

for k=2:total_time_steps
    
    z(:,:,k)   = Measured_data(:,:,k);  % measurement vector 
    
    % z(:,:,k)=0; % it is just for test; comment this line
    % generalized velocity
    Omega(:,k-1) = Dt*inv(Mass_matrix)*Adjoint_EKF(inv(X_pos(:,:,k-1)))*mu_vector; 
    % matrix D for linearization
    D(:,:,k-1)   = DMatrix(X_pos(:,:,k-1),Dt,Mass_matrix,mu_vector); 
    % matrix F for linearization
    F(:,:,k-1)   = Adjoint_EKF(expm(wedge_EKF(-Omega(:,k-1))))...
                   +Right_Jacob_SE(Omega(:,k-1))*D(:,:,k-1); 
   
    % EKF propagation step
    % the prior estimate of state
    X_pri(:,:,k) = X_pos(:,:,k-1)*expm(wedge_EKF(Omega(:,k-1))); 
    % the prior covariance matrix
    P_pri(:,:,k) = F(:,:,k-1)*P_pos(:,:,k-1)*F(:,:,k-1)'...
                   +Right_Jacob_SE(Omega(:,k-1))*Q*Right_Jacob_SE(Omega(:,k-1))'; 

    % EKF update step
    % innovation vector
    nu(:,k)        = LieAlgebra(inv(h_fcn(X_pri(:,:,k)))*z(:,:,k)); 

    % H matrix; since h(X)=X
    H(:,:,k) = [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1];

    % EKF gain
    K(:,:,k)     = P_pri(:,:,k)*H(:,:,k)'*inv(H(:,:,k)*P_pri(:,:,k)*H(:,:,k)'+R); 

    % posterior estimate
    X_pos(:,:,k) = X_pri(:,:,k)*expm(wedge_EKF( K(:,:,k)*nu(:,k))) ; 
    % posterior covariance
    P_pos(:,:,k) = Right_Jacob_SE(K(:,:,k)*nu(:,k))*(eye(6)-K(:,:,k)*H(:,:,k))...
                   *P_pri(:,:,k)*Right_Jacob_SE(K(:,:,k)*nu(:,k))'; 
end

for i=1:total_time_steps
    EE(i)=norm((X_pos(:,:,i)\X_gt(:,:,i))-eye(4));
    EE_noise(i)=norm((Measured_data(:,:,i)\X_gt(:,:,i))-eye(4));
    alpha(i)=acos(trace(X_gt(1:3,1:3,i))/2-1/2);
    alpha_noise(i)=acos(trace(Measured_data(1:3,1:3,i))/2-1/2);
    alpha_pos(i)=acos(trace(X_pos(1:3,1:3,i))/2-1/2);
end

%%
g_It00=[eye(3) rho_It;0 0 0 1];
Dt=0.1;
Dt_s=0.1;
noise_st=0.000005;