% Make the inertial parameters uncertain


% End-Efffector
m0_unc=robot.Bodies{1,n+1}.Mass.*(0.8+rand*0.4); %up to 20% uncertainty

% arm + spacecraft in reverse
mm_unc=zeros(1,n);
mm_unc(n)=10.*(0.8+rand*0.4);%kg spacecraft

for i=2:n
    mm_unc(n+1-i)=robot.Bodies{1,i}.Mass.*(0.8+rand*0.4);
end

% ----------------- inertia matrices
I0_unc(1:3,1:3)=[robot.Bodies{1,8}.Inertia(1) robot.Bodies{1,8}.Inertia(6) robot.Bodies{1,8}.Inertia(5); ...
    robot.Bodies{1,8}.Inertia(6) robot.Bodies{1,8}.Inertia(2) robot.Bodies{1,8}.Inertia(4); ...
    robot.Bodies{1,8}.Inertia(5) robot.Bodies{1,8}.Inertia(4) robot.Bodies{1,8}.Inertia(3)].*(m0_unc/robot.Bodies{1,n+1}.Mass);

%spacecraft
Im_unc(1:3,1:3,n)=eye(3)*1.66667.*(mm_unc(n)/10); %for cube with density 0.8
%arm
for i=2:n
    Im_unc(1:3,1:3,n-i+1)=[robot.Bodies{1,i}.Inertia(1) robot.Bodies{1,i}.Inertia(6) robot.Bodies{1,i}.Inertia(5); ...
    robot.Bodies{1,i}.Inertia(6) robot.Bodies{1,i}.Inertia(2) robot.Bodies{1,i}.Inertia(4); ...
    robot.Bodies{1,i}.Inertia(5) robot.Bodies{1,i}.Inertia(4) robot.Bodies{1,i}.Inertia(3)].*(mm_unc(n+1-i)/robot.Bodies{1,n+1-i}.Mass);
end


% calculate Inertia matrices in the joint frames
[M_curly0_unc,M_curlym_unc]=M_curly_ee(m0_unc,I0_unc,mm_unc,Im_unc,Ad_gcm_inv,Ad_gcm0_inv);

% ************** Initiate forces

% f_0=[0;0;0;0;0;0];
% f_m=[0;0;0;0;0;0;0];
% f_e=[0;0;0;0;0;0];



% Form the math matrix diagonalized in the base frame
clear M_frak
[M_frak0_unc,M_frak_unc] =M_frak(M_curly0_unc,M_curlym_unc,Ad_gbar_inv);

diag_M_unc =diagonalize(M_frak0_unc, M_frak_unc);


m0_frak=diag_M_unc(1:6,1:6);
m_diag=diag_M_unc(7:end,7:end);

mu_unc=zeros(6,1);
mu_t_unc=zeros(6,1);

g_I0=eye(4);

% *************************** Set Target parameters

mt_unc=10*(0.8+rand*0.4);
M_t_body_unc=[eye(3)*mt_unc zeros(3); zeros(3) eye(3)*0.416667].*(mt_unc/10);

M_t_unc=(((Adjoint(g_It))'\M_t_body_unc)/(Adjoint(g_It)));
mu_t=inv((Adjoint(g_It))')*M_t_body_unc*[V_It;w_t];

M_temp_unc=sim('M0_initiate',0.01);
M0_unc=M_temp_unc.M0.data(:,:,1);%iota0'*M_frak0*iota0;
M0m_unc=M_temp_unc.M0m.data(:,:,1);
% M0m=M_temp.M0m.data(:,:,1);
P0_unc=M0_unc*V_I0+(M0_unc\M0m_unc)*q_dot_m;
mu=P0_unc;

