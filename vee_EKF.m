function w=vee_EKF(W)
% This function gets a skew symmetric matrix W and applies vee operator 
% on it to generate the vector w


    p=W(3,2);
    q=W(1,3);
    r=W(2,1);

    w=[p;q;r];

