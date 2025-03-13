function W=skew_EKF(w)
% This function takes a 3*1 vector and make a 3*3 skew symmetric matrix


     w1=w(1,1);
     w2=w(2,1);
     w3=w(3,1);

     W=[0    -w3     w2
        w3    0     -w1
       -w2    w1     0];
 