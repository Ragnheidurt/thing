%% Question 1
% Four states:
% State 1: Both engines work
% State 2: Engine 2 out but 1 still works
% State 3: Engine 1 out but 2 still works
% State 4: Both engines out

%% Question 2

syms lambda1 lambda2 mu1 mu2
intensity1 = [-(lambda1+lambda2) lambda2 lambda1 0;
    3*mu2 -(lambda1+mu2*3) 0 lambda1;
    mu1*3 0 -(lambda2+3*mu1) lambda2;
    0 3*mu1 0*mu2 -(3*mu1+0*mu2)];

intensity2 = [-(lambda1+lambda2) lambda2 lambda1 0;
    3*mu2 -(lambda1+mu2*3) 0 lambda1;
    mu1*3 0 -(lambda2+3*mu1) lambda2;
    0 2*mu1 1*mu2 -(2*mu1+1*mu2)];

intensity3 = [-(lambda1+lambda2) lambda2 lambda1 0;
    3*mu2 -(lambda1+mu2*3) 0 lambda1;
    mu1*3 0 -(lambda2+3*mu1) lambda2;
    0 1*mu1 2*mu2 -(1*mu1+2*mu2)];

intensity4 = [-(lambda1+lambda2) lambda2 lambda1 0;
    3*mu2 -(lambda1+mu2*3) 0 lambda1;
    mu1*3 0 -(lambda2+3*mu1) lambda2;
    0 0*mu1 3*mu2 -(0*mu1+3*mu2)];

%% Question 3
% Answered on paper

%% Question 4

syms pi1 pi2 pi3 pi4
pi = [pi1; pi2; pi3; pi4];
solve(intensity1*pi == 0)


%% Question 5



