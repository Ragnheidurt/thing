birthdate = 19980422;  % Write the birth date on format yyyymmdd for oldest member
format compact
[lambda1,lambda2,mu1,mu2,V1,V2,V] = getFerrydata(birthdate);  
h=0.001; % Discretization step

%% ANALYTIC SOLUTION
% Question 1 should be answered in the report

% Question 2 should be answered in the report, and submitted below

Qi = [-(lambda1+lambda2) lambda2 lambda1 0;
    3*mu2 -(3*mu2+lambda1) 0 lambda1;
    3*mu1 0 -(3*mu1+lambda2) lambda2;
    0 3*mu1 0 -3*mu1];
Qieig=sort(eig(Qi)); % We compare the eigenvalues
% To make sure that the order of your states will not change the result
Qii = [-(lambda1+lambda2) lambda2 lambda1 0;
    3*mu2 -(3*mu2+lambda1) 0 lambda1;
    3*mu1 0 -(3*mu1+lambda2) lambda2;
    0 2*mu1 1*mu2 -(2*mu1+1*mu2)];
Qiieig=sort(eig(Qii));

Qiii = [-(lambda1+lambda2) lambda2 lambda1 0;
    3*mu2 -(3*mu2+lambda1) 0 lambda1;
    3*mu1 0 -(3*mu1+lambda2) lambda2;
    0 1*mu1 2*mu2 -(1*mu1+2*mu2)];
Qiiieig=sort(eig(Qiii));

Qiv = [-(lambda1+lambda2) lambda2 lambda1 0;
    3*mu2 -(3*mu2+lambda1) 0 lambda1;
    3*mu1 0 -(3*mu1+lambda2) lambda2;
    0 0 3*mu2 -3*mu2];
Qiveig=sort(eig(Qiv));

% Question 3 should be answered in the report

% Question 4 should be answered in the report, describe how you do it, and  

vect = sym('pi', [1 4]);


PIi_v = vpasolve([vect*Qi == 0, sum(vect) == 1],vect) % Stationary Distributions Case (i)
PIi = double([PIi_v.pi1, PIi_v.pi2, PIi_v.pi3, PIi_v.pi4]);
PIisort=sort(PIi) % We compare the sorted vectors

PIii_v = vpasolve([vect*Qii == 0, sum(vect) == 1],vect) % Stationary Distributions Case (ii)
PIii = double([PIii_v.pi1, PIii_v.pi2, PIii_v.pi3, PIii_v.pi4]);
PIiisort=sort(PIii);

PIiii_v = vpasolve([vect*Qiii == 0, sum(vect) == 1],vect) % Stationary Distributions Case (iii)
PIiii = double([PIiii_v.pi1, PIiii_v.pi2, PIiii_v.pi3, PIiii_v.pi4]);
PIiiisort=sort(PIiii);

PIiv_v = vpasolve([vect*Qiv == 0, sum(vect) == 1],vect) % Stationary Distributions Case (iv)
PIiv = double([PIiv_v.pi1, PIiv_v.pi2, PIiv_v.pi3, PIiv_v.pi4]);
PIivsort=sort(PIiv)

% Question 5 should be answered in the report, describe how you do it, and  

velo_vect = [V, V1, V2, 0];

AVi = velo_vect(1)*double(PIi_v.pi1) + velo_vect(2)*double(PIi_v.pi2) + velo_vect(3)*double(PIi_v.pi3) + velo_vect(4)*double(PIi_v.pi4); % Average Ferry Velocity for Case (i) 
AVii = velo_vect(1)*double(PIii_v.pi1) + velo_vect(2)*double(PIii_v.pi2) + velo_vect(3)*double(PIii_v.pi3) + velo_vect(4)*double(PIii_v.pi4); % Average Ferry Velocity for Case (ii) 
AViii = velo_vect(1)*double(PIiii_v.pi1) + velo_vect(2)*double(PIiii_v.pi2) + velo_vect(3)*double(PIiii_v.pi3) + velo_vect(4)*double(PIiii_v.pi4); % Average Ferry Velocity for Case (iii) 
AViv = velo_vect(1)*double(PIiv_v.pi1) + velo_vect(2)*double(PIiv_v.pi2) + velo_vect(3)*double(PIiv_v.pi3) + velo_vect(4)*double(PIiv_v.pi4); % Average Ferry Velocity for Case (iv) 

AV = [AVi AVii AViii AViv]
class(AV)

%% CONTINUOUS TIME SIMULATION
% Question 6a should be answered in the report
% Start in State V (Velocity = V, both engines working)
% From State V, we can either jump to state V1 (Engine 2 breaks) or V2
% (Engine 1 breaks)
% The dwell times can be determined from the for loop
% The total jump time is the sum of dwell times and probabilities
    
T = 50000;
S = 1; % initial State = 1
C = [0 0 0 0]; % Counter for all states

Qnew = zeros(4,1);
for i=1:4
    for j=1:4
        if j~=i
            Qnew(i,1) = Qnew(i,1) + Qi(i,j);
        end
    end
end
Qi
Qnew

Pnew = [Qi(1,2)/Qnew(1,1); Qi(2,1)/Qnew(2,1); Qi(3,1)/Qnew(3,1); Qi(4,2)/Qnew(4,1)]
transition = [2 3; 1 4; 1 4; 2 3];
Ti = zeros(4,1);
i = 0;
while i<T
    Ti = exprnd(1./Qnew);
    prop = rand();
    if prop<Pnew(S)
        S = transition(S,1);
    else
        S = transition(S,2); 
    end
    C_t = 0;
    C(S) = C(S) + Ti(S);
    i = i + Ti(S);
end

Ti
C

AV_cont = [];

for i = 1:4
    AV_cont(i) = C(i)/T*velo_vect(i)
end

Bones = sum(AV_cont)

%     
% 
% 
% % Question 7a should be answered in the report, describe how you do it, and check that the result agrees with the analytic result
% 
% % Question 8a should be answered in the report, describe how you do it, do the calculations and enter results below
%%

S2 = 1;
Ti2 = zeros(4,1);
C2 = [0 0 0 0];

C_t = 0;
i=0;
while i<T
    Ti2 = unifrnd(zeros(4,1),2*1./Qnew);
    prop = rand();
    if prop<Pnew(S2,1)
        S2 = transition(S2,1);
    else
        S2 = transition(S2,2); 
    end
    C2(S2) = C2(S2) +Ti2(S2);
    i = i+Ti2(S2);
end

Ti2
C2
AV_cont = [];

for i = 1:4
    AV_cont(i) = C2(i)/T*velo_vect(i)
end

AVuni = sum(AV_cont)

%%

Ti3 = zeros(4,1);
S3 = 1;
C3 = [0 0 0 0];
C_t3 = 0;
i = 0;
while i<T
    Ti3 = 1./Qnew;
    prop = rand();
    if prop<Pnew(S3)
        S3 = transition(S3,1);
    else
        S3 = transition(S3,2); 
    end
    C3(S3) = C3(S3) + Ti3(S3);
    i = i + Ti3(S3);
end

Ti3
C3

AV_cont = [];

for i = 1:4
    AV_cont(i) = C3(i)/T*velo_vect(i)
end

AVdet = sum(AV_cont)


%% DISCRETE TIME SIMULATION
% Question 6b should be answered in the report, describe how you determine the transition matrix and enter below

% Question 7b should be answered in the report, describe how you do it, and check that the result agrees with the analytic result

% Question 8b should be answered in the report, describe how you do it, do the calculations and enter results below

I = eye(4);

P = I+h*Qi

transition = [2 3; 1 4; 1 4; 2 3];
T = 10000;
Ti = zeros(4,1);
S4 = 1;
C4 = [0 0 0 0];
C_t4 = 0; 
C5 = []; % Counter for steps until state 4 is reached, provided that it was not in state 4 before
C6 = 0; % Counter for times system changed to step 4
C7 = 0; % Counter for the probability of failure after 10 steps
Probfail10 = 0; % Initialization of Probfail10
for t = 0:0.001:T
    Random = rand();
    if Random<P(S4,S4)
        S4 = S4;
    elseif Random<P(S4,transition(S4,1))+P(S4,S4)
        S4 = transition(S4,1);
        if S4 == 4
%             C5 = C5 + C_t4/0.001;
%             C_t4 = 0;
%             C6 = C6 + 1;
             C5(end+1) = C_t4/0.001; % Counter for 9b
             C_t4 = 0; % Counter for 9b
        end
    else
        S4 = transition(S4,2);
        if S4 == 4
%             C5 = C5 + C_t4/0.001;
%             C_t4 = 0;
%             C6 = C6 + 1;
             C5(end+1) = C_t4/0.001;
             C_t4 = 0;
        end
    end
    C4(S4) = C4(S4) + 0.001;
    C_t4 = C_t4 + 0.001;
end

C4
AV_disc = [];

for i = 1:4
    AV_disc(i) = C4(i)/T*velo_vect(i)
end

AVsup = sum(AV_disc)


% Question 9b should be answered in the report, describe how you do it, do the calculations and enter results below

Probfail10 = P(1,2)*P(2,4)+P(1,3)*P(3,4)
% ETtTF = mean(C5) % Implementation see lines 263ff. and 272ff.

T4 = [1 -lambda2/(lambda1+lambda2) -lambda1/(lambda1+lambda2) 1/(lambda1+lambda2);
    3*mu2/(3*mu2+lambda1) -1 0 -1/(3*mu2+lambda1);
    3*mu1/(3*mu1+lambda2) 0 -1 -1/(3*mu1+lambda2)];

T4i = rref(T4);
ETtTF = T4i(1,4)/h

V
V1
V2