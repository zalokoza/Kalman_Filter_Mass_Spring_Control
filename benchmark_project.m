function [kalman_gains,P] = custom_kalman(A,B1,C,Rw,Rv,M,N)
kalman_gains = cell(N,1);
     for i=1:N
     G = M*C'/(C*M*C'+Rv);
     P = M - G*C*M;
     M = A*P*A'+B1*Rw*B1';
     kalman_gains{i} = G;
     end
end

function K = optimal_gain(Q,R,A,B,N) % Output optimal gain Matrix K
    K = cell(N-1,1);
    P = Q;
    for i = N :-1:1
        gain = (B'*P*B + R)^-1*B'*P*A;
        P = A'*P*(A-B*gain) + Q;
        K{i} = gain;
    end
    disp('The final P matrix is: ')
    disp(P)
end

function [nr2,dr2] = nice(old,z)
    % NICE takes the original expression and makes it look nice
    % Input: (a*z^2 + b*z + c)/(d*z^2 + e*z + f)
    % Output: (a1*z^2 + b1*z +c1)/(z^2 + e1*z + f1)
    % old contains symbol of z
    % rearrange the expression to collect coefficients
    new = collect(old,z);
    % obtain the coefficients of numerator and denominator
    dr = fliplr(coeffs(feval(symengine,'denom',new)));
    nr = fliplr(coeffs(feval(symengine,'numer',new)));
    % normalize the coefficients with the coeff of dr z^2 term
    dr2 = double(dr/dr(1)); nr2 = double(nr/dr(1));
end

clear

% Continuous, No Control, Disturbance Added

k =  2;
m1 = 1;
m2 = 1;
b = .04;
x0 = [0;0;0;0];
q0 = [0;0;0;0];
%% Continuous, No Control, MISO State Space

A = [0 1 0 0
    -k/m1 -b/m1 k/m1 b/m1
    0 0 0 1
    k/m2 b/m2 -k/m2 -b/m2];

B1 = [0; 1/m1; 0; 0];
B2 = [0; 0; 1; 1/m2];
B = [B1 B2];
C = [0 0 1 0]; % Collocated
% C = [0 0 1 0;
%     0 0 0 1]; % Cheating hehe
% C = eye(4); % Even more cheating
D = 0;
K = [0 0 0 0] % No control effort i.e. uncompensated system
%%
% Discretize the state space model. Valid for a ZoH applied to controller

T = .01;
B = (eye(4)*T + A*T^2/2 + A^2*T^3/6 + A^3*T^4/24)*B;
A = expm(A*T);

% Calculate linear quadratic optimal feedback gain depending on cost function
% Note: use B1 instead of B because we want state feedback control law for
% the linear actuator, NOT for the input channel meant to simulate
% disturbance

Q = 30*eye(4);
R = 1;
N = 25;
K = optimal_gain(Q, R, A, B1, N);
K = K{1} % Pick the optimal gain for k = 1 and use it for all times t = kT

% Calculate Kalman Gain Matrix

B1 = [1;1;1;1]; % Plant uncertainty is entering all of the states
Rw = 1; % Covariance of 1
Rv = 1; % Covariance of 1
M = eye(4); % Assume initial M is identity
N = 25; % See if Kalman Gain Matrix converges in 25 steps

[kalman_gains, P] = custom_kalman(A,B1,C,Rw,Rv,M,N);
G = kalman_gains{25};

% For curiousity, compare previous results with results obtained from
% included MATLAB command "kalman"

sys = ss(A,B,C,D, T,Name = 'Discretized Plant');
[kalmf, L, P] = kalman(sys,Rw, Rv, 0, 'current');

%% Model Redesign: SISO

A = [0 1 0 0
    -k/m1 -b/m1 k/m1 b/m1
    0 0 0 1
    k/m2 b/m2 -k/m2 -b/m2];

B = [0; 1/m1; 0; 0];
C = [0 0 1 0]; % Collocated
% C = [0 0 1 0;
%     0 0 0 1]; % Cheating hehe
% C = eye(4); % Even more cheating
D = 0;

%%
% Discretize the state space model. Valid for a ZoH applied to controller

T = .01;
B = (eye(4)*T + A*T^2/2 + A^2*T^3/6 + A^3*T^4/24)*B;
A = expm(A*T);
sys = ss(A,B,C,D, T,Name = 'Discretized Plant');

% Calculate linear quadratic optimal feedback gain depending on cost function

Q = 10*eye(4);
R = 1;
N = 25;
K = optimal_gain(Q, R, A, B, N);
K = K{1} % Pick the optimal gain for k = 1 and use it for all times t = kT

% Calculate Kalman Gain Matrix

B1 = [1;1;1;1]; % Plant uncertainty is entering all of the states
Rw = 40; % Plant Noise Covariance of 1
Rv = 40; % Sensor Noise Covariance of 1
M = eye(4); % Assume initial M is identity
N = 200; % See if Kalman Gain Matrix converges in 200 steps

[kalman_gains, P] = custom_kalman(A,B,C,Rw,Rv,M,N);
G = kalman_gains{N};

% For curiousity, compare previous results with results obtained from
% included MATLAB command "kalman"

[kalmf, L, P] = kalman(sys,Rw, Rv, 0, 'current');

%% Perform OL Analysis to calculate Stability Margins
z = tf('z', T);
G_z = C/(z*eye(4)-A)*B;
[num, den, T] = tfdata(sys)
G2 = tf(num,den, T);
%G = minreal(G)
D_LQG = z*K/(z*eye(4)-A + G*C*A + B*K - G*C*B*K)*G;
syms z, T = .01;
D_LQG2 = z*K/(z*eye(4)-A + G*C*A + B*K - G*C*B*K)*G;
[n, d] = nice(D_LQG2, z);
D_LQG2 = tf(n, d, T);
%D_LQG = minreal(D_LQG)
figure;
margin(D_LQG2*G2);

%% Check Robustness of Model by varying parameters

% Change k without changing the kalman or feedback gains

k = 50; % Test robustness by changing k here

A = [0 1 0 0
    -k/m1 -b/m1 k/m1 b/m1
    0 0 0 1
    k/m2 b/m2 -k/m2 -b/m2];

B = [0; 1/m1; 0; 0];

% Don't forget to re-discretize it

B = (eye(4)*T + A*T^2/2 + A^2*T^3/6 + A^3*T^4/24)*B;

A = expm(A*T); 

