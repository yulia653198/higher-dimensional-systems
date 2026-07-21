%% Compute the feedback matrix K for the sixth-order agent model
clc;
clear;

%% Sixth-order system matrices
damping = 0.25;

Omega = [ 0     0.08 -0.04;
         -0.08  0     0.06;
          0.04 -0.06  0   ];

A = [zeros(3), eye(3);
     zeros(3), Omega-damping*eye(3)];

B = [zeros(3);
     eye(3)];

%% CARE weighting matrices
Q = eye(6);
R = eye(3);

%% Controllability check
n = size(A,1);
controllabilityRank = rank(ctrb(A,B));

fprintf('rank(ctrb(A,B)) = %d (system order n = %d)\n', ...
    controllabilityRank,n);

if controllabilityRank < n
    error('The pair (A,B) is not controllable.');
end

%% Solve the continuous-time algebraic Riccati equation
% A''P + PA - PBR^{-1}B''P + Q = 0.
[P,closedLoopPolesFromCare,Kcare] = care(A,B,Q,R);

% Since R=I_3, this is also K=B''P, as used in the manuscript.
K = R\(B'*P);

%% Independent closed-loop check
Acl = A-B*K;
closedLoopEigenvalues = eig(Acl);

disp('A =');
disp(A);

disp('B =');
disp(B);

disp('P =');
disp(P);

disp('K = R^{-1}B^T P =');
disp(K);

disp('eig(A-BK) =');
disp(closedLoopEigenvalues);

disp('Closed-loop poles returned by care =');
disp(closedLoopPolesFromCare);

disp('Feedback matrix returned by care =');
disp(Kcare);

if all(real(closedLoopEigenvalues) < 0)
    fprintf('The matrix A-BK is Hurwitz.\n');
else
    warning('The matrix A-BK is not Hurwitz.');
end
