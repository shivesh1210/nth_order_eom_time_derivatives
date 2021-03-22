%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author: Andreas Mueller
% last modified: August 4, 2020
% remark: this code is by no means optimized

function [Q,Qd,Q2d] = InvDyn(q,qd,q2d,q3d,q4d)

global Chain;
global Param;

global n; % DOF, number of joints
global g; % gravity vector


%% Forward kinematics recursion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first body is treated separately since there is no predecessor
Chain(1).f = SE3Exp(Param(1).Y,q(1));
Chain(1).C = Chain(1).f*Param(1).A;
Chain(1).S = SE3AdjMatrix(Chain(1).f)*Param(1).Y;
Chain(1).V = Chain(1).S*qd(1);
Chain(1).Sd = SE3adMatrix(Chain(1).V)*Chain(1).S;
Chain(1).Vd = Chain(1).S*q2d(1) + Chain(1).Sd*qd(1) ; %+ [[0;0;0];g];
Chain(1).S2d = SE3adMatrix(Chain(1).Vd)*Chain(1).S + SE3adMatrix(Chain(1).V)*SE3adMatrix(Chain(1).V)*Chain(1).S;
Chain(1).V2d = Chain(1).S*q3d(1) + 2*Chain(1).Sd*q2d(1) + Chain(1).S2d*qd(1);
Chain(1).S3d = SE3adMatrix(Chain(1).V2d)*Chain(1).S + 2*SE3adMatrix(Chain(1).Vd)*SE3adMatrix(Chain(1).V)*Chain(1).S ...
    + SE3adMatrix(Chain(1).V)*SE3adMatrix(Chain(1).Vd)*Chain(1).S ...
    + SE3adMatrix(Chain(1).V)*SE3adMatrix(Chain(1).V)*SE3adMatrix(Chain(1).V)*Chain(1).S;
Chain(1).V3d = Chain(1).S*q4d(1) + 3*Chain(1).Sd*q3d(1) + 3*Chain(1).S2d*q2d(1) + Chain(1).S3d*qd(1);

% the forward recursion run
for i=2:n;
    Chain(i).f = Chain(i-1).f*SE3Exp(Param(i).Y,q(i));
    Chain(i).C = Chain(i).f*Param(i).A;
    Chain(i).S = SE3AdjMatrix(Chain(i).f)*Param(i).Y;
    Chain(i).V = Chain(i-1).V + Chain(i).S*qd(i);
    Chain(i).Sd = SE3adMatrix(Chain(i).V)*Chain(i).S;
    Chain(i).Vd = Chain(i-1).Vd + Chain(i).S*q2d(i) + Chain(i).Sd*qd(i);
%     Chain(i).S2d = SE3adMatrix(Chain(i).Vd)*Chain(i).S + SE3adMatrix(Chain(i).V)*SE3adMatrix(Chain(i).V)*Chain(i).S;
    Chain(i).S2d = (SE3adMatrix(Chain(i).Vd) + SE3adMatrix(Chain(i).V)*SE3adMatrix(Chain(i).V))*Chain(i).S;
    Chain(i).V2d = Chain(i-1).V2d + Chain(i).S*q3d(i) + 2*Chain(i).Sd*q2d(i) + Chain(i).S2d*qd(i);
%     Chain(i).S3d = SE3adMatrix(Chain(i).V2d)*Chain(i).S + 2*SE3adMatrix(Chain(i).Vd)*SE3adMatrix(Chain(i).V)*Chain(i).S ...
%         + SE3adMatrix(Chain(i).V)*SE3adMatrix(Chain(i).Vd)*Chain(i).S ...
%         + SE3adMatrix(Chain(i).V)*SE3adMatrix(Chain(i).V)*SE3adMatrix(Chain(i).V)*Chain(i).S;
    Chain(i).S3d = (SE3adMatrix(Chain(i).V2d) + 2*SE3adMatrix(Chain(i).Vd)*SE3adMatrix(Chain(i).V) ...
        + SE3adMatrix(Chain(i).V)*SE3adMatrix(Chain(i).Vd) ...
        + SE3adMatrix(Chain(i).V)*SE3adMatrix(Chain(i).V)*SE3adMatrix(Chain(i).V))*Chain(i).S;

    Chain(i).V3d = Chain(i-1).V3d + Chain(i).S*q4d(i) + 3*Chain(i).Sd*q3d(i) + 3*Chain(i).S2d*q2d(i) + Chain(i).S3d*qd(i);
end;

%% Inverse dynamics recursion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initializations of auxiliary variables
M = zeros(6,6);
AdInvC = zeros(6,6);
Pi = zeros(6,1);
Pid = zeros(6,1);
Pi2d = zeros(6,1);
Pi3d = zeros(6,1);
adV = zeros(6,6);
adVd = zeros(6,6);
adV2d = zeros(6,6);

% the backward recursion run
for i=n:-1:1
    AdInvC = SE3AdjInvMatrix(Chain(i).C);
    adV = SE3adMatrix(Chain(i).V);
    adVd = SE3adMatrix(Chain(i).Vd);
    adV2d = SE3adMatrix(Chain(i).V2d);
    adV3d = SE3adMatrix(Chain(i).V3d);

    M = AdInvC'*Param(i).Mb*AdInvC;
    Pi = M*Chain(i).V;
    Pid = M*Chain(i).Vd - adV'*Pi;
    Pi2d = M*(Chain(i).V2d - adV*Chain(i).Vd) ...
           -2*adV'*Pid - adVd'*Pi - adV'*adV'*Pi;
    Pi3d = M*(Chain(i).V3d - 2*adV*Chain(i).V2d + adV*adV*Chain(i).Vd) ...
           -3*adV'*Pi2d - 3*(adVd + adV*adV)'*Pid ...
           -(adV2d + 2*adVd*adV + adV*adVd + adV*adV*adV)'*Pi;
    if i == n
        Chain(i).W = Pid;
        Chain(i).Wd = Pi2d;
        Chain(i).W2d = Pi3d;
    else
        Chain(i).W = Chain(i+1).W + Pid;
        Chain(i).Wd = Chain(i+1).Wd + Pi2d;
        Chain(i).W2d = Chain(i+1).W2d + Pi3d;
    end
    Chain(i).Q = Chain(i).S'*Chain(i).W;
    Chain(i).Qd = Chain(i).S'*Chain(i).Wd + Chain(i).Sd'*Chain(i).W;
    Chain(i).Q2d = Chain(i).S'*Chain(i).W2d + Chain(i).S2d'*Chain(i).W + 2*Chain(i).Sd'*Chain(i).Wd;
    Q(i) = Chain(i).Q;
    Qd(i) = Chain(i).Qd;
    Q2d(i) = Chain(i).Q2d;
end;

end
