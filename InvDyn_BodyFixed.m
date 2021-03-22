%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author: Andreas Mueller
% last modified: August 4, 2020
% remark: this code is by no means optimized

function [Q,Qd,Q2d] = InvDyn_BodyFixed(q,qd,q2d,q3d,q4d)

global Chain;
global Param;

global n; % DOF, number of joints
global g; % gravity vector

%% Forward kinematics recursion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first body is treated separately since there is no predecessor
Chain(1).f = SE3Exp(Param(1).Y,q(1));
Chain(1).C = Chain(1).f*Param(1).A;
Chain(1).Crel = SE3Inv(Chain(1).C);
Chain(1).AdCrel = SE3AdjMatrix(Chain(1).Crel);
Chain(1).V = Param(1).X*qd(1);
Chain(1).adX = SE3adMatrix(Param(1).X);
Chain(1).Vd = Param(1).X*q2d(1) - qd(1)*Chain(1).adX*Chain(1).V; %+ SE3AdjInvMatrix(Chain(1).C)*[[0;0;0];g];
Chain(1).V2d = Param(1).X*q3d(1) - q2d(1)*Chain(1).adX*Chain(1).V ...
                - 2*qd(1)*Chain(1).adX*Chain(1).Vd ...
                - qd(1)*qd(1)*Chain(1).adX*Chain(1).adX*Chain(1).V;
Chain(1).V3d = q4d(1)*Param(1).X ...
        - Chain(1).adX*(q3d(1)*Chain(1).V + 3*q2d(1)*Chain(1).Vd + 2*qd(1)*Chain(1).V2d ... 
                        + Chain(1).adX*(qd(1)*qd(1)*Chain(1).Vd + 2*qd(1)*q2d(1)*Chain(1).V));

% the forward recursion run
for i=2:n
    Chain(i).f = Chain(i-1).f*SE3Exp(Param(i).Y,q(i));
    Chain(i).C = Chain(i).f*Param(i).A;
    Chain(i).Crel = SE3Inv(Chain(i).C)*Chain(i-1).C;
    Chain(i).AdCrel = SE3AdjMatrix(Chain(i).Crel);
    Chain(i).V = Chain(i).AdCrel*Chain(i-1).V + qd(i)*Param(i).X;
    Chain(i).adX = SE3adMatrix(Param(i).X);
    Chain(i).Vd = Chain(i).AdCrel*Chain(i-1).Vd - qd(i)*Chain(i).adX*Chain(i).V + q2d(i)*Param(i).X;
    Chain(i).V2d = Chain(i).AdCrel*Chain(i-1).V2d + Param(i).X*q3d(i) ...
                    - q2d(i)*Chain(i).adX*Chain(i).V ...
                    - 2*qd(i)*Chain(i).adX*Chain(i).Vd - qd(i)*qd(i)*Chain(i).adX*Chain(i).adX*Chain(i).V;
    Chain(i).V3d = Chain(i).AdCrel*Chain(i-1).V3d + q4d(i)*Param(i).X ...
        - qd(i)*Chain(i).adX*Chain(i).AdCrel*Chain(i-1).V2d ...
        - Chain(i).adX*(q3d(i)*Chain(i).V + 3*q2d(i)*Chain(i).Vd + 2*qd(i)*Chain(i).V2d ... 
                        + Chain(i).adX*(qd(i)*qd(i)*Chain(i).Vd + 2*qd(i)*q2d(i)*Chain(i).V));
end

%% Inverse dynamics recursion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the backward recursion run
for i=n:-1:1
    AdInvC = SE3AdjInvMatrix(Chain(i).C);
    adV = SE3adMatrix(Chain(i).V);
    adVd = SE3adMatrix(Chain(i).Vd);
    adV2d = SE3adMatrix(Chain(i).V2d);
    adV3d = SE3adMatrix(Chain(i).V3d);

    if i == n
        Chain(i).W = Param(i).Mb*Chain(i).Vd - SE3adMatrix(Chain(i).V)'*Param(i).Mb*Chain(i).V;
        Chain(i).Wd = Param(i).Mb*Chain(i).V2d - SE3adMatrix(Chain(i).V)'*Param(i).Mb*Chain(i).Vd - SE3adMatrix(Chain(i).Vd)'*Param(i).Mb*Chain(i).V;
        Chain(i).W2d = Param(i).Mb*Chain(i).V3d - SE3adMatrix(Chain(i).V)'*Param(i).Mb*Chain(i).V2d ...
                        - SE3adMatrix(Chain(i).V2d)'*Param(i).Mb*Chain(i).V ...
                        - 2*SE3adMatrix(Chain(i).Vd)'*Param(i).Mb*Chain(i).Vd;
    else
        Chain(i).W = SE3AdjMatrix(Chain(i+1).Crel)'*Chain(i+1).W + Param(i).Mb*Chain(i).Vd - SE3adMatrix(Chain(i).V)'*Param(i).Mb*Chain(i).V;
        Chain(i).Wd = SE3AdjMatrix(Chain(i+1).Crel)'*(Chain(i+1).Wd - qd(i+1)*SE3adMatrix(Param(i+1).X)'*Chain(i+1).W) ...
                        + Param(i).Mb*Chain(i).V2d - SE3adMatrix(Chain(i).V)'*Param(i).Mb*Chain(i).Vd ...
                        - SE3adMatrix(Chain(i).Vd)'*Param(i).Mb*Chain(i).V;
        Chain(i).W2d = SE3AdjMatrix(Chain(i+1).Crel)'*(Chain(i+1).W2d ...
                        - 2*qd(i+1)*SE3adMatrix(Param(i+1).X)'*Chain(i+1).Wd ...
                        + (qd(i+1)*qd(i+1)*SE3adMatrix(Param(i+1).X)'*SE3adMatrix(Param(i+1).X)' ...
                        - q2d(i+1)*SE3adMatrix(Param(i+1).X)')*Chain(i+1).W) ...
                        + Param(i).Mb*Chain(i).V3d - SE3adMatrix(Chain(i).V)'*Param(i).Mb*Chain(i).V2d ...
                        - SE3adMatrix(Chain(i).V2d)'*Param(i).Mb*Chain(i).V - 2*SE3adMatrix(Chain(i).Vd)'*Param(i).Mb*Chain(i).Vd;
    end
    Chain(i).Q = Param(i).X'*Chain(i).W;
    Chain(i).Qd = Param(i).X'*Chain(i).Wd;
    Chain(i).Q2d = Param(i).X'*Chain(i).W2d;
    Q(i) = Chain(i).Q;
    Qd(i) = Chain(i).Qd;
    Q2d(i) = Chain(i).Q2d;
end

end
