%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author: Shivesh Kumar
% last modified: October 26, 2020
% remark: this code is by no means optimized

function [Qnd] = NthOrder_RecurisveInvDyn_BodyFixed(qnd)

dim = size(qnd);
m = dim(2); % Order of derivatives required

q = qnd(:,1);
qd = qnd(:,2);
q2d = qnd(:,3);
q3d = qnd(:,4);
q4d = qnd(:,5);

global Chain;
global Param;

global n; % DOF, number of joints
global g; % gravity vector

%% Forward kinematics recursion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note the last dimension is used for storing the order of the derivative.

% first body is treated separately since there is no predecessor
Chain(1).f(:,:,1) = SE3Exp(Param(1).Y,q(1));
Chain(1).C(:,:,1) = Chain(1).f*Param(1).A;
Chain(1).Crel(:,:,1) = SE3Inv(Chain(1).C);
Chain(1).AdCrel(:,:,1) = SE3AdjMatrix(Chain(1).Crel);
Chain(1).DAdCrel(:,:,1) = Chain(1).AdCrel(:,:,1);

Chain(1).V(:,1) = Param(1).X*qd(1);
Chain(1).adX = SE3adMatrix(Param(1).X);

Chain(1).V(:,2) = Param(1).X*q2d(1); % this should go into the (for i=2:m ...) loop below 
Chain(1).V(:,3) = Param(1).X*q3d(1);
Chain(1).V(:,4) = Param(1).X*q4d(1);

Chain(1).B(:,1,1) = SE3AdjMatrix(Chain(1).Crel)*Param(1).X; % 
Chain(1).Bsf(:,1,1) = Chain(1).B(:,1,1)*qd(1);  % eqn 30
for i=2:m
    Chain(1).B(:,1,i) = zeros(6,1); % Lie bracket zero
end

% Compute the basic kinematic quantities
for i=2:n   % loop through bodies
    Chain(i).f = Chain(i-1).f*SE3Exp(Param(i).Y,q(i));
    Chain(i).C = Chain(i).f*Param(i).A;
    for j = 1:i
        Chain(i).Crel(:,:,j) = SE3Inv(Chain(i).C)*Chain(j).C;
        Chain(i).AdCrel(:,:,j) = SE3AdjMatrix(Chain(i).Crel(:,:,j));
        Chain(i).B(:,j,1) = SE3AdjMatrix(Chain(i).Crel(:,:,j))*Param(j).X;
    end
    Chain(i).DAdCrel(:,:,1) = Chain(i).AdCrel(:,:,i-1);
end

for r = 2:m-1   % loop through derivatives
    for i=1:n   % loop through bodies
        temp3 = zeros(6,6);
        real_r = r - 1;
        for p = 0:real_r-1
            order_qnd = real_r - p;
            temp3 = temp3 + nchoosek(real_r - 1, p)*Chain(i).DAdCrel(:,:,p+1)*qnd(i,order_qnd+1);
        end
        Chain(i).DAdCrel(:,:,r) =  - SE3adMatrix(Param(i).X)*temp3;
    end
end

% rth-order forward kinematics run 
for r = 1:m-1   % loop through derivatives
    for i=2:n   % loop through bodies
        real_r = r - 1;
        for j = i:-1:1
            temp = zeros(6,1);
            temp2 = zeros(6,1);
            if r > 1  % actual derivative
                for l = 0 : real_r - 1
                    order_B = l;
                    order_Bsf = real_r - l - 1;
                    temp = temp + nchoosek(real_r-1,l)*SE3adMatrix(Chain(i).B(:,j,order_B+1))*Chain(i).Bsf(:,j,order_Bsf + 1);
                end
                Chain(i).B(:,j,r) = temp;
            end
            for k = j : i
                for l = 0 : real_r
                    order_B = l;
                    order_q = real_r-l+1;
                    temp2 = temp2 + nchoosek(real_r,l)*Chain(i).B(:,k,order_B + 1)*qnd(k,order_q+1);
                end
                Chain(i).Bsf(:,j,r) = temp2;
            end
        end
        Chain(i).V(:,r) = Chain(i).Bsf(:,1,r);
    end
end

Vnd = [];
for i = 1:n
    Vnd = [Vnd; Chain(i).V];
end


%% Inverse Dynamics Recursion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1 : m - 2
    
    real_k = k - 1;
    temp = zeros(6,1);
    for p = 0:real_k
        temp = temp - nchoosek(real_k,p)*SE3adMatrix(Chain(n).V(:,p+1))'*Param(n).Mb*Chain(n).V(:,real_k-p+1);
    end
    Chain(n).W(:,k) = Param(n).Mb*Chain(n).V(:,k+1) + temp;
    Chain(n).Q(k) = Param(n).X'*Chain(n).W(:,k);

    for i = n-1:-1:1
        term = zeros(6,1);
        for p = 0:real_k
            term = term - nchoosek(real_k,p)*SE3adMatrix(Chain(i).V(:,p+1))'*Param(i).Mb*Chain(i).V(:,real_k-p+1) ...
            + nchoosek(real_k,p)*Chain(i+1).DAdCrel(:,:,p+1)'*Chain(i+1).W(:,real_k-p+1);
        end
        Chain(i).W(:,k) = Param(i).Mb*Chain(i).V(:,k+1) + term;
        Chain(i).Q(k) = Param(i).X'*Chain(i).W(:,k);
    end
end

Qnd = [];
for i = 1:n
    Qnd = [Qnd; Chain(i).Q];
end
