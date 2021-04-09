%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author: Shivesh Kumar
% last modified: October 26, 2020
% remark: this code is by no means optimized

function [Qnd] = NthOrder_ClosedFormInvDyn_BodyFixed(qnd)

dim = size(qnd);
m = dim(2); % Order of derivatives required

q = qnd(:,1);
qd = qnd(:,2);
q2d = qnd(:,3);
q3d = qnd(:,4);

global Chain;
global Param;

global n; % DOF, number of joints
global g; % gravity vector

% Initialization for the first body
FK(1).f = SE3Exp(Param(1).Y,q(1));
FK(1).C = FK(1).f*Param(1).A;
% Compute FK for each body
for i = 2:n
    FK(i).f = FK(i-1).f*SE3Exp(Param(i).Y,q(i));
    FK(i).C = FK(i).f*Param(i).A;    
end
% Block diagonal matrix A (6n x 6n) of the Adjoint of body frame
A = []; 
for i=1:n
    A = blkdiag(A,eye(6,6));
    for j=1:i-1
        Crel = SE3Inv(FK(i).C)*FK(j).C;
        AdCrel = SE3AdjMatrix(Crel); 
        r = 6*(i-1)+1;
        c = 6*(j-1)+1;
        A(r:r+5,c:c+5) = AdCrel;       
    end
end

% Block diagonal matrix X (6n x n) of the screw coordinate vector associated to all joints in the body frame (Constant)
X = []; 
for i=1:n
X = blkdiag(X,Param(i).X);
end

% System level Jacobian 
J = A*X;

% System twist (6n x 1)
V = J*qd;

% Block diagonal matrix a (6n x 6n)
a = []; 
for i=1:n
% ad = SE3adMatrix(Param(i).X); 
a = blkdiag(a,SE3adMatrix(Param(i).X)*qd(i));
end

% System acceleration (6n x 1)
Ad = A*a - A*a*A;

% System acceleration (6n x 1)
Vd = J*q2d - A*a*V;
% Vd'
% Block Diagonal Mb (6n x 6n) Mass inertia matrix in body frame (Constant)
Mb = []; 
for i=1:n
Mb = blkdiag(Mb,Param(i).Mb);
end

% Block diagonal matrix b (6n x 6n) used in Coriolis matrix
b = []; 
for i=1:n
b = blkdiag(b,SE3adMatrix(V(6*i-5:6*i)));
end

% Block diagonal matrix Cb (6n x 6n)
Cb = -Mb*A*a - b'*Mb;

% Lets setup the Equations of Motion

% Mass inertia matrix in joint space (n x n)
M = J'*Mb*J;

% Coriolis-Centrifugal matrix in joint space (n x n)
C = J'*Cb*J;

% Gravity Term 
Utemp = zeros(6*n,6);
Utemp(1:6,1:6) = eye(6);
U = A*Utemp;
Vd_0 = zeros(6,1);
Vd_0(4:6) = g;
Qgrav = J'*Mb*U*Vd_0;

% Generalized forces Q
Q = M*q2d + C*qd + Qgrav;     % with gravity

%%%%%%%%%% For verification purpose only %%%%%%%%%%%%%%%%%%%%%%%%
% First time derivative of Block diagonal matrix a (6n x 6n)
ad = []; 
for i=1:n
ad = blkdiag(ad,q2d(i)*SE3adMatrix(Param(i).X));
end
% First time derivative of Block diagonal matrix b (6n x 6n) used in Coriolis matrix
bd = []; 
for i=1:n
bd = blkdiag(bd,SE3adMatrix(Vd(6*i-5:6*i)));
end
% Second time derivative of Block diagonal matrix a (6n x 6n)
a2d = []; 
for i=1:n
a2d = blkdiag(a2d,SE3adMatrix(Param(i).X)*q3d(i));
end
% Second time derivative of Block diagonal matrix b (6n x 6n) used in Coriolis matrix
% Third order Forward Kinematics
V2d = J*q3d - A*ad*V - 2*A*a*Vd - A*a*a*V; 
b2d = []; 
for i=1:n
b2d = blkdiag(b2d,SE3adMatrix(V2d(6*i-5:6*i)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Nth Order Derivatives of EOM %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% nth order derivative of matrix a (6n x 6n)
% a(qd), ad(qdd), a2d(q2d) etc.
for i = 1 : m - 1
    atemp = [];
    for j = 1 : n
        atemp = blkdiag(atemp,SE3adMatrix(Param(j).X)*qnd(j,i+1));
    end
    and(:,:,i) = atemp;
end

% nth order derivative of matrix A
% Init nth order derivative of A
And(:,:,1) = A;
And(:,:,2) = Ad;
for i = 2 : m -1
    first_term = zeros(6*n,6*n);
    second_term = zeros(6*n,6*n);
    im = i + 1; % matlab index i
    % Compute the first term
    for k = 0 : i-1
        order_A = i-1-k;
        order_a = k;
        first_term = first_term + nchoosek(i-1,k)*And(:,:,order_A+1)*and(:,:,order_a+1); 
    end
    
    % Compute 2nd term
    for k = 0 : i-1
        third_term = zeros(6*n,6*n);
        for j = 0 : k
            order_a = k - j;
            order_A = j;
            third_term = third_term + nchoosek(k,j)*and(:,:,order_a+1)*And(:,:,order_A+1);
        end
        order_A = i-1-k;
        second_term = second_term + nchoosek(i-1,k)*And(:,:,order_A+1)*third_term; 
    end
    
    And(:,:,im) = first_term - second_term;
end

% nth order derivative of matrix system Jacobian J
Jnd = zeros(6*n, n, m);
for i = 1 : m
    Jnd(:,:,i) = And(:,:,i)*X;
end

% nth order derivative of system Velocity V (first index stores V and subsequently the rest of derivatives - verified numerically)
Vnd = zeros(6*n, m-1);
% Vnd(:,1) = V;
for i = 0 : m-2     
    V_temp = zeros(6*n,1);
    im = i + 1; % matlab index
    for k = 0 : i
        order_J = i-k;
        order_qnd = k;
        V_temp = V_temp + nchoosek(i,k)*Jnd(:,:,order_J+1)*qnd(:,order_qnd+2); 
        % qnd_order is shifted by 1 for matlab indexing and 1 for the fact V = J*qd starts with first order derivative
    end
    Vnd(:,im) = V_temp;
end

% implementation of mass matrix derivative
Mnd = zeros(n,n,m); % projected to joint space
for i = 0:m-1
    term = zeros(n,n);
    for k = 0:i
        order_Jnd = i - k;
        order_JndT = k;
        term = term + nchoosek(i,k)*Jnd(:,:,order_Jnd+1)'*Mb*Jnd(:,:,order_JndT + 1);
    end
    Mnd(:,:,i+1) = term;
end

% Only for verification
% % First time derivative of Mass inertia matrix in joint space (n x n)
% Mbd = -Mb*A*a - (Mb*A*a)';
% Md = J'*Mbd*J;
% % Second time derivative of Mass inertia matrix in joint space (n x n)
% Mb2d = - Mb*A*ad - (Mb*A*ad)' + 2*Mb*A*a*A*a + 2*(Mb*A*a*A*a)' + 2*a'*A'*Mb*A*a - Mb*A*a*a - (Mb*A*a*a)';   
% M2d = J'*Mb2d*J;
% Mnd(:,:,1) - M
% Mnd(:,:,2) - Md
% Mnd(:,:,3) - M2d


% nth order derivative of matrix b (6n x 6n) used in Coriolis matrix
% b(V), bd(Vd), b2d(V2d), and so on.
bnd = []; 
for i = 1 : m - 1
    btemp = [];
    for j = 1 : n
        btemp = blkdiag(btemp,SE3adMatrix(Vnd(6*j-5:6*j,i)));
    end
    bnd(:,:,i) = btemp;
end

%%%%%%%%%%%%%%%%%%%%%%% Coriolis-Centrifugal Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% implementation of coriolis-centrifugal matrix
Cpnd = zeros(6*n,6*n,m-1); % projected to joint space
for i = 0:m-2
    term = zeros(6*n,6*n);
    for k = 0:i
        order_And = i - k;
        order_and = k;
        term = term + nchoosek(i,k)*And(:,:,order_And+1)*and(:,:,order_and + 1);
    end
    Cpnd(:,:,i+1) = - Mb*term - bnd(:,:,i+1)'*Mb;
end

Cnd = zeros(n,n,m-1); % projected to joint space
for i = 0:m-2
    temp1 = zeros(n,n);
    for k = 0:i
        temp2 = zeros(6*n,n);
        for j = 0:k
            order_Cpnd = k-j;
            order_Jnd_inner = j;
            temp2 = temp2 + nchoosek(k,j)*Cpnd(:,:,order_Cpnd + 1)*Jnd(:,:,order_Jnd_inner + 1);
        end
        order_Jnd_outer = i - k; 
        temp1 = temp1 + nchoosek(i,k)*Jnd(:,:,order_Jnd_outer + 1)'*temp2; 
    end
    Cnd(:,:,i+1) = temp1;
end

% only for verification
% % First time derivative of Coriolis-Centrifugal matrix in joint space (n x n)
% Cbd = Mb*A*a*A*a - Mb*A*a*a - Mb*A*ad - bd'*Mb - Cb*A*a - a'*A'*Cb; 
% Cd = J'*Cbd*J;
% Cnd(:,:,2) - Cd
% % Second time derivative of Coriolis-Centrifugal matrix in joint space (n x n)
% Cddot = - Mb*A*a2d - 3*Mb*A*a*ad - Mb*A*a*a*a - b2d'*Mb + Mb*A*ad*A*a + ...
% Mb*A*a*a*A*a + 2*Mb*A*a*A*ad + 2*Mb*A*a*A*a*a - 2*Mb*A*a*A*a*A*a;
% Cb2d = Cddot - (Cbd + a'*A'*Cb)*A*a - a'*A'*(Cbd + Cb*A*a) - Cb*A*ad - ...
% ad'*A'*Cb - Cb*A*a*a - a'*a'*A'*Cb - Cbd*A*a - a'*A'*Cbd; 
% C2d = J'*Cb2d*J;
% Cnd(:,:,3) - C2d

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%% Gravity Forces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% implementation of gravity term (verified until any order)
Qgravnd = zeros(n, m);
for i = 0 : m - 1
    temp = zeros(n,6*n);
    for k = 0:i
        order_Jnd = i - k;
        order_And = k;
        temp = temp + nchoosek(i,k)*Jnd(:,:,order_Jnd + 1)'*Mb*And(:,:,order_And + 1);
    end
    Qgravnd(:,i+1) = temp*Utemp*Vd_0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% nth order derivative of generalized forces 
Qnd(:,1) = Q;   % init
for i = 1:m-3
    im = i+1;
    term = zeros(n,1);
    for k = 2:i+1
        order_Mnd = i - k + 2;
        order_Cnd = i - k + 1;
        Pk = nchoosek(i,k-2)*Mnd(:,:,order_Mnd+1) + nchoosek(i,k-1)*Cnd(:,:,order_Cnd+1);
        term = term + Pk*qnd(:,k+1);
    end
    Qnd(:,im) = Mnd(:,:,1)*qnd(:,im+2) + term + Cnd(:,:,im)*qnd(:,2) + Qgravnd(:,im);
end
