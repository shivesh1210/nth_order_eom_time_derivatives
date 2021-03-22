%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author: Shivesh Kumar
% last modified: October 18, 2020
% remark: this code is by no means optimized

function [Q,Qd,Q2d] = ClosedFormInvDyn_BodyFixed(q,qd,q2d,q3d,q4d, WEE, WDEE, W2DEE)

if nargin == 5   % if the number of inputs equals 5
  WEE = zeros(6,1); 
  WDEE = zeros(6,1);
  W2DEE = zeros(6,1);
end

if nargin == 6   % if the number of inputs equals 6
  WDEE = zeros(6,1);
  W2DEE = zeros(6,1);
end

if nargin == 7   % if the number of inputs equals 7
  W2DEE = zeros(6,1);
end

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

% External Wrench
Wext = zeros(6*n,1);
Wext(end - 5:end) = WEE;    % WEE (t) is the time varying wrench on the EE link. 
Qext = J'*Wext;

% Generalized forces Q
% Q = M*q2d + C*qd;   % without gravity
Q = M*q2d + C*qd + Qgrav + Qext;     % with gravity

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%First Order Derivatives of EOM%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First time derivative of Block diagonal matrix a (6n x 6n)
ad = []; 
for i=1:n
ad = blkdiag(ad,q2d(i)*SE3adMatrix(Param(i).X));
end

% Third order Forward Kinematics
V2d = J*q3d - A*ad*V - 2*A*a*Vd - A*a*a*V; 

% First time derivative of Block diagonal matrix b (6n x 6n) used in Coriolis matrix
bd = []; 
for i=1:n
bd = blkdiag(bd,SE3adMatrix(Vd(6*i-5:6*i)));
end

% First time derivative of Mass inertia matrix in joint space (n x n)
Mbd = -Mb*A*a - (Mb*A*a)';
Md = J'*Mbd*J;

% First time derivative of Coriolis-Centrifugal matrix in joint space (n x n)
Cbd = Mb*A*a*A*a - Mb*A*a*a - Mb*A*ad - bd'*Mb - Cb*A*a - a'*A'*Cb; 
Cd = J'*Cbd*J;

% First time derivative of gravity force
Qdgrav = J'*Mbd*U*Vd_0;

% First time derivative of External Wrench
Wdext = zeros(6*n,1);
Wdext(end - 5:end) = WDEE;    % WDEE (t) is the time derivative of wrench on the EE link. 
Qdext = J'*(Wdext - (A*a)'*Wext);

% Qd = M*q3d + (Md + C)*q2d + Cd*qd;  % without gravity
Qd = M*q3d + (Md + C)*q2d + Cd*qd + Qdgrav + Qdext; % with gravity

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Second Order Derivatives of EOM%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Second time derivative of Block diagonal matrix a (6n x 6n)
a2d = []; 
for i=1:n
a2d = blkdiag(a2d,SE3adMatrix(Param(i).X)*q3d(i));
end

% Second time derivative of Block diagonal matrix b (6n x 6n) used in Coriolis matrix
b2d = []; 
for i=1:n
b2d = blkdiag(b2d,SE3adMatrix(V2d(6*i-5:6*i)));
end

% Second time derivative of Mass inertia matrix in joint space (n x n)
Mb2d = - Mb*A*ad - (Mb*A*ad)' + 2*Mb*A*a*A*a + 2*(Mb*A*a*A*a)' + 2*a'*A'*Mb*A*a - Mb*A*a*a - (Mb*A*a*a)';   
M2d = J'*Mb2d*J;

% Second time derivative of Coriolis-Centrifugal matrix in joint space (n x n)
Cddot = - Mb*A*a2d - 3*Mb*A*a*ad - Mb*A*a*a*a - b2d'*Mb + Mb*A*ad*A*a + ...
Mb*A*a*a*A*a + 2*Mb*A*a*A*ad + 2*Mb*A*a*A*a*a - 2*Mb*A*a*A*a*A*a;
Cb2d = Cddot - (Cbd + a'*A'*Cb)*A*a - a'*A'*(Cbd + Cb*A*a) - Cb*A*ad - ad'*A'*Cb - Cb*A*a*a - a'*a'*A'*Cb - Cbd*A*a - a'*A'*Cbd; 
C2d = J'*Cb2d*J;


% Cb2d = Mb*A*a2d - b2d'*Mb + 2*bd'*Mb*A*a + 2*a'*A'*bd'*Mb + 3*Mb*A*ad*A*a ...
% + 2*Mb*A*a*A*ad + 2*a'*A'*Mb*A*ad - Cb*A*ad - ad'*A'*Cb  + 2*Cb*A*a*A*a ...
% + 2*a'*A'*Cb*A*a + 2*a'*A'*a'*A'*Cb - 4*Mb*A*a*A*a*A*a - 2*a'*A'*Mb*A*a*A*a ... 
% + 2*Mb*A*a*a*A*a + 2*a'*A'*Mb*A*a*a + 2*Mb*A*a*A*a*a - 3*Mb*A*a*ad ...
% - Mb*A*a*a*a + Mb*A*a*a*A*a - Cb*A*a*a - a'*a'*A'*Cb;
% C2d = J'*Cb2d*J;

% Second time derivative of gravity force
Q2dgrav = J'*Mb2d*U*Vd_0;

% Second time derivative of External Wrench
W2dext = zeros(6*n,1);
W2dext(end - 5:end) = W2DEE;    % WDEE (t) is the time derivative of wrench on the EE link. 
Q2dext = J'*(W2dext - 2*(A*a)'*Wdext  + (2*(A*a*A*a)' - (A*ad)' - (A*a*a)')*Wext);

% Second time derivative of generalized forces
% Q2d = M*q4d + (2*Md + C)*q3d + (M2d + 2*Cd)*q2d + C2d*qd;   % without gravity
Q2d = M*q4d + (2*Md + C)*q3d + (M2d + 2*Cd)*q2d + C2d*qd + Q2dgrav + Q2dext;   % with gravity and external forces


