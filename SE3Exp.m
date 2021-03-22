function [C] = SE3Exp(XX,t)
X = XX';
xi = X(1:3);
eta = X(4:6);
xihat = [0,    -X(3),  X(2);
     X(3),    0,  -X(1);
    -X(2), X(1),     0];
R = eye(3) + sin(t)*xihat + (1-cos(t))*(xihat*xihat);
p = (eye(3)-R)*(xihat*eta') + xi'*(xi*eta')*t;
[C] = [R,p;[0,0,0],[1]];
end

