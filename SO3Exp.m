function [R] = SO3Exp(x,t)
xhat = [0,    -x(3),  x(2);
     x(3),    0,  -x(1);
    -x(2), x(1),     0];
[R] = eye(3) + sin(t)*xhat + (1-cos(t))*(xhat*xhat);
end

