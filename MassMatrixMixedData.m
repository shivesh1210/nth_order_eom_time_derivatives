function [M] = MassMatrixMixedData(m,Theta,COM)
[M] = [Theta(1,1),Theta(1,2),Theta(1,3),0,(-COM(3))*m,COM(2)*m;
Theta(1,2),Theta(2,2),Theta(2,3),COM(3)*m,0,(-COM(1))*m;
Theta(1,3),Theta(2,3),Theta(3,3),(-COM(2))*m,COM(1)*m,0;
0,COM(3)*m,(-COM(2))*m,m,0,0;
(-COM(3))*m,0,COM(1)*m,0,m,0;
COM(2)*m,(-COM(1))*m,0,0,0,m];
end

