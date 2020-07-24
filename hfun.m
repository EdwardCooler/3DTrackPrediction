%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 观测方程函数
 function [dd,alpha,beta]=hfun(X1,Station)
if nargin < 2
    error('Not enough input arguments.');
end

dd = sqrt((X1(1,1)-Station(1,1))^2+(X1(2,1)-Station(2,1))^2+(X1(3,1)-Station(3,1))^2);
alpha = atan((X1(2,1)-Station(2,1))/(X1(1,1)-Station(1,1)));
beta = atan((X1(3,1)-Station(3,1))/sqrt((X1(1,1)-Station(1,1))^2+(X1(2,1)-Station(2,1))^2));

