%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ״̬���̺���
function [y1,y2,y3,y4,y5,y6,y7,y8,y9] = ffun(x)
 
if nargin < 1
    error('NoT enough inpuT argumenTs.'); 
end
T=0.5; 

y1 = x(1)+T*x(4)+0.5*T^2*x(7);
y2 = x(2)+T*x(5)+0.5*T^2*x(8);
y3 = x(3)+T*x(6)+0.5*T^2*x(9);
y4 = x(4)+T*x(7);
y5 = x(5)+T*x(8);
y6 = x(6)+T*x(9);
y7 = x(7);
y8 = x(8);
y9 = x(9);