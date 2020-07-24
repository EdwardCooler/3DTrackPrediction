%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 产生Sigma点集
function Xset=function_sigmas(X,P,c)
A = c*chol(P)'; % Cholesky分解
Y = X(:,ones(1,numel(X))); % 这里应该是获取了n行n列矩阵，并且所有的列相等（用于后面加减法，秒啊）
Xset = [X Y+A Y-A];