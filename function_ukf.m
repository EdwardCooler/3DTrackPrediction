function [X,P]=function_ukf(Station,X,P,Z,Q,R)
L=numel(X);                           % 状态维数（行数）
m=numel(Z);                           % 观测维数（行数）
alpha=1e-2;                           % 默认系数，参照UT变换，下同
ki=0;                                 % 默认系数
beta=2;                               % 默认系数
lambda=alpha^2*(L+ki)-L;              % 默认系数
c=L+lambda;                           % 默认系数
Wm=[lambda/c 0.5/c+zeros(1,2*L)];     % 权值
Wc=Wm;
Wc(1)=Wc(1)+(1-alpha^2+beta);         % 权值
c=sqrt(c);
% 第一步，获得一组Sigma点集
% Sigma点集，在状态X附近的点集，X是6*13矩阵，每列为1样本
Xsigmaset=function_sigmas(X,P,c);              
% 第二、三、四步：对Sigma点集进行一步预测，得到均值X1means和协方差P1和新sigma点集X1
[X1means,X1,P1,X2]=function_ut(Station,1,Xsigmaset,Wm,Wc,L,Q);   
% 第五、六步：得到观测预测，Z1为X1集合的预测，Zpre为Z1的均值，Pzz为协方差
[Zpre,Z1,Pzz,Z2]=function_ut(Station,0,X1,Wm,Wc,m,R);
% 第七步，计算Pxz
Pxz=X2*diag(Wc)*Z2';
% 计算卡尔曼增益
K=Pxz*inv(Pzz);
% 第八步，状态和协方差更新
X=X1means+K*(Z-Zpre);       
P=P1-K*Pxz';