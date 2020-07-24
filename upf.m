%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 用UKF改进的粒子滤波算法--EPF

% 用UKF产生建议分布
% 输入参数说明：
%    Xiset是上t-1时刻的粒子集合，Z是t时刻的观测
%    Pin对应Xiset粒子集合的方差
% 输出参数说明：
%    Xo是upf算法最终的估计结果
%    Xoset是k时刻的粒子集合，其均值就是Xo
%    Pout是Xoset对应的方差
function [Xo,Xoset,Pout]=upf(Xiset,Z,n,Pin,N,R,Qukf,Rukf,Station)
 
% 重采样策略函数
resamplingScheme=1;

% 中间变量初始化
Zpre=ones(3,N);      % 观测预测 
Xset_pre=ones(n,N);  % 粒子集合预测
w = ones(1,N);     % 权值初始化
Xo=zeros(n,1);

Pout=ones(n,n*N);    % 协方差预测，这个Pout也是一个9x(9xN)的矩阵，需要修改
% Xukf=ones(1,N);      % UKF估计结果 
% Xukf_pre=ones(1,N);  % EKF的一步预测

% 第一步：根据UKF计算得到的结果进行采样
for i=1:N
     % 利用UKF计算得到均值和方差，第i个粒子的均值和方差
    [Xukf,Pout(:,9*(i-1)+1:9*(i-1)+9)]=function_ukf(Station,Xiset(:,i),Pin(:,9*(i-1)+1:9*(i-1)+9),Z,Qukf,Rukf);
    % [X,P]=function_ukf(Station,X,P,Z,Q,R)
    % 用上面得到的均值和方差来为粒子集合采样
    Xset_pre(i) = Xukf + sqrtm(Pout(:,9*(i-1)+1:9*(i-1)+9)) * randn(n,1);
end

% 第二步：计算权重
for i=1:N
    % 观测预测
    [dd,alpha,beta]=feval('hfun',Xset_pre(:,i),Station);
    Zpre(:,i) =[dd,alpha,beta]';
    z1 = Z-Zpre(:,i);
    % 计算权重，1e-99为最小非0数字，防止变0
    lik = inv(sqrt(2*pi*det(R)))*exp(-.5*(z1)'*inv(R)*(z1))+ 1e-99;
    % prior = ((Xset_pre(i)-Xiset(i))^(g1-1)) * exp(-g2*(Xset_pre(i)-Xiset(i)));
    % proposal = inv(sqrt(Pout(i))) * exp(-0.5*inv(Pout(i)) *((Xset_pre(i)-Xukf(i))^(2)));
    % w(i) = lik*prior/proposal;
    w(i) = lik;
end

% 权值归一化 
w = w./sum(w);

% 第三步：重采样 
if resamplingScheme == 1
    outIndex = residualR(1:N,w');        
elseif resamplingScheme == 2
    outIndex = systematicR(1:N,w');      
else
    outIndex = multinomialR(1:N,w');     
end

% 第四步，更新粒子集合 
Xoset = Xset_pre(:,outIndex); 
% 更新粒子方差
Pm = ones(n,n*N);      % 协方差预测
for i = 1:N
    Pm(:,9*(i-1)+1:9*(i-1)+9) = Pout(:,9*(outIndex(i)-1)+1:9*(outIndex(i)-1)+9);
end
Pout = Pm;
% 第五步，得到本次计算的粒子滤波估计值
target=[mean(Xoset(1,:)),mean(Xoset(2,:)),mean(Xoset(3,:)),mean(Xoset(4,:)),mean(Xoset(5,:)),mean(Xoset(6,:)),mean(Xoset(7,:)),mean(Xoset(8,:)),mean(Xoset(9,:))]';
Xo(:,1)=target;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


