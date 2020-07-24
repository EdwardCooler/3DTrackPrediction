%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 基本粒子滤波算法
% 输入：Xiset为二维数组
% 输出：Xo为nx1矩阵，Xoset为nxN矩阵
function [Xo,Xoset,Neff]=pf(Xiset,Z,N,n,R,Q,Station)
 
tic
% 重采样策略
resamplingScheme=1;

% 中间变量初始化
Zpre=ones(3,N);     % 观测预测   
Xsetpre=ones(n,N);  % 粒子集合预测
w = ones(1,N);      % 权值初始化
Xo=zeros(n,1);

% 第一步，根据每一个粒子对先验分布采样 
for i=1:N
    [y1,y2,y3,y4,y5,y6,y7,y8,y9] = ffun(Xiset(:,i));
    Xsetpre(:,i) = [y1,y2,y3,y4,y5,y6,y7,y8,y9]' + 10*sqrtm(Q)*randn(n,1);
end

% 第二步，计算粒子权重
for i=1:N
    [dd,alpha,beta]=feval('hfun',Xsetpre(:,i),Station);
    Zpre(:,i) =[dd,alpha,beta]';
    z1 = Z-Zpre(:,i);
    % w(i) = inv(sqrtm(R)) * exp(-0.5*inv(R)*((Z-Zpre(:,i))^(2))) + 1e-99; 
    w(i) = inv(sqrt(2*pi*det(R)))*exp(-.5*(z1)'*inv(R)*(z1))+ 1e-99;%权值计算，这里其实可以把inv(sqrt(2*pi*det(R)))去掉
end

% 第三步，权重归一化
w = w./sum(w);   

% 求有效粒子数
Neff = 1/sum(w.^2);

% 将随机采样之后的图给画出来

% 第四步，重采样
if resamplingScheme == 1
    outIndex = residualR(1:N,w');       
elseif resamplingScheme == 2
    outIndex = systematicR(1:N,w');  
else
    outIndex = multinomialR(1:N,w');  
end

% 第五步，更新粒子集合 
Xoset = Xsetpre(:,outIndex); 
% 第六步，得到本次计算的粒子滤波估计值
target=[mean(Xoset(1,:)),mean(Xoset(2,:)),mean(Xoset(3,:)),mean(Xoset(4,:)),mean(Xoset(5,:)),mean(Xoset(6,:)),mean(Xoset(7,:)),mean(Xoset(8,:)),mean(Xoset(9,:))]';
Xo(:,1)=target;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


