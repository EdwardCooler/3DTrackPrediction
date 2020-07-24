%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 用EKF产生建议分布
%  输入参数说明
%  Xiset：上t-1时刻的粒子集合
%  Z：t时刻的观测
%  Pin：Xiset粒子集合的方差集合

% 输出参数说明
% Xo：EPF算法最终的估计结果
% Xoset：k时刻的粒子结合，其均值就是Xo
% Pout：Xoset对应的方差
function [Xo,Xoset,Pout,Neff]=epf(Xiset,Z,n,Pin,N,R,Qekf,Rekf,Station)
% 重采样策略函数
resamplingScheme=1;
% 中间变量初始化
Zpre=ones(3,N);    % 观测预测     
Xsetpre=ones(n,N); % 粒子集合预测
w = ones(1,N);     % 权值初始化
Xo=zeros(n,1);

Pout=ones(n,n*N);    % 协方差预测，这个Pout也是一个9x(9xN)的矩阵，需要修改
% Xekf=ones(n,1);    % EKF估计结果
% Xekf_pre=ones(1,N);% EKF的一步预测

% 第一步：根据EKF计算得到的结果进行采样
for i=1:N
    % 利用EKF计算得到每个粒子的均值和协方差方差，注意这里放进来的粒子集合Xiset是(9xN)的，协方差Pin是9x(9xN)的，后者需分开为N个9x9的矩阵，怎么分开呢？
    % Pin(:,9.*(i-1)+1:9.*(i-1)+9)
    [Xekf,Pout(:,9*(i-1)+1:9*(i-1)+9)]=ekf(Xiset(:,i),Z,Pin(:,9*(i-1)+1:9*(i-1)+9),Qekf,Rekf,Station);
    %[Xekf,Pout]=ekf(Xin,Z,Pin,Qekf,Rekf,Station)
    % 用上面得到的均值和方差来为粒子集合采样
    Xsetpre(:,i)=Xekf + sqrtm(Pout(:,9*(i-1)+1:9*(i-1)+9)) * randn(n,1);
end

% 第二步：计算权重
for i=1:N
    % 观测预测
    [dd,alpha,beta]=feval('hfun',Xsetpre(:,i),Station);
    Zpre(:,i) =[dd,alpha,beta]';
    z1 = Z-Zpre(:,i);
    % 计算权重，1e-99为最小非0数字，防止变0
    lik = inv(sqrt(2*pi*det(R)))*exp(-.5*(z1)'*inv(R)*(z1))+ 1e-99;
    % 先验
    %prior = ((Xsetpre(i)-Xiset(i))^(g1-1)) * exp(-g2*(Xsetpre(i)-Xiset(i)));
    % 建议分布
    %proposal = inv(sqrt(Pout(i))) * exp(-0.5*inv(Pout(i)) *((Xsetpre(i)-Xekf(i))^(2)));
    %w(i) = lik*prior/proposal;
    w(i) = lik;
end
% 权值归一化
w= w./sum(w);

summ = 0;
alfa = 0.2;
for i=1:N
    summ = summ + w(i)^alfa; 
end
for i=1:N
    w(i) =  w(i)^alfa / summ; 
end

% 求有效粒子数
Neff = 1/sum(w.^2);

% 第三步：重采样
if resamplingScheme == 1
    outIndex = residualR(1:N,w');   
elseif resamplingScheme == 2
    outIndex = systematicR(1:N,w');     
else
    outIndex = multinomialR(1:N,w');    
end

% 第四步，更新粒子集合 
Xoset = Xsetpre(:,outIndex); 
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


