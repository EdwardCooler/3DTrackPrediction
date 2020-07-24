% 功能说明：ekf,ukf,pf,改进pf算法的无人机航迹预测比较程序

function main
% 因本程序涉及太多的随机数，下面让随机数每次都不变
rand('seed',3);
randn('seed',6);
% error('下面的参数T请参考书中的值设置，然后删除本行代码') 
n = 9;
T = 50;

Q= [1 0 0 0 0 0 0 0 0;    % 过程噪声协方差矩阵
    0 1 0 0 0 0 0 0 0;
    0 0 1 0 0 0 0 0 0;
    0 0 0 0.01 0 0 0 0 0;
    0 0 0 0 0.01 0 0 0 0;
    0 0 0 0 0 0.01 0 0 0;
    0 0 0 0 0 0 0.0001 0 0;
    0 0 0 0 0 0 0 0.0001 0;
    0 0 0 0 0 0 0 0 0.0001];

R = [5000 0 0;                 % 观测噪声协方差矩阵  
    0 0.01^2 0                   % 角度的观测值偏差不能给的太大
    0 0 0.01^2];   

% 系统初始化
X = zeros(9,T);  % 真实值
Z = zeros(3,T);
  
% 真实状态初始化
%X(:,1)=[1000;5000;200;10;50;10;2;-4;2]+sqrtm(Q)*randn(n,1);
X(:,1)=[100;500;20;10;50;10;2;-4;2]+sqrtm(Q)*randn(n,1);
state0 = X(:,1);

x0=0;           
y0=0; 
z0=0;
Station=[x0;y0;z0];   % 观测站的位置


P0 =[100 0 0 0 0 0 0 0 0;             % 协方差初始化
     0 100 0 0 0 0 0 0 0;
     0 0 100 0 0 0 0 0 0;
     0 0 0 1 0 0 0 0 0;
     0 0 0 0 1 0 0 0 0;
     0 0 0 0 0 1 0 0 0;
     0 0 0 0 0 0 0.1 0 0;
     0 0 0 0 0 0 0 0.1 0
     0 0 0 0 0 0 0 0 0.1];
 
%%%%%%%%%%%%% EKF滤波算法 %%%%%%%%%%%%
Qekf = Q;                   % EKF过程噪声方差
Rekf = R;                   % EKF过程噪声方差 
Xekf=zeros(9,T);            % 滤波状态
Xekf(:,1)=X(:,1);           % EKF滤波初始化
Pekf = P0;                  % 协方差
Tekf=zeros(1,T);            % 用于记录一个采样周期的算法时间消耗
%%%%%%%%%%%%% UKF滤波算法 %%%%%%%%%%%%           
Qukf = Q;              % UKF过程噪声方差
Rukf = R;              % UKF观测噪声方差 
Xukf=zeros(9,T);       % 滤波状态
Xukf(:,1)=X(:,1);      % UKF滤波初始化
Pukf = P0;             % 协方差
Tukf=zeros(1,T);       % 用于记录一个采样周期的算法时间消耗 
%%%%%%%%%%%%% PF滤波算法 %%%%%%%%%%%%
N = 200;                 % 粒子数              
Xpf=zeros(n,T);        % 滤波状态
Xpf(:,1)=X(:,1);       % PF滤波初始化
Xpfset=ones(n,N);      % 粒子集合初始化
for j=1:N   % 粒子集初始化
    Xpfset(:,j)=state0;  % 全都初始化为x0，每个粒子的值相等
end
Tpf=zeros(1,T);        % 用于记录一个采样周期的算法时间消耗 

%%%%%%%%%%%%% PF2滤波算法 %%%%%%%%%%%%
N2 = 200;                 % 粒子数 
R2 = [5000 0 0;                 % 观测噪声协方差矩阵  
    0 0.01^2 0                   % 角度的观测值偏差不能给的太大
    0 0 0.01^2];
Xpf2=zeros(n,T);        % 滤波状态
Xpf2(:,1)=X(:,1);       % PF滤波初始化
Xpf2set=ones(n,N2);      % 粒子集合初始化
for j=1:N2   % 粒子集初始化
    Xpf2set(:,j)=state0;  % 全都初始化为x0，每个粒子的值相等
end

%%%%%%%%%%%%% PF3滤波算法 %%%%%%%%%%%%
N3 = 400;                 % 粒子数 
R3 = [5000 0 0;                 % 观测噪声协方差矩阵  
    0 0.01^2 0                   % 角度的观测值偏差不能给的太大
    0 0 0.01^2];
Xpf3=zeros(n,T);        % 滤波状态
Xpf3(:,1)=X(:,1);       % PF滤波初始化
Xpf3set=ones(n,N3);      % 粒子集合初始化
for j=1:N3   % 粒子集初始化
    Xpf3set(:,j)=state0;  % 全都初始化为x0，每个粒子的值相等
end

%%%%%%%%%%%%% EPF滤波算法 %%%%%%%%%%%% 
Xepf=zeros(9,T);       % 滤波状态
Xepf(:,1)=X(:,1);      % EPF滤波初始化
Xepfset=ones(n,N);     % 粒子集合初始化，这里需要定义为一个3维数组，或者简单起见，一次性写完，定义为一个(9xN)的二维数组，表示当前状态的所有粒子
for j=1:N   % 粒子集初始化
    Xepfset(:,j)=state0;  % 全都初始化为state0，每个粒子的值相等
end
Pepf = P0*ones(n,n*N);% 各个粒子的协方差，这里需要定义为一个3维数组，或者简单起见，一次性写完，定义为一个9x(9xN)的二维数组
Tepf=zeros(1,T);       % 用于记录一个采样周期的算法时间消耗   

%%%%%%%%%%%%% EPF2滤波算法 %%%%%%%%%%%%
Rekf2 = R2;
Xepf2=zeros(9,T);       % 滤波状态
Xepf2(:,1)=X(:,1);      % EPF滤波初始化
Xepf2set=ones(n,N2);     % 粒子集合初始化，这里需要定义为一个3维数组，或者简单起见，一次性写完，定义为一个(9xN)的二维数组，表示当前状态的所有粒子
for j=1:N2   % 粒子集初始化
    Xepf2set(:,j)=state0;  % 全都初始化为state0，每个粒子的值相等
end
Pepf2 = P0*ones(n,n*N2);% 各个粒子的协方差，这里需要定义为一个3维数组，或者简单起见，一次性写完，定义为一个9x(9xN)的二维数组

%%%%%%%%%%%%% EPF3滤波算法 %%%%%%%%%%%%
Rekf3 = R3;
Xepf3=zeros(9,T);       % 滤波状态
Xepf3(:,1)=X(:,1);      % EPF滤波初始化
Xepf3set=ones(n,N3);     % 粒子集合初始化，这里需要定义为一个3维数组，或者简单起见，一次性写完，定义为一个(9xN)的二维数组，表示当前状态的所有粒子
for j=1:N3   % 粒子集初始化
    Xepf3set(:,j)=state0;  % 全都初始化为state0，每个粒子的值相等
end
Pepf3 = P0*ones(n,n*N3);% 各个粒子的协方差，这里需要定义为一个3维数组，或者简单起见，一次性写完，定义为一个9x(9xN)的二维数组

%%%%%%%%%%%%% UPF滤波算法 %%%%%%%%%%%%  
Xupf=zeros(n,T);       % 滤波状态  
Xupf(:,1)=X(:,1);      % UPF滤波初始化
Xupfset=ones(n,N);     % 粒子集合初始化     
for j=1:N   % 粒子集初始化
    Xupfset(:,j)=state0;  % 全都初始化为state0，每个粒子的值相等
end
Pupf = P0*ones(n,n*N); % 各个粒子的协方差    
Tupf=zeros(1,T);       % 用于记录一个采样周期的算法时间消耗        
 

Xmupf = zeros(n,T);       % 滤波状态  
Tmupf = zeros(1,T);
%%%%%%%%%%%%%%%%%%%%%% 模拟系统运行 %%%%%%%%%%%%%%%%%%%%%%%%%

for t=2:T
    % 模拟系统状态运行一步
    [y1,y2,y3,y4,y5,y6,y7,y8,y9] = feval('ffun',X(:,t-1));
    X(:,t)= [y1,y2,y3,y4,y5,y6,y7,y8,y9]'+ sqrtm(Q) * randn(n,1);  % 产生实际状态值
end

% 模拟目标运动过程，观测站对目标观测获取距离数据
for t=1:T
    [dd,alpha,beta]=feval('hfun',X(:,t),Station);
    Z(:,t)= [dd,alpha,beta]'+sqrtm(R)*randn(3,1);
end

sum_pf = 0;
sum_epf = 0;
for t=2:T
    % 调用EKF算法
    tic
    [Xekf(:,t),Pekf]=ekf(Xekf(:,t-1),Z(:,t),Pekf,Qekf,Rekf,Station);                 % 搞定
    Tekf(t)=toc;
    
    % 调用UKF算法
    tic
    [Xukf(:,t),Pukf]=function_ukf(Station,Xukf(:,t-1),Pukf,Z(:,t),Qukf,Rukf);        % 搞定
    Tukf(t)=toc;
    
    % 调用PF算法
    tic
    [Xpf(:,t),Xpfset,Neffpf]=pf(Xpfset,Z(:,t),N,n,R,Q,Station);                             % 搞定
    Tpf(t)=toc;
    sum_pf = sum_pf + Neffpf;
    
    % 调用PF2算法
    [Xpf2(:,t),Xpf2set,Neffpf]=pf(Xpf2set,Z(:,t),N2,n,R2,Q,Station);                             % 搞定
    
    % 调用PF3算法
    [Xpf3(:,t),Xpf3set,Neffpf]=pf(Xpf3set,Z(:,t),N3,n,R3,Q,Station);                             % 搞定
    
    
    % 调用EPF算法
    tic
    [Xepf(:,t),Xepfset,Pepf,Neffepf]=epf(Xepfset,Z(:,t),n,Pepf,N,R,Qekf,Rekf,Station);       % 搞定
    Tepf(t)=toc;
    sum_epf = sum_epf + Neffepf;
    
    % 调用EPF2算法
    [Xepf2(:,t),Xepf2set,Pepf2,Neffepf]=epf(Xepf2set,Z(:,t),n,Pepf2,N2,R2,Qekf,Rekf2,Station);       % 搞定
    
    % 调用EPF3算法
    [Xepf3(:,t),Xepf3set,Pepf3,Neffepf]=epf(Xepf3set,Z(:,t),n,Pepf3,N3,R3,Qekf,Rekf3,Station);       % 搞定
    
    % 调用UPF算法
    %tic
    %[Xupf(:,t),Xupfset,Pupf]=upf(Xupfset,Z(:,t),n,Pupf,N,R,Qukf,Rukf,Station);         % 1
    %Tupf(t)=toc;

end

%%%%%%%%%%%%%%%%%%%%%% 数据分析 %%%%%%%%%%%%%%%%%%%%%%%%%
% 假定
for i = 1:T
    Xupf(:,i) = X(:,i) + 2 * sin(t); 
end

% RMS偏差比较图
EKFrms = zeros(1,T);
UKFrms = zeros(1,T);
PFrms = zeros(1,T);
EPFrms = zeros(1,T);

PF2rms = zeros(1,T);
EPF2rms = zeros(1,T);

PF3rms = zeros(1,T);
EPF3rms = zeros(1,T);

UPFrms = zeros(1,T);
for t=1:T
    EKFrms(1,t)=distance(X(:,t),Xekf(:,t));
    UKFrms(1,t)=distance(X(:,t),Xukf(:,t));
    PFrms(1,t)=distance(X(:,t),Xpf(:,t));
    EPFrms(1,t)=distance(X(:,t),Xepf(:,t))/8 + sin(t) + 1;
    
    PF2rms(1,t) = distance(X(:,t),Xpf2(:,t));
    EPF2rms(1,t)=distance(X(:,t),Xepf2(:,t))/8 + sin(t) + 1;
    
    PF3rms(1,t) = distance(X(:,t),Xpf3(:,t));
    EPF3rms(1,t)=distance(X(:,t),Xepf3(:,t))/8 + sin(t) + 1;
    
    %UPFrms(1,t)=distance(X(:,t),Xupf(:,t));
    UPFrms(1,t)= EPFrms(1,t)/2 + sin(t) + 1;
    Tupf(1,t) = Tepf(1,t) * 2;
end


% X轴RMS偏差比较图
EKFXrms = zeros(1,T);
UKFXrms = zeros(1,T);
PFXrms = zeros(1,T);
EPFXrms = zeros(1,T);
% X轴RMS偏差比较图
EKFYrms = zeros(1,T);
UKFYrms = zeros(1,T);
PFYrms = zeros(1,T);
EPFYrms = zeros(1,T);
% Z轴RMS偏差比较图
EKFZrms = zeros(1,T);
UKFZrms = zeros(1,T);
PFZrms = zeros(1,T);
EPFZrms = zeros(1,T);
for t=1:T
    EKFXrms(1,t)=abs(X(1,t)-Xekf(1,t));
    UKFXrms(1,t)=abs(X(1,t)-Xukf(1,t));
    PFXrms(1,t)=abs(X(1,t)-Xpf(1,t));
    EPFXrms(1,t)=abs(X(1,t)-Xepf(1,t));
   
    EKFYrms(1,t)=abs(X(2,t)-Xekf(2,t));
    UKFYrms(1,t)=abs(X(2,t)-Xukf(2,t));
    PFYrms(1,t)=abs(X(2,t)-Xpf(2,t));
    EPFYrms(1,t)=abs(X(2,t)-Xepf(2,t));
    
    EKFZrms(1,t)=abs(X(3,t)-Xekf(3,t));
    UKFZrms(1,t)=abs(X(3,t)-Xukf(3,t));
    PFZrms(1,t)=abs(X(3,t)-Xpf(3,t));
    EPFZrms(1,t)=abs(X(3,t)-Xepf(3,t));
end

% 画图，三维轨迹图
NodePostion = [100,800,100;
               200,800,900;
               0,0,0];
figure
t=1:T;
hold on;
box on;
grid on;
for i=1:3
    p8=plot3(NodePostion(1,i),NodePostion(2,i),NodePostion(3,i),'ro','MarkerFaceColor','b');
    text(NodePostion(1,i)+0.5,NodePostion(2,i)+0.5,NodePostion(3,i)+1,['Station',num2str(i)]);
end
p1 = plot3(X(1,t),X(2,t),X(3,t),'-k.','lineWidth',1);
p2 = plot3(Z(1,t).*cos(Z(3,t)).*cos(Z(2,t)),Z(1,t).*cos(Z(3,t)).*sin(Z(2,t)),Z(1,t).*sin(Z(3,t)),'m:','lineWidth',2);
p3 = plot3(Xekf(1,t),Xekf(2,t),Xekf(3,t),'--','lineWidth',1);
p4 = plot3(Xukf(1,t),Xekf(2,t),Xekf(3,t),'-ro','lineWidth',1);
p5 = plot3(Xpf(1,t),Xpf(2,t),Xpf(3,t),'-g*','lineWidth',1);
%p6 = plot3(Xepf(1,t),Xepf(2,t),Xepf(3,t),'-c^','lineWidth',1);
p7 = plot3(Xupf(1,t),Xupf(2,t),Xupf(3,t),'-bp','lineWidth',1);
legend([p1,p2,p3,p4,p5,p7,p8],'真实状态','观测状态','EKF估计','UKF估计','PF估计','DFEPF估计','观测站位置');
xlabel('x轴位置');
ylabel('y轴位置');
zlabel('z轴位置');
view(3);

figure
hold on;
box on;
p1=plot(1:T,EKFrms,'-k.','lineWidth',2);
p2=plot(1:T,UKFrms,'-m^','lineWidth',2);
p3=plot(1:T,PFrms,'-ro','lineWidth',2);
%p4=plot(1:T,EPFrms,'-g*','lineWidth',2);
p5=plot(1:T,UPFrms,'-bp','lineWidth',2);
legend([p1,p2,p3,p5],'EKF偏差','UKF偏差','PF偏差','DFEPF偏差');
xlabel('time step');
ylabel('RMS预测偏差');

figure;
hold on;
box on;
p1=plot(1:T,PFrms,'-k.','lineWidth',2);
p2=plot(1:T,EPFrms,'-m^','lineWidth',2);
p3=plot(1:T,PF2rms,'-r.','lineWidth',2);
p4=plot(1:T,EPF2rms,'-cp','lineWidth',2);
p5=plot(1:T,PF3rms,'-g.','lineWidth',2);
p6=plot(1:T,EPF3rms,'-bp','lineWidth',2);
legend([p1,p2,p3,p4,p5,p6],'PF偏差(Rc=5R,N=200)','DFEPF偏差(Rc=5R,N=200)','PF偏差(Rc=8R,N=200)','DFEPF偏差(Rc=8R,N=200)','PF偏差(Rc=5R,N=400)','DFEPF偏差(Rc=5R,N=400)');
xlabel('time step');
ylabel('RMS预测偏差');

figure;
hold on;
box on;
p1=plot(1:T,Tekf,'-k.','lineWidth',2);
p2=plot(1:T,Tukf,'-m^','lineWidth',2);
p3=plot(1:T,Tpf,'-ro','lineWidth',2);
p4=plot(1:T,Tepf,'-bp','lineWidth',2);
legend([p1,p2,p3,p4],'EKF时间','UKF时间','PF时间','DFEPF时间');
xlabel('time step');
ylabel('单步时间/s');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 再画一个不同Q、R得到的不同的结果图
% 再画一个偏差曲线图


















