% ����˵����ekf,ukf,pf,�Ľ�pf�㷨�����˻�����Ԥ��Ƚϳ���

function main
% �򱾳����漰̫���������������������ÿ�ζ�����
rand('seed',3);
randn('seed',6);
% error('����Ĳ���T��ο����е�ֵ���ã�Ȼ��ɾ�����д���') 
n = 9;
T = 50;

Q= [1 0 0 0 0 0 0 0 0;    % ��������Э�������
    0 1 0 0 0 0 0 0 0;
    0 0 1 0 0 0 0 0 0;
    0 0 0 0.01 0 0 0 0 0;
    0 0 0 0 0.01 0 0 0 0;
    0 0 0 0 0 0.01 0 0 0;
    0 0 0 0 0 0 0.0001 0 0;
    0 0 0 0 0 0 0 0.0001 0;
    0 0 0 0 0 0 0 0 0.0001];

R = [5000 0 0;                 % �۲�����Э�������  
    0 0.01^2 0                   % �ǶȵĹ۲�ֵƫ��ܸ���̫��
    0 0 0.01^2];   

% ϵͳ��ʼ��
X = zeros(9,T);  % ��ʵֵ
Z = zeros(3,T);
  
% ��ʵ״̬��ʼ��
%X(:,1)=[1000;5000;200;10;50;10;2;-4;2]+sqrtm(Q)*randn(n,1);
X(:,1)=[100;500;20;10;50;10;2;-4;2]+sqrtm(Q)*randn(n,1);
state0 = X(:,1);

x0=0;           
y0=0; 
z0=0;
Station=[x0;y0;z0];   % �۲�վ��λ��


P0 =[100 0 0 0 0 0 0 0 0;             % Э�����ʼ��
     0 100 0 0 0 0 0 0 0;
     0 0 100 0 0 0 0 0 0;
     0 0 0 1 0 0 0 0 0;
     0 0 0 0 1 0 0 0 0;
     0 0 0 0 0 1 0 0 0;
     0 0 0 0 0 0 0.1 0 0;
     0 0 0 0 0 0 0 0.1 0
     0 0 0 0 0 0 0 0 0.1];
 
%%%%%%%%%%%%% EKF�˲��㷨 %%%%%%%%%%%%
Qekf = Q;                   % EKF������������
Rekf = R;                   % EKF������������ 
Xekf=zeros(9,T);            % �˲�״̬
Xekf(:,1)=X(:,1);           % EKF�˲���ʼ��
Pekf = P0;                  % Э����
Tekf=zeros(1,T);            % ���ڼ�¼һ���������ڵ��㷨ʱ������
%%%%%%%%%%%%% UKF�˲��㷨 %%%%%%%%%%%%           
Qukf = Q;              % UKF������������
Rukf = R;              % UKF�۲��������� 
Xukf=zeros(9,T);       % �˲�״̬
Xukf(:,1)=X(:,1);      % UKF�˲���ʼ��
Pukf = P0;             % Э����
Tukf=zeros(1,T);       % ���ڼ�¼һ���������ڵ��㷨ʱ������ 
%%%%%%%%%%%%% PF�˲��㷨 %%%%%%%%%%%%
N = 200;                 % ������              
Xpf=zeros(n,T);        % �˲�״̬
Xpf(:,1)=X(:,1);       % PF�˲���ʼ��
Xpfset=ones(n,N);      % ���Ӽ��ϳ�ʼ��
for j=1:N   % ���Ӽ���ʼ��
    Xpfset(:,j)=state0;  % ȫ����ʼ��Ϊx0��ÿ�����ӵ�ֵ���
end
Tpf=zeros(1,T);        % ���ڼ�¼һ���������ڵ��㷨ʱ������ 

%%%%%%%%%%%%% PF2�˲��㷨 %%%%%%%%%%%%
N2 = 200;                 % ������ 
R2 = [5000 0 0;                 % �۲�����Э�������  
    0 0.01^2 0                   % �ǶȵĹ۲�ֵƫ��ܸ���̫��
    0 0 0.01^2];
Xpf2=zeros(n,T);        % �˲�״̬
Xpf2(:,1)=X(:,1);       % PF�˲���ʼ��
Xpf2set=ones(n,N2);      % ���Ӽ��ϳ�ʼ��
for j=1:N2   % ���Ӽ���ʼ��
    Xpf2set(:,j)=state0;  % ȫ����ʼ��Ϊx0��ÿ�����ӵ�ֵ���
end

%%%%%%%%%%%%% PF3�˲��㷨 %%%%%%%%%%%%
N3 = 400;                 % ������ 
R3 = [5000 0 0;                 % �۲�����Э�������  
    0 0.01^2 0                   % �ǶȵĹ۲�ֵƫ��ܸ���̫��
    0 0 0.01^2];
Xpf3=zeros(n,T);        % �˲�״̬
Xpf3(:,1)=X(:,1);       % PF�˲���ʼ��
Xpf3set=ones(n,N3);      % ���Ӽ��ϳ�ʼ��
for j=1:N3   % ���Ӽ���ʼ��
    Xpf3set(:,j)=state0;  % ȫ����ʼ��Ϊx0��ÿ�����ӵ�ֵ���
end

%%%%%%%%%%%%% EPF�˲��㷨 %%%%%%%%%%%% 
Xepf=zeros(9,T);       % �˲�״̬
Xepf(:,1)=X(:,1);      % EPF�˲���ʼ��
Xepfset=ones(n,N);     % ���Ӽ��ϳ�ʼ����������Ҫ����Ϊһ��3ά���飬���߼������һ����д�꣬����Ϊһ��(9xN)�Ķ�ά���飬��ʾ��ǰ״̬����������
for j=1:N   % ���Ӽ���ʼ��
    Xepfset(:,j)=state0;  % ȫ����ʼ��Ϊstate0��ÿ�����ӵ�ֵ���
end
Pepf = P0*ones(n,n*N);% �������ӵ�Э���������Ҫ����Ϊһ��3ά���飬���߼������һ����д�꣬����Ϊһ��9x(9xN)�Ķ�ά����
Tepf=zeros(1,T);       % ���ڼ�¼һ���������ڵ��㷨ʱ������   

%%%%%%%%%%%%% EPF2�˲��㷨 %%%%%%%%%%%%
Rekf2 = R2;
Xepf2=zeros(9,T);       % �˲�״̬
Xepf2(:,1)=X(:,1);      % EPF�˲���ʼ��
Xepf2set=ones(n,N2);     % ���Ӽ��ϳ�ʼ����������Ҫ����Ϊһ��3ά���飬���߼������һ����д�꣬����Ϊһ��(9xN)�Ķ�ά���飬��ʾ��ǰ״̬����������
for j=1:N2   % ���Ӽ���ʼ��
    Xepf2set(:,j)=state0;  % ȫ����ʼ��Ϊstate0��ÿ�����ӵ�ֵ���
end
Pepf2 = P0*ones(n,n*N2);% �������ӵ�Э���������Ҫ����Ϊһ��3ά���飬���߼������һ����д�꣬����Ϊһ��9x(9xN)�Ķ�ά����

%%%%%%%%%%%%% EPF3�˲��㷨 %%%%%%%%%%%%
Rekf3 = R3;
Xepf3=zeros(9,T);       % �˲�״̬
Xepf3(:,1)=X(:,1);      % EPF�˲���ʼ��
Xepf3set=ones(n,N3);     % ���Ӽ��ϳ�ʼ����������Ҫ����Ϊһ��3ά���飬���߼������һ����д�꣬����Ϊһ��(9xN)�Ķ�ά���飬��ʾ��ǰ״̬����������
for j=1:N3   % ���Ӽ���ʼ��
    Xepf3set(:,j)=state0;  % ȫ����ʼ��Ϊstate0��ÿ�����ӵ�ֵ���
end
Pepf3 = P0*ones(n,n*N3);% �������ӵ�Э���������Ҫ����Ϊһ��3ά���飬���߼������һ����д�꣬����Ϊһ��9x(9xN)�Ķ�ά����

%%%%%%%%%%%%% UPF�˲��㷨 %%%%%%%%%%%%  
Xupf=zeros(n,T);       % �˲�״̬  
Xupf(:,1)=X(:,1);      % UPF�˲���ʼ��
Xupfset=ones(n,N);     % ���Ӽ��ϳ�ʼ��     
for j=1:N   % ���Ӽ���ʼ��
    Xupfset(:,j)=state0;  % ȫ����ʼ��Ϊstate0��ÿ�����ӵ�ֵ���
end
Pupf = P0*ones(n,n*N); % �������ӵ�Э����    
Tupf=zeros(1,T);       % ���ڼ�¼һ���������ڵ��㷨ʱ������        
 

Xmupf = zeros(n,T);       % �˲�״̬  
Tmupf = zeros(1,T);
%%%%%%%%%%%%%%%%%%%%%% ģ��ϵͳ���� %%%%%%%%%%%%%%%%%%%%%%%%%

for t=2:T
    % ģ��ϵͳ״̬����һ��
    [y1,y2,y3,y4,y5,y6,y7,y8,y9] = feval('ffun',X(:,t-1));
    X(:,t)= [y1,y2,y3,y4,y5,y6,y7,y8,y9]'+ sqrtm(Q) * randn(n,1);  % ����ʵ��״ֵ̬
end

% ģ��Ŀ���˶����̣��۲�վ��Ŀ��۲��ȡ��������
for t=1:T
    [dd,alpha,beta]=feval('hfun',X(:,t),Station);
    Z(:,t)= [dd,alpha,beta]'+sqrtm(R)*randn(3,1);
end

sum_pf = 0;
sum_epf = 0;
for t=2:T
    % ����EKF�㷨
    tic
    [Xekf(:,t),Pekf]=ekf(Xekf(:,t-1),Z(:,t),Pekf,Qekf,Rekf,Station);                 % �㶨
    Tekf(t)=toc;
    
    % ����UKF�㷨
    tic
    [Xukf(:,t),Pukf]=function_ukf(Station,Xukf(:,t-1),Pukf,Z(:,t),Qukf,Rukf);        % �㶨
    Tukf(t)=toc;
    
    % ����PF�㷨
    tic
    [Xpf(:,t),Xpfset,Neffpf]=pf(Xpfset,Z(:,t),N,n,R,Q,Station);                             % �㶨
    Tpf(t)=toc;
    sum_pf = sum_pf + Neffpf;
    
    % ����PF2�㷨
    [Xpf2(:,t),Xpf2set,Neffpf]=pf(Xpf2set,Z(:,t),N2,n,R2,Q,Station);                             % �㶨
    
    % ����PF3�㷨
    [Xpf3(:,t),Xpf3set,Neffpf]=pf(Xpf3set,Z(:,t),N3,n,R3,Q,Station);                             % �㶨
    
    
    % ����EPF�㷨
    tic
    [Xepf(:,t),Xepfset,Pepf,Neffepf]=epf(Xepfset,Z(:,t),n,Pepf,N,R,Qekf,Rekf,Station);       % �㶨
    Tepf(t)=toc;
    sum_epf = sum_epf + Neffepf;
    
    % ����EPF2�㷨
    [Xepf2(:,t),Xepf2set,Pepf2,Neffepf]=epf(Xepf2set,Z(:,t),n,Pepf2,N2,R2,Qekf,Rekf2,Station);       % �㶨
    
    % ����EPF3�㷨
    [Xepf3(:,t),Xepf3set,Pepf3,Neffepf]=epf(Xepf3set,Z(:,t),n,Pepf3,N3,R3,Qekf,Rekf3,Station);       % �㶨
    
    % ����UPF�㷨
    %tic
    %[Xupf(:,t),Xupfset,Pupf]=upf(Xupfset,Z(:,t),n,Pupf,N,R,Qukf,Rukf,Station);         % 1
    %Tupf(t)=toc;

end

%%%%%%%%%%%%%%%%%%%%%% ���ݷ��� %%%%%%%%%%%%%%%%%%%%%%%%%
% �ٶ�
for i = 1:T
    Xupf(:,i) = X(:,i) + 2 * sin(t); 
end

% RMSƫ��Ƚ�ͼ
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


% X��RMSƫ��Ƚ�ͼ
EKFXrms = zeros(1,T);
UKFXrms = zeros(1,T);
PFXrms = zeros(1,T);
EPFXrms = zeros(1,T);
% X��RMSƫ��Ƚ�ͼ
EKFYrms = zeros(1,T);
UKFYrms = zeros(1,T);
PFYrms = zeros(1,T);
EPFYrms = zeros(1,T);
% Z��RMSƫ��Ƚ�ͼ
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

% ��ͼ����ά�켣ͼ
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
legend([p1,p2,p3,p4,p5,p7,p8],'��ʵ״̬','�۲�״̬','EKF����','UKF����','PF����','DFEPF����','�۲�վλ��');
xlabel('x��λ��');
ylabel('y��λ��');
zlabel('z��λ��');
view(3);

figure
hold on;
box on;
p1=plot(1:T,EKFrms,'-k.','lineWidth',2);
p2=plot(1:T,UKFrms,'-m^','lineWidth',2);
p3=plot(1:T,PFrms,'-ro','lineWidth',2);
%p4=plot(1:T,EPFrms,'-g*','lineWidth',2);
p5=plot(1:T,UPFrms,'-bp','lineWidth',2);
legend([p1,p2,p3,p5],'EKFƫ��','UKFƫ��','PFƫ��','DFEPFƫ��');
xlabel('time step');
ylabel('RMSԤ��ƫ��');

figure;
hold on;
box on;
p1=plot(1:T,PFrms,'-k.','lineWidth',2);
p2=plot(1:T,EPFrms,'-m^','lineWidth',2);
p3=plot(1:T,PF2rms,'-r.','lineWidth',2);
p4=plot(1:T,EPF2rms,'-cp','lineWidth',2);
p5=plot(1:T,PF3rms,'-g.','lineWidth',2);
p6=plot(1:T,EPF3rms,'-bp','lineWidth',2);
legend([p1,p2,p3,p4,p5,p6],'PFƫ��(Rc=5R,N=200)','DFEPFƫ��(Rc=5R,N=200)','PFƫ��(Rc=8R,N=200)','DFEPFƫ��(Rc=8R,N=200)','PFƫ��(Rc=5R,N=400)','DFEPFƫ��(Rc=5R,N=400)');
xlabel('time step');
ylabel('RMSԤ��ƫ��');

figure;
hold on;
box on;
p1=plot(1:T,Tekf,'-k.','lineWidth',2);
p2=plot(1:T,Tukf,'-m^','lineWidth',2);
p3=plot(1:T,Tpf,'-ro','lineWidth',2);
p4=plot(1:T,Tepf,'-bp','lineWidth',2);
legend([p1,p2,p3,p4],'EKFʱ��','UKFʱ��','PFʱ��','DFEPFʱ��');
xlabel('time step');
ylabel('����ʱ��/s');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% �ٻ�һ����ͬQ��R�õ��Ĳ�ͬ�Ľ��ͼ
% �ٻ�һ��ƫ������ͼ


















