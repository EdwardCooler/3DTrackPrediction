%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ��EKF��������ֲ�
%  �������˵��
%  Xiset����t-1ʱ�̵����Ӽ���
%  Z��tʱ�̵Ĺ۲�
%  Pin��Xiset���Ӽ��ϵķ����

% �������˵��
% Xo��EPF�㷨���յĹ��ƽ��
% Xoset��kʱ�̵����ӽ�ϣ����ֵ����Xo
% Pout��Xoset��Ӧ�ķ���
function [Xo,Xoset,Pout,Neff]=epf(Xiset,Z,n,Pin,N,R,Qekf,Rekf,Station)
% �ز������Ժ���
resamplingScheme=1;
% �м������ʼ��
Zpre=ones(3,N);    % �۲�Ԥ��     
Xsetpre=ones(n,N); % ���Ӽ���Ԥ��
w = ones(1,N);     % Ȩֵ��ʼ��
Xo=zeros(n,1);

Pout=ones(n,n*N);    % Э����Ԥ�⣬���PoutҲ��һ��9x(9xN)�ľ�����Ҫ�޸�
% Xekf=ones(n,1);    % EKF���ƽ��
% Xekf_pre=ones(1,N);% EKF��һ��Ԥ��

% ��һ��������EKF����õ��Ľ�����в���
for i=1:N
    % ����EKF����õ�ÿ�����ӵľ�ֵ��Э����ע������Ž��������Ӽ���Xiset��(9xN)�ģ�Э����Pin��9x(9xN)�ģ�������ֿ�ΪN��9x9�ľ�����ô�ֿ��أ�
    % Pin(:,9.*(i-1)+1:9.*(i-1)+9)
    [Xekf,Pout(:,9*(i-1)+1:9*(i-1)+9)]=ekf(Xiset(:,i),Z,Pin(:,9*(i-1)+1:9*(i-1)+9),Qekf,Rekf,Station);
    %[Xekf,Pout]=ekf(Xin,Z,Pin,Qekf,Rekf,Station)
    % ������õ��ľ�ֵ�ͷ�����Ϊ���Ӽ��ϲ���
    Xsetpre(:,i)=Xekf + sqrtm(Pout(:,9*(i-1)+1:9*(i-1)+9)) * randn(n,1);
end

% �ڶ���������Ȩ��
for i=1:N
    % �۲�Ԥ��
    [dd,alpha,beta]=feval('hfun',Xsetpre(:,i),Station);
    Zpre(:,i) =[dd,alpha,beta]';
    z1 = Z-Zpre(:,i);
    % ����Ȩ�أ�1e-99Ϊ��С��0���֣���ֹ��0
    lik = inv(sqrt(2*pi*det(R)))*exp(-.5*(z1)'*inv(R)*(z1))+ 1e-99;
    % ����
    %prior = ((Xsetpre(i)-Xiset(i))^(g1-1)) * exp(-g2*(Xsetpre(i)-Xiset(i)));
    % ����ֲ�
    %proposal = inv(sqrt(Pout(i))) * exp(-0.5*inv(Pout(i)) *((Xsetpre(i)-Xekf(i))^(2)));
    %w(i) = lik*prior/proposal;
    w(i) = lik;
end
% Ȩֵ��һ��
w= w./sum(w);

summ = 0;
alfa = 0.2;
for i=1:N
    summ = summ + w(i)^alfa; 
end
for i=1:N
    w(i) =  w(i)^alfa / summ; 
end

% ����Ч������
Neff = 1/sum(w.^2);

% ���������ز���
if resamplingScheme == 1
    outIndex = residualR(1:N,w');   
elseif resamplingScheme == 2
    outIndex = systematicR(1:N,w');     
else
    outIndex = multinomialR(1:N,w');    
end

% ���Ĳ����������Ӽ��� 
Xoset = Xsetpre(:,outIndex); 
% �������ӷ���
Pm = ones(n,n*N);      % Э����Ԥ��
for i = 1:N
    Pm(:,9*(i-1)+1:9*(i-1)+9) = Pout(:,9*(outIndex(i)-1)+1:9*(outIndex(i)-1)+9);
end
Pout = Pm;
% ���岽���õ����μ���������˲�����ֵ
target=[mean(Xoset(1,:)),mean(Xoset(2,:)),mean(Xoset(3,:)),mean(Xoset(4,:)),mean(Xoset(5,:)),mean(Xoset(6,:)),mean(Xoset(7,:)),mean(Xoset(8,:)),mean(Xoset(9,:))]';
Xo(:,1)=target;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


