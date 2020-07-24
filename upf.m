%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ��UKF�Ľ��������˲��㷨--EPF

% ��UKF��������ֲ�
% �������˵����
%    Xiset����t-1ʱ�̵����Ӽ��ϣ�Z��tʱ�̵Ĺ۲�
%    Pin��ӦXiset���Ӽ��ϵķ���
% �������˵����
%    Xo��upf�㷨���յĹ��ƽ��
%    Xoset��kʱ�̵����Ӽ��ϣ����ֵ����Xo
%    Pout��Xoset��Ӧ�ķ���
function [Xo,Xoset,Pout]=upf(Xiset,Z,n,Pin,N,R,Qukf,Rukf,Station)
 
% �ز������Ժ���
resamplingScheme=1;

% �м������ʼ��
Zpre=ones(3,N);      % �۲�Ԥ�� 
Xset_pre=ones(n,N);  % ���Ӽ���Ԥ��
w = ones(1,N);     % Ȩֵ��ʼ��
Xo=zeros(n,1);

Pout=ones(n,n*N);    % Э����Ԥ�⣬���PoutҲ��һ��9x(9xN)�ľ�����Ҫ�޸�
% Xukf=ones(1,N);      % UKF���ƽ�� 
% Xukf_pre=ones(1,N);  % EKF��һ��Ԥ��

% ��һ��������UKF����õ��Ľ�����в���
for i=1:N
     % ����UKF����õ���ֵ�ͷ����i�����ӵľ�ֵ�ͷ���
    [Xukf,Pout(:,9*(i-1)+1:9*(i-1)+9)]=function_ukf(Station,Xiset(:,i),Pin(:,9*(i-1)+1:9*(i-1)+9),Z,Qukf,Rukf);
    % [X,P]=function_ukf(Station,X,P,Z,Q,R)
    % ������õ��ľ�ֵ�ͷ�����Ϊ���Ӽ��ϲ���
    Xset_pre(i) = Xukf + sqrtm(Pout(:,9*(i-1)+1:9*(i-1)+9)) * randn(n,1);
end

% �ڶ���������Ȩ��
for i=1:N
    % �۲�Ԥ��
    [dd,alpha,beta]=feval('hfun',Xset_pre(:,i),Station);
    Zpre(:,i) =[dd,alpha,beta]';
    z1 = Z-Zpre(:,i);
    % ����Ȩ�أ�1e-99Ϊ��С��0���֣���ֹ��0
    lik = inv(sqrt(2*pi*det(R)))*exp(-.5*(z1)'*inv(R)*(z1))+ 1e-99;
    % prior = ((Xset_pre(i)-Xiset(i))^(g1-1)) * exp(-g2*(Xset_pre(i)-Xiset(i)));
    % proposal = inv(sqrt(Pout(i))) * exp(-0.5*inv(Pout(i)) *((Xset_pre(i)-Xukf(i))^(2)));
    % w(i) = lik*prior/proposal;
    w(i) = lik;
end

% Ȩֵ��һ�� 
w = w./sum(w);

% ���������ز��� 
if resamplingScheme == 1
    outIndex = residualR(1:N,w');        
elseif resamplingScheme == 2
    outIndex = systematicR(1:N,w');      
else
    outIndex = multinomialR(1:N,w');     
end

% ���Ĳ����������Ӽ��� 
Xoset = Xset_pre(:,outIndex); 
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


