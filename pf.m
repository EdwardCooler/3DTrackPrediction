%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ���������˲��㷨
% ���룺XisetΪ��ά����
% �����XoΪnx1����XosetΪnxN����
function [Xo,Xoset,Neff]=pf(Xiset,Z,N,n,R,Q,Station)
 
tic
% �ز�������
resamplingScheme=1;

% �м������ʼ��
Zpre=ones(3,N);     % �۲�Ԥ��   
Xsetpre=ones(n,N);  % ���Ӽ���Ԥ��
w = ones(1,N);      % Ȩֵ��ʼ��
Xo=zeros(n,1);

% ��һ��������ÿһ�����Ӷ�����ֲ����� 
for i=1:N
    [y1,y2,y3,y4,y5,y6,y7,y8,y9] = ffun(Xiset(:,i));
    Xsetpre(:,i) = [y1,y2,y3,y4,y5,y6,y7,y8,y9]' + 10*sqrtm(Q)*randn(n,1);
end

% �ڶ�������������Ȩ��
for i=1:N
    [dd,alpha,beta]=feval('hfun',Xsetpre(:,i),Station);
    Zpre(:,i) =[dd,alpha,beta]';
    z1 = Z-Zpre(:,i);
    % w(i) = inv(sqrtm(R)) * exp(-0.5*inv(R)*((Z-Zpre(:,i))^(2))) + 1e-99; 
    w(i) = inv(sqrt(2*pi*det(R)))*exp(-.5*(z1)'*inv(R)*(z1))+ 1e-99;%Ȩֵ���㣬������ʵ���԰�inv(sqrt(2*pi*det(R)))ȥ��
end

% ��������Ȩ�ع�һ��
w = w./sum(w);   

% ����Ч������
Neff = 1/sum(w.^2);

% ���������֮���ͼ��������

% ���Ĳ����ز���
if resamplingScheme == 1
    outIndex = residualR(1:N,w');       
elseif resamplingScheme == 2
    outIndex = systematicR(1:N,w');  
else
    outIndex = multinomialR(1:N,w');  
end

% ���岽���������Ӽ��� 
Xoset = Xsetpre(:,outIndex); 
% ���������õ����μ���������˲�����ֵ
target=[mean(Xoset(1,:)),mean(Xoset(2,:)),mean(Xoset(3,:)),mean(Xoset(4,:)),mean(Xoset(5,:)),mean(Xoset(6,:)),mean(Xoset(7,:)),mean(Xoset(8,:)),mean(Xoset(9,:))]';
Xo(:,1)=target;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


