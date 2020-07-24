%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UT�任�Ӻ���
% ���룺funΪ���������XsigmaΪ��������Wm��WcΪȨֵ��nΪ״̬ά����COVΪЭ����
% �����XmeansΪ��ֵ��PΪЭ���XdivΪԤ���ȥ��ֵ
function [Xmeans,Xsigma_pre,P,Xdiv]=function_ut(Station,flag,Xsigma,Wm,Wc,n,COV)
LL=size(Xsigma,2);      % �õ�Xsigma����������size(Xsigma,2)����˼�Ƿ��ؾ���Xsigma������
Xmeans=zeros(n,1);      % ��ֵ
Xsigma_pre=zeros(n,LL); % һ��Ԥ��
for k=1:LL 
    % ������Ҫ�ж�һ��Ӧ�õ�����һ������
    if flag==1  % ����ffun
        [y1,y2,y3,y4,y5,y6,y7,y8,y9] = feval('ffun',Xsigma(:,k));
        Xsigma_pre(:,k) = [y1,y2,y3,y4,y5,y6,y7,y8,y9]';  % һ��Ԥ�⣬�����fun������ǰ���f����h����
    else % ����hfun
        [dd,alpha,beta] = feval('hfun',Xsigma(:,k),Station);
        Xsigma_pre(:,k) = [dd,alpha,beta]';  % һ��Ԥ�⣬�����fun������ǰ���f����h����
    end
    Xmeans=Xmeans+Wm(k)*Xsigma_pre(:,k);
end
% Xmeans(:,ones(1,LL))��Xmeans��չ��n*LL����ÿһ�ж����
Xdiv=Xsigma_pre-Xmeans(:,ones(1,LL));  % Ԥ���ȥ��ֵ
P=Xdiv*diag(Wc)*Xdiv'+COV;             % Э����