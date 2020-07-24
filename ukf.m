%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% �޼��������˲��㷨
function [Xout,Pout]=ukf(Xin,Z,Pin,Qukf,Rukf,t) 
% �޼��任�Ĳ���
alpha = 1;
beta  = 0;
kappa = 2;

states       = size(Xin(:),1);
observations = size(Z(:),1);
vNoise       = size(Qukf,2);
wNoise       = size(Rukf,2);
noises       = vNoise+wNoise;

% ������������״̬����������չ 
% ע�⣬�˴���Ϊ��������ģ�ͣ�Ϊ���ǽ�Լ���㿪֧
if (noises)
    N=[Qukf zeros(vNoise,wNoise); zeros(wNoise,vNoise) Rukf];
    PQ=[Pin zeros(states,noises);zeros(noises,states) N];
    xQ=[Xin;zeros(noises,1)];
else
    PQ=Pin;
    xQ=Xin;
end

% ͨ��UT�任������sigma�㼯�����ǵ�Ȩֵ
[xSigmaPts, wSigmaPts, nsp] = scaledSymmetricSigmaPoints(xQ, PQ, alpha, beta, kappa);
 
% Ϊ�˼������У�����ִ��Ч�ʣ�����wSigmaPts���Ƶ�������
wSigmaPts_xmat = repmat(wSigmaPts(:,2:nsp),states,1);
wSigmaPts_zmat = repmat(wSigmaPts(:,2:nsp),observations,1);

% ����sigma�㼯Ԥ��ֵ
xPredSigmaPts = feval('ffun',xSigmaPts(1:states,:),t)+xSigmaPts(states+1:states+vNoise,:);
zPredSigmaPts = feval('hfun',xPredSigmaPts,t)+xSigmaPts(states+vNoise+1:states+noises,:);

% �����ֵ��Ԥ��ֵ
xPred = sum(wSigmaPts_xmat .* (xPredSigmaPts(:,2:nsp) - repmat(xPredSigmaPts(:,1),1,nsp-1)),2);
zPred = sum(wSigmaPts_zmat .* (zPredSigmaPts(:,2:nsp) - repmat(zPredSigmaPts(:,1),1,nsp-1)),2);
xPred=xPred+xPredSigmaPts(:,1);
zPred=zPred+zPredSigmaPts(:,1);
 
% ���㷽���Э���ע���һ�е�Ȩֵ���ֵ��Ȩ�ǲ�һ����
% �˴���Ҫ�����޼��仯�Ĺ���
exSigmaPt = xPredSigmaPts(:,1)-xPred;
ezSigmaPt = zPredSigmaPts(:,1)-zPred;

PPred   = wSigmaPts(nsp+1)*exSigmaPt*exSigmaPt';
PxzPred = wSigmaPts(nsp+1)*exSigmaPt*ezSigmaPt';
S       = wSigmaPts(nsp+1)*ezSigmaPt*ezSigmaPt';

exSigmaPt = xPredSigmaPts(:,2:nsp) - repmat(xPred,1,nsp-1);
ezSigmaPt = zPredSigmaPts(:,2:nsp) - repmat(zPred,1,nsp-1);
PPred     = PPred + (wSigmaPts_xmat .* exSigmaPt) * exSigmaPt';
S         = S + (wSigmaPts_zmat .* ezSigmaPt) * ezSigmaPt';
PxzPred   = PxzPred + exSigmaPt * (wSigmaPts_zmat .* ezSigmaPt)';

% ���㿨��������
K  = PxzPred / S;

% ������Ϣ
inovation = Z - zPred;

% ��״̬���� 
Xout = xPred + K*inovation;

% �Է������ 
Pout = PPred - K*S*K';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

