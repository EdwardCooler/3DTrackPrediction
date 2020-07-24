%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 无迹卡尔曼滤波算法
function [Xout,Pout]=ukf(Xin,Z,Pin,Qukf,Rukf,t) 
% 无迹变换的参数
alpha = 1;
beta  = 0;
kappa = 2;

states       = size(Xin(:),1);
observations = size(Z(:),1);
vNoise       = size(Qukf,2);
wNoise       = size(Rukf,2);
noises       = vNoise+wNoise;

% 用噪声向量对状态向量进行扩展 
% 注意，此处简化为加性噪声模型，为的是节约计算开支
if (noises)
    N=[Qukf zeros(vNoise,wNoise); zeros(wNoise,vNoise) Rukf];
    PQ=[Pin zeros(states,noises);zeros(noises,states) N];
    xQ=[Xin;zeros(noises,1)];
else
    PQ=Pin;
    xQ=Xin;
end

% 通过UT变换，计算sigma点集和它们的权值
[xSigmaPts, wSigmaPts, nsp] = scaledSymmetricSigmaPoints(xQ, PQ, alpha, beta, kappa);
 
% 为了加速运行（代码执行效率），将wSigmaPts复制到矩阵中
wSigmaPts_xmat = repmat(wSigmaPts(:,2:nsp),states,1);
wSigmaPts_zmat = repmat(wSigmaPts(:,2:nsp),observations,1);

% 计算sigma点集预测值
xPredSigmaPts = feval('ffun',xSigmaPts(1:states,:),t)+xSigmaPts(states+1:states+vNoise,:);
zPredSigmaPts = feval('hfun',xPredSigmaPts,t)+xSigmaPts(states+vNoise+1:states+noises,:);

% 计算均值的预测值
xPred = sum(wSigmaPts_xmat .* (xPredSigmaPts(:,2:nsp) - repmat(xPredSigmaPts(:,1),1,nsp-1)),2);
zPred = sum(wSigmaPts_zmat .* (zPredSigmaPts(:,2:nsp) - repmat(zPredSigmaPts(:,1),1,nsp-1)),2);
xPred=xPred+xPredSigmaPts(:,1);
zPred=zPred+zPredSigmaPts(:,1);
 
% 计算方差和协方差，注意第一列的权值与均值的权是不一样的
% 此处主要是做无迹变化的工作
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

% 计算卡尔曼增益
K  = PxzPred / S;

% 计算新息
inovation = Z - zPred;

% 对状态更新 
Xout = xPred + K*inovation;

% 对方差更新 
Pout = PPred - K*S*K';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

