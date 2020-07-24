%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UT变换子函数
% 输入：fun为函数句柄，Xsigma为样本集，Wm和Wc为权值，n为状态维数，COV为协方差
% 输出：Xmeans为均值，P为协方差，Xdiv为预测减去均值
function [Xmeans,Xsigma_pre,P,Xdiv]=function_ut(Station,flag,Xsigma,Wm,Wc,n,COV)
LL=size(Xsigma,2);      % 得到Xsigma样本个数，size(Xsigma,2)的意思是返回矩阵Xsigma的列数
Xmeans=zeros(n,1);      % 均值
Xsigma_pre=zeros(n,LL); % 一步预测
for k=1:LL 
    % 这里需要判断一下应该调用哪一个函数
    if flag==1  % 调用ffun
        [y1,y2,y3,y4,y5,y6,y7,y8,y9] = feval('ffun',Xsigma(:,k));
        Xsigma_pre(:,k) = [y1,y2,y3,y4,y5,y6,y7,y8,y9]';  % 一步预测，这里的fun就是最前面的f或者h函数
    else % 调用hfun
        [dd,alpha,beta] = feval('hfun',Xsigma(:,k),Station);
        Xsigma_pre(:,k) = [dd,alpha,beta]';  % 一步预测，这里的fun就是最前面的f或者h函数
    end
    Xmeans=Xmeans+Wm(k)*Xsigma_pre(:,k);
end
% Xmeans(:,ones(1,LL))将Xmeans扩展成n*LL矩阵，每一列都相等
Xdiv=Xsigma_pre-Xmeans(:,ones(1,LL));  % 预测减去均值
P=Xdiv*diag(Wc)*Xdiv'+COV;             % 协方差