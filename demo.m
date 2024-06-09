warning off
clear all;
addpath(genpath('./'));
addpath('LibADMM-master/proximal_operators');
ds = {'yale_newdouble'};
di = 1;
load([ds{di},'.mat']);

numOfViews = length(X);
numOfSamples = size(X{1}, 2);
cls_num = length(unique(gt));

for i=1:numOfViews
    X{i} = X{i}./repmat(sqrt(sum(X{i}.^2,1)),size(X{i},1),1);  %normalized
    [L(:,:,i)] = solveF(X{i});
end

paras.alpha=0.01;
paras.beta=2;
paras.omega = [10 20 15];  %1*V

[S,obj,resofeachiter]=SGTL(L,cls_num,gt,paras.alpha,paras.beta,paras.omega);
S=double(S);

res=[];
for retry=1:20
    results=clustering8(abs(S)+abs(S'),cls_num, gt);
    res = [res;results];
end
mres=mean(res)
