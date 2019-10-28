% close all;
% clear all;
clc;
load('stimuli2.mat');
x2=stimuli2;
N=size(x2,2);
for i=1:N
    x2(i,:)=filtering(x2(i,:));
    [psd1 freq]=pwelch(x1(i,:));
    output(:,i)=psd1(:,1);
    [a1, b1] = dwt(x1(i,:),'db4');
    cA1(i,:) = a1;
    cD1(i,:) = b1;
    output(:,i)=cD1(i,:);
    [hjorth_params(i,1) hjorth_params(i,2) hjorth_params(i,3)] = hjorth(x1(i,:));
    ar_coeffs1(i,:) = aryule(x1(i,:),10);
    
end;
