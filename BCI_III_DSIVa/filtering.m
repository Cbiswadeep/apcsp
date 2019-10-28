% This function returns the filtered output of a given (input) signal
% Fs = 250 Hz, Fpass1 = 3 Hz, Fpass2 = 13 Hz, Apass = 1, Order = 10

function [y] = filtering(x)
%x=data';
BandPassSpecObj = fdesign.bandpass('N,Fp1,Fp2,Ap', 10, 8, 12, 1, 500);
BandPassFilt = design(BandPassSpecObj, 'cheby1');
y = filter(BandPassFilt, x);
plot(y)
hold on
return


