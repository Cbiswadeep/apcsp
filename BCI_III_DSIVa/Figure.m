load('data_set_IVa_al.mat')
cnt= 0.1*double(cnt);
cnt= eegButterFilter(cnt, 8, 30, 5);
[M_al,N_al]=size(cnt);
hil_al = hilbert(cnt);
phase_al= atan2(imag(hil_al), real(hil_al));
load('data_set_IVa_aa.mat')
cnt= 0.1*double(cnt);
cnt= eegButterFilter(cnt, 8, 30, 5);
[M_aa,N_aa]=size(cnt);
hil_aa = hilbert(cnt);
phase_aa= atan2(imag(hil_aa), real(hil_aa));
load('data_set_IVa_av.mat')
cnt= 0.1*double(cnt);
cnt= eegButterFilter(cnt, 8, 30, 5);
[M_av,N_av]=size(cnt);
hil_av = hilbert(cnt);
phase_av= atan2(imag(hil_av), real(hil_av));
load('data_set_IVa_aw.mat')
cnt= 0.1*double(cnt);
cnt= eegButterFilter(cnt, 8, 30, 5);
[M_aw,N_aw]=size(cnt);
hil_aw = hilbert(cnt);
phase_aw= atan2(imag(hil_aw), real(hil_aw));
load('data_set_IVa_ay.mat')
cnt= 0.1*double(cnt);
cnt= eegButterFilter(cnt, 8, 30, 5);
[M_ay,N_ay]=size(cnt);
hil_ay = hilbert(cnt);
phase_ay= atan2(imag(hil_ay), real(hil_ay));
M= [M_al M_aa M_av M_aw M_ay];
M_R= min(M);
% M_R= 1000;
rowToDelete = M_R:M_al;
phase_al(rowToDelete, :) = [];
rowToDelete = M_R:M_aa;
phase_aa(rowToDelete, :) = [];
rowToDelete = M_R:M_av;
phase_av(rowToDelete, :) = [];
rowToDelete = M_R:M_aw;
phase_aw(rowToDelete, :) = [];
rowToDelete = M_R:M_ay;
phase_ay(rowToDelete, :) = [];
plot(phase_al(:,1));
hold on;
plot(phase_aa(:,1));
plot(phase_av(:,1));
plot(phase_aw(:,1));
plot(phase_ay(:,1));