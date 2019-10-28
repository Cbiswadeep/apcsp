load('data_set_IVa_al.mat')
cnt= 0.1*double(cnt);
[M,N]=size(cnt);
phase=[];
for i=1:N
    phase(:,i) = hilbert(cnt(:,i));
end 
phase1= atan2(imag(phase), real(phase));
X=phase1;
% plot(X(:,1),X(:,2), '.');
plot(X(:,1));
