
[M,N,P]=size(trainingEEGSignals{1,1}.x);
for j= 1:P
for i=1:N
    phase(:,i,j) = hilbert(trainingEEGSignals{1,1}.x(:,i));
end
end
phase1= atan2(imag(phase), real(phase));
phase1=reshape(phase1,[],3);
X=phase1;
plot3(X(:,1),X(:,2),X(:,3),'.');