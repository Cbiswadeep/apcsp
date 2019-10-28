load('trainingEEGSignals.mat')

[M,N,P]=size(trainingEEGSignals{1,1}.x);
for j= 1:P
for i=1:N
    phase(:,i,j) = hilbert(trainingEEGSignals{1,1}.x(:,i));
end
end
phase1= atan2(imag(phase), real(phase));

amp_phase1 = (trainingEEGSignals{1,1}.x).* exp(1i*phase1);


CSP_amp = learnCSPLagrangian(trainingEEGSignals{1,1});

CSP_phase = learnCSPLagrangianp(trainingEEGSignals{1,1}, phase1);

CSP_phamp = learnCSPLagrangianphaamp(trainingEEGSignals{1,1},amp_phase1);
 
% acc1 = evaluateRCSPAlgo('CSP', 'BCI_III_DSIVa', 0)
% 
% acc2 = evaluateRCSPAlgo('CSPPhase', 'BCI_III_DSIVa', 0)
% 
% acc3 = evaluateRCSPAlgo('CSPPhaseAmp', 'BCI_III_DSIVa', 0)
%
%acc4 = evaluateRCSPAlgo('TR_CSP', 'BCI_III_DSIVa', 0)
amp=[];
phase86=[];
ampphase=[];
schur=[];
x_count = 0;
y_count = 0;
z_count = 0;
w_count = 0;
for i=1:10
    x = evaluateRCSPAlgo_var('CSP', 'BCI_III_DSIVa', 1,i);
    y = evaluateRCSPAlgo_var('CSPPhase', 'BCI_III_DSIVa', 0,i);
    z = evaluateRCSPAlgo_var('CSPPhaseAmp', 'BCI_III_DSIVa', 0,i);
    w = evaluateRCSPAlgo_var('CSPPhaseAmpJ', 'BCI_III_DSIVa', 0,i);
    amp= [amp x];
    phase86= [phase86 y];
    ampphase= [ampphase z];
    schur = [schur w];
    if x>y & x>z & x>w
        x_count = x_count+1;
    elseif y>x & y>z &y>w
        y_count = y_count+1;
    elseif z>x & z>y & z>w
        z_count = z_count+1;
    elseif w>x & w>y & w>z
        w_count = w_count+1;
    end
end
x , y, z , w;
fprintf('CSP: %d.\n',x_count);
fprintf('CSP_Phase: %d.\n',y_count);
fprintf('CSP_PhaseAmp: %d.\n',z_count);
fprintf('CSP_PhaseAmp: %d.\n',w_count);
duration= 1:10;
plot( duration, amp, duration, phase86, 'b', duration, ampphase, 'g', duration, schur,'k');