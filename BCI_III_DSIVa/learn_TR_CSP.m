function CSPMatrix = learn_TR_CSP(EEGSignals,amp_phase1, alpha)
%
% <The RCSP (Regularized Common Spatial Pattern) Toolbox.>
%     Copyright (C) 2010  Fabien LOTTE
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%this function learn a Regularized CSP (Common Spatial Patterns) filters, using Tikhonov Regularization,
%to discriminate two mental states in EEG signals.
%
%Input:
%EEGSignals: the training EEG signals, composed of 2 classes. These signals
%are a structure such that:
%   EEGSignals.x: the EEG signals as a [Ns * Nc * Nt] Matrix where
%       Ns: number of EEG samples per trial
%       Nc: number of channels (EEG electrodes)
%       nT: number of trials
%   EEGSignals.y: a [1 * Nt] vector containing the class labels for each trial
%   EEGSignals.s: the sampling frequency (in Hz)
%alpha: regularization parameter (must be >= 0)
%
%Output:
%CSPMatrix: the learnt CSP filters (a [(2*Nc) * Nc] matrix with the spatial filters as rows)
%
%by Fabien LOTTE (fprlotte@i2r.a-star.edu.sg)
%created: 03/03/2010
%last revised: 03/03/2010
%
%See also: extractCSPFeatures

%check and initializations
nbChannels = size(amp_phase1,2);
nbTrials = size(amp_phase1,3);
classLabels = unique(EEGSignals.y);
nbClasses = length(classLabels);
if nbClasses ~= 2
    disp('ERROR! CSP can only be used for two classes');
    return;
end
covMatrices = cell(nbClasses,1); %the covariance matrices for each class

%computing the normalized covariance matrices for each trial
trialCov = zeros(nbChannels,nbChannels,nbTrials);
for t=1:nbTrials
    E = EEGSignals.x(:,:,t)';
    EE = E * E';
    trialCov(:,:,t) = EE ./ trace(EE);
end
clear E;
clear EE;

%computing the covariance matrix for each class
for c=1:nbClasses      
    covMatrices{c} = mean(trialCov(:,:,EEGSignals.y == classLabels(c)),3);  
end

%obtaining the first set of filters (class 1 / class 2)
M1 = inv(covMatrices{2} + alpha*eye(size(covMatrices{2}))) * covMatrices{1};
%eigen value decomposition of matrix M1
[U1 D1] = eig(M1);
eigenvalues = diag(D1);
[~, egIndex] = sort(eigenvalues, 'descend');
U1 = U1(:,egIndex);

%obtaining the second set of filters (class 2 / class 1)
M2 = inv(covMatrices{1} + alpha*eye(size(covMatrices{1}))) * covMatrices{2};
[U2 D2] = eig(M2);
eigenvalues = diag(D2);
[~, egIndex] = sort(eigenvalues, 'ascend');
U2 = U2(:,egIndex);

CSPMatrix = [U1 U2]';