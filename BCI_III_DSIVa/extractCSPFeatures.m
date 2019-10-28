function features = extractCSPFeatures(EEGSignals, CSPMatrix, nbFilterPairs)
% extractCSPFeatures('BCI_III_DSIVa','CSPPhaseAmp', 3)
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
%extract features from an EEG data set using the Common Spatial Patterns (CSP) algorithm
%
%Input:
%EEGSignals: the EEGSignals from which extracting the CSP features. These signals
%are a structure such that:
%   EEGSignals.x: the EEG signals as a [Ns * Nc * Nt] Matrix where
%       Ns: number of EEG samples per trial
%       Nc: number of channels (EEG electrodes)
%       nT: number of trials
%   EEGSignals.y: a [1 * Nt] vector containing the class labels for each trial
%   EEGSignals.s: the sampling frequency (in Hz)
%CSPMatrix: the CSP projection matrix, learnt previously (see function learnCSP)
%nbFilterPairs: number of pairs of CSP filters to be used. The number of
%   features extracted will be twice the value of this parameter. The
%   filters selected are the one corresponding to the lowest and highest
%   eigenvalues
%
%Output:
%features: the features extracted from this EEG data set 
%   as a [Nt * (nbFilterPairs*2 + 1)] matrix, with the class labels as the
%   last column   
%
%by Fabien LOTTE (fprlotte@i2r.a-star.edu.sg)
%created: 02/03/2010
%last revised: 02/03/2010

%initializations
nbTrials = size(EEGSignals.x,3);
features = zeros(nbTrials, 2*nbFilterPairs+1);
Filter = CSPMatrix([1:nbFilterPairs (end-nbFilterPairs+1):end],:);

%extracting the CSP features from each trial
for t=1:nbTrials    
    %projecting the data onto the CSP filters    
    projectedTrial = Filter * EEGSignals.x(:,:,t)';    
    
    %generating the features as the log variance of the projected signals
    variances = var(projectedTrial,0,2);    
    for f=1:length(variances)
        features(t,f) = log(variances(f));
    end
    features(t,end) = EEGSignals.y(t);    
end