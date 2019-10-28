function [bestAlpha bestScore] = TR_CSPWithBestParams(trainingEEGSignals, alphaList, k, nbFilterPairs)
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
%identify the best hyperparameter for the TR_CSP(Common Spatial Pattern with Tikhonov Regularization) spatial filter.
%This hyperparameter (regularization parameter Alpha) is
%selected using K-fold stratified Cross Validation (CV) and a LDA (Linear Discriminant
%Analysis) as classifier
%
%Input:
%trainingEEGSignals: all training EEG signals. It is an EEG data set
%   structure such that:
%   EEGSignals.x: the EEG signals as a [Ns * Nc * Nt] Matrix where
%       Ns: number of EEG samples per trial
%       Nc: number of channels (EEG electrodes)
%       nT: number of trials
%   EEGSignals.y: a [1 * Nt] vector containing the class labels for each trial
%   EEGSignals.s: the sampling frequency (in Hz)
%alphaList: the list (vector) of candidate values for the regularization parameter
%k: number of fold in the k-fold cross-validation
%nbFilterPairs: the number of filter pairs to be considered for feature
%   extraction
%
%Output:
%bestAlpha: the best value for the regularization parameter
%bestScore: the k-fold CV accuracy obtained using the best hyperparameters
%CSPMatrix: the CSP Matrix obtained for the best hyperparameters 
%   (the spatial filters are the rows of this matrix)
%
%by Fabien LOTTE (fprlotte@i2r.a-star.edu.sg)
%created: 10/03/2009
%last revised: 10/03/2010

%identifying the k subtraining set and validation set
disp('generating the subtraining and validation sets for each fold...');
subTrainingSets = cell(k,1);
validationSets = cell(k,1);
classLabels = unique(trainingEEGSignals.y);
%dividing the data according to each class
trainingSetClass1 = getSubDataSet(trainingEEGSignals, trainingEEGSignals.y == classLabels(1));
trainingSetClass2 = getSubDataSet(trainingEEGSignals, trainingEEGSignals.y == classLabels(2));
nbTrials1 = length(trainingSetClass1.y);
nbTrials2 = length(trainingSetClass2.y);
sizeChunk1 = floor(nbTrials1/k);
sizeChunk2 = floor(nbTrials2/k);
clear trainingEEGSignals;

%generating the different training/testing sets of the cross validation
for iter=1:k            
    subTrainingSetClass1 = getSubDataSet(trainingSetClass1, [1:(iter-1)*sizeChunk1 (iter*sizeChunk1+1):nbTrials1]);
    subTrainingSetClass2 = getSubDataSet(trainingSetClass2, [1:(iter-1)*sizeChunk2 (iter*sizeChunk2+1):nbTrials2]);
    validationSetClass1 = getSubDataSet(trainingSetClass1, ((iter-1)*sizeChunk1+1):(iter*sizeChunk1));
    validationSetClass2 = getSubDataSet(trainingSetClass2, ((iter-1)*sizeChunk2+1):(iter*sizeChunk2));
    subTrainingSets{iter} = concat2Signals(subTrainingSetClass1, subTrainingSetClass2);
    clear subTrainingSetClass1; clear subTrainingSetClass2;
    validationSets{iter} = concat2Signals(validationSetClass1, validationSetClass2);
    clear validationSetClass1; clear validationSetClass2;
end

%evaluating the performances of the difference potential hyperparameters
bestScore = 0;

%computing the k-fold cross validation accuracy for each hyperparameter
%values using a LDA as classifier
for alpha = alphaList        
    score = 0;
    for iter=1:k
        %learning CSP spatial filters
        CSPMatrix = learn_TR_CSP(subTrainingSets{iter}, alpha);
        featureTrain = extractCSPFeatures(subTrainingSets{iter}, CSPMatrix, nbFilterPairs);
        LDAModel = LDA_Train(featureTrain);            
        featureTest = extractCSPFeatures(validationSets{iter}, CSPMatrix, nbFilterPairs);
        results = LDA_Test(featureTest, LDAModel);
        localScore = results.accuracy;
        score = score + localScore;
    end
    score = score / k;  
    disp(['alpha=' num2str(alpha) ' => ' num2str(score) '%']);
    if score > bestScore
        bestScore = score;
        bestAlpha = alpha;
    end    
end         
disp(['bestAlpha=' num2str(bestAlpha)]);