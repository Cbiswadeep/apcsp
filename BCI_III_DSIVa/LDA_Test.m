function result = LDA_Test(testData, ldaParams)
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
%Classify a set of data using a trained LDA classifier and provide the
%resulting score (LDA output) and classification accuracy (%)
%
%Input:
%testData: a matrix representing the feature vectors to be classified by LDA
%this matrix is a matrix with size [nT,nF+1] where
%          nT          = number of training data,
%          nF          = number of features (input dimension),
%          lastcol     = class labels (one per training data) 
%ldaParams: the LDA hyperplane params (obtained after training with LDA_Train)
%
%Output: 
%result: the classification results where:
%   result.output: the predicted classification output for each feature vector
%   result.classes: the predicted class label for each feature vector
%   result.accuracy: the accuracy obtained (%)
%
% by Fabien Lotte (fprlotte@i2r.a-star.edu.sg)
% created: 02/03/2010
%last revised: 02/03/2010
%
% see also: LDA_Train

nbData = size(testData,1);
nbFeaturesPlus1 = size(testData,2);
trueLabels = testData(:,nbFeaturesPlus1);
result.output = zeros(nbData,1);
result.classes = zeros(nbData,1);

%classifying the input data
for i=1:size(testData,1)
        inputVec = testData(i,1:(nbFeaturesPlus1-1))';
        score = ldaParams.a0 + ldaParams.a1N' * inputVec;
        if score >= 0       
           result.classes(i) = ldaParams.classLabels(1);
        else
           result.classes(i) = ldaParams.classLabels(2);
        end        
        result.output(i) = score;
end

%computing the classification accuracy
nbErrors = sum(trueLabels~=result.classes);
result.accuracy = ((nbData - nbErrors)/nbData) * 100;