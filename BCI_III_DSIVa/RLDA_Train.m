function ldaParams = RLDA_Train(trainingData)
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
%function to train a LRLDA (Ledoix's Regularized linear discriminant analysis) classifier using a
%set of training data for two classes. It is basically the classical LDA
%for which the covariance matrix is regularized according to Ledoix and
%Wolf' method (2004):
%O. Ledoit and M.Wolf. A well-conditioned estimator for large-dimensional covariance matrices. Journal
%of Multivariate Analysis, 88(2):365–411, 2004.
%
%Input:
%trainingData: matrix of feature vectors for training the LDA
%this matrix is a matrix with size [nT,nF+1] where
%          nT          = number of training data,
%          nF          = number of features (input dimension),
%          lastcol     = class labels (one per training data) 
%
%Output:
%ldaParams: the parameters of the LDA discriminant hyperplane with
%   ldaParams.a0: bias of the discriminant hyperplane    
%   ldaParams.a1N: slope of the discriminant hyperplane
%   ldaParams.classLabels: class labels for this data set
%the decision function is then a0 + a1N' * v where v is the input feature vector
%
%by Fabien Lotte (fprlotte@i2r.a-star.edu.sg)
%created: 02/03/2010
%last revised: 02/03/2010
%
%See also: LDA_Test

%checking that the number of training classes is 2
if size(unique(trainingData(:,size(trainingData,2))),1) ~= 2
    display('ERROR: The LDA classifier can only deal with two classes !');
    return;
end

%dividing data into two classes (and removing the label)
ldaParams.classLabels = unique(trainingData(:,end));
class1Data = trainingData(trainingData(:,end)==ldaParams.classLabels(1),1:(end-1));
class2Data = trainingData(trainingData(:,end)==ldaParams.classLabels(2),1:(end-1));

%mean vector estimation for each class
mu1 = mean(class1Data);
mu2 = mean(class2Data);

%covariance matrix estimation (using Ledoit and wolf's regularized covariance
%matrix estimation)
sigma1 = cov1para(class1Data);
sigma2 = cov1para(class2Data);
sigma = (sigma1 + sigma2)/2;

%computing the discriminant hyperplane coefficients
sigmaInv = pinv(sigma);
a0 = - (1/2) * (mu1 + mu2) * sigmaInv * (mu1 - mu2)';
a1N = sigmaInv * (mu1 - mu2)';

ldaParams.a0 = a0;
ldaParams.a1N = a1N;