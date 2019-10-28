function pooledFeatureSet = concatFeatureSets(featureSets)
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
%concatenante the feature vector sets contained in the cell array of data
%set passed as parameter into a single feature vector set
%
%Input:
%featureSets: a cell array with a feature vector set per cell
%
%Output:
%pooledFeatureSet: the concatenation of all feature vectors sets contained
%in the input parameter
%
%by Fabien LOTTE (fprlotte@i2r.a-star.edu.sg)
%created: 17/06/2009
%last revised: 17/06/2009

nbDataSets = length(featureSets);

%computing the size of the concatenation
nbVectors = 0;
for d=1:nbDataSets
    nbVectors = nbVectors + size(featureSets{d},1);
end

pooledFeatureSet = zeros(nbVectors, size(featureSets{1},2));

currentIndex = 1;
for d=1:nbDataSets
    nbVectorLocal = size(featureSets{d},1);
    pooledFeatureSet(currentIndex:(currentIndex+nbVectorLocal-1),:) = featureSets{d};
    currentIndex = currentIndex+nbVectorLocal;
end

