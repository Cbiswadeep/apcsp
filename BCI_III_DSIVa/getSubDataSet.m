function subDataSet = getSubDataSet(dataSet, trialIndexes)
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
%This function returns a subset of a given EEG, as specified by a list of
%trial indexes
%
%Input:
%dataSet: an EEG data set
%trialIndexes: an array with the indexes of the trials to be kept in the
%   sub data set
%
%Output
%subDataSet: an EEG data set containing only the selected trials 
%
%by Fabien LOTTE (fprlotte@i2r.a-star.edu.sg)
%created: 01/10/2009
%last revised: 01/10/2009

subDataSet.x = dataSet.x(:,:,trialIndexes);
subDataSet.y = dataSet.y(trialIndexes);
subDataSet.c = dataSet.c;
subDataSet.s = dataSet.s;