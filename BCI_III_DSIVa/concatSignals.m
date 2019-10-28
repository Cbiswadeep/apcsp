function concatenatedSignals = concatSignals(eegDataSet)
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
%concatenate the input EEG data sets into a single, bigger data set
%
%input:
%eegDataSet: a cell array, with each cell containing an eegDataSet
%
%output:
%concatenatedSignals: the resulting eeg data set, being formed by the
%   concatenation of all input EEG data sets
%
%by Fabien LOTTE (fprlotte@i2r.a-star.edu.sg)
%created: 28/05/2009
%last revised: 28/05/2009

concatenatedSignals = eegDataSet{1};

for d=2:length(eegDataSet)
    concatenatedSignals.x = cat(3,concatenatedSignals.x, eegDataSet{d}.x);
    concatenatedSignals.y = [concatenatedSignals.y eegDataSet{d}.y];
end
