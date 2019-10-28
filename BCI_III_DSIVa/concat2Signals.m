function concatenatedSignals = concat2Signals(eegDataSet_A, eegDataSet_B)
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
%concatenate the 2 input EEG data sets into a single, bigger data set
%
%input:
%eegDataSet_A: the first EEG data set
%eegDataSet_N: the second EEG data set
%
%output:
%concatenatedSignals: the resulting eeg data set, being formed by the
%   concatenation of the 2 input EEG data sets
%
%by Fabien LOTTE (fprlotte@i2r.a-star.edu.sg)
%created: 28/05/2009
%last revised: 01/10/2009

concatenatedSignals.c = eegDataSet_A.c;
concatenatedSignals.s = eegDataSet_A.s;
concatenatedSignals.x = cat(3,eegDataSet_A.x, eegDataSet_B.x);
concatenatedSignals.y = [eegDataSet_A.y eegDataSet_B.y];