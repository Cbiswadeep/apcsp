function filteredEEGData = eegButterFilter(EEGData, low, high, order)
% eegButterFilter(testingEEGSignals{1,1}, 8, 30, 5);
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
%band-pass filter a set of EEG signals using a butterworth filter
%
%input params:
%EEGData: extracted EEG signals
%low:low cutoff frequency
%high: high cutoff frequency
%order: filter order
%
%output
%filteredEEGData: the EEG data band-pass filtered in the specified
%   frequency band
%
%by Fabien Lotte (fprlotte@i2r.a-star.edu.sg)
%created: 17/02/2009
%last revised: 17/02/2009

%identifying various constants
nbSamples = size(EEGData.x,1);
nbChannels = size(EEGData.x,2);
nbTrials = size(EEGData.x,3);

%preparing the output data
filteredEEGData.x = zeros(nbSamples, nbChannels, nbTrials);
filteredEEGData.c = EEGData.c;
filteredEEGData.s = EEGData.s;
filteredEEGData.y = EEGData.y;

%designing the butterworth band-pass filter 
lowFreq = low * (2/EEGData.s);
highFreq = high * (2/EEGData.s);
[B A] = butter(order, [lowFreq highFreq]);

%filtering all channels in this frequency band, for the training data
for i=1:nbTrials %all trials
    for j=1:nbChannels %all channels
        filteredEEGData.x(:,j,i) = filter(B,A,EEGData.x(:,j,i));
    end
end
   