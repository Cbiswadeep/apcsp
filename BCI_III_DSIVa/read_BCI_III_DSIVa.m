function read_BCI_III_DSIVa(passBand)
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
%read data sets from BCI competition III data set IVa
%and convert them to the I2R format
%
%Input:
%passBand: an optional structure defining a frequency band in which
%   filtering all the EEG signals this structure is such that:
%       passBand.low= low cut-off frequency (in Hz)
%       passBand.high= high cut-off frequency (in Hz)
%   by default no filtering is done
%
%This file is supposed to be used from the Evaluations directory
%
%by Fabien Lotte (fprlotte@i2r.a-star.edu.sg)

%names of the files with the original data
filenames = ['data_set_IVa_aa.mat'];
         
trueLabelsFiles = ['true_labels_aa.mat'];
         
nbSubjects = 1;
trainingEEGSignals = cell(nbSubjects,1);
testingEEGSignals = cell(nbSubjects,1);

%some constants
fs = 100; %sampling rate
startEpoch = 0.5; %an epoch starts 0.5s after the cue
endEpoch = 2.5;%an epoch ends 2.5s after the cue
nbSamplesPerTrial = ceil((endEpoch - startEpoch) * fs) + 1;

for s=1:nbSubjects    
    
    disp(['Reading data from subject ' num2str(s)]);
    
    %reading the data from this subject
    disp('reading files...');
    load(filenames(s,:));
    load(trueLabelsFiles(s,:));
    disp('...done!');
    
    %conversion to uV values
    cnt= 0.1*double(cnt);
    
    %perform band-pass filtering with a butterworth filter if needed
    if exist('passBand','var')
        disp(['will band-pass filter all signals in ' num2str(passBand.low) '-' num2str(passBand.high) 'Hz']); 
        order = 5; 
        lowFreq = passBand.low * (2/fs);
        highFreq = passBand.high * (2/fs);
        [B A] = butter(order, [lowFreq highFreq]);
        cnt = filter(B,A,cnt);
    end
    
    nbChannels = size(cnt,2);
    
    %identifying the training and testing trials
    labels = mrk.y;
    cues = mrk.pos;
    trueLabels = true_y;
    trainingIndexes = find(~isnan(labels));
    testingIndexes = find(isnan(labels));
    
    nbTrainTrials = length(trainingIndexes);
    disp(['nbTrainTrials = ' num2str(nbTrainTrials)]);
    nbTestTrials = length(testingIndexes);
    disp(['nbTestTrials =  ' num2str(nbTestTrials)]);
    
    %initializing structures
    disp('initializing structures...');
    trainingEEGSignals{s}.x = zeros(nbSamplesPerTrial, nbChannels, nbTrainTrials);
    trainingEEGSignals{s}.y = labels(trainingIndexes)-1;
    trainingEEGSignals{s}.s = fs;
    trainingEEGSignals{s}.c = nfo.clab;
    testingEEGSignals{s}.x = zeros(nbSamplesPerTrial, nbChannels, nbTestTrials);
    testingEEGSignals{s}.y = trueLabels(testingIndexes)-1;
    testingEEGSignals{s}.s = fs;
    testingEEGSignals{s}.c = nfo.clab;
    disp('...done!');
    
    %assigning data to the corresponding structure
    disp('assigning data to the corresponding structure...');
    
    %training set    
    for trial=1:nbTrainTrials    
        %getting the cue
        cueIndex = cues(trainingIndexes(trial));
        %getting the data
        epoch = cnt((cueIndex + round(startEpoch*fs)):(cueIndex + round(endEpoch*fs)),:);
        %disp(size(epoch));
        %disp(size(trainingEEGSignals{s}.x(:,:,trial)));
        trainingEEGSignals{s}.x(:,:,trial) = epoch;
    end
    
    %testing set
    for trial=1:nbTestTrials    
        %getting the cue
        cueIndex = cues(testingIndexes(trial));
        %getting the data
        testingEEGSignals{s}.x(:,:,trial) = cnt((cueIndex + round(startEpoch*fs)):(cueIndex + round(endEpoch*fs)),:);
    end
    
    disp('...done!');
end

%saving the results to the appropriate matlab files
save('trainingEEGSignals.mat','trainingEEGSignals');
save('testingEEGSignals.mat','testingEEGSignals');
    
    
    
    
    
    