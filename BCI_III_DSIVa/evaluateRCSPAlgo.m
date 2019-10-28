function accuracy  = evaluateRCSPAlgo(algo, dataSet, printTopo)
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
% Assess the accuracy of various Regularized CSP (RCSP) algorithms
% on data from two possible BCI competition data sets: BCI competition III
% (2005), data set IVa (foot and right hand motor imagery),
% BCI competition III, data set IIIa (left and right hand motor imagery),
% or BCI competition IV, data set IIa (left and right hand motor imagery)
% The method is evaluated on the available testing sets.
% evaluateRCSPAlgo('CSP', 'BCI_III_DSIVa', 0)
%Input:
%algo: a string describing the kind of CSP algorithm to use for evaluation.
%   This string can be:
%       'CSP': regular CSP algorithm 
%           (Lagrangian optimization assuming invertible covariance matrices)
%       'CSPJointDiag': regular CSP algorithm (optimization by Joint Diagonalization)
%       'SR_CSP': Spatially Regularized CSP
%       'DL_SR_CSP': Spatially Regularized CSP with Diagional Loading
%       'DL_CSP_auto': regularized CSP with automatique diagonal loading (Ledoit & Wolf method)
%       'DL_CSP': regularized CSP with Diagonal Loading (cross validation)
%       'DL_CSP_diff': regularized CSP with Diagonal Loading (cross
%          validation) but different regularization parameter value for each class
%       'TR_CSP': CSP with Tikhonov Regularization
%       'WTR_CSP': CSP with Weight Tikhonov Regularization
%       'DL_TR_CSP': CSP with Tikhonov Regularization and Diagonal Loading
%       'SSR_CSP': CSP with Subject-to-Subject Regularization (i.e., using
%          EEG data from other subjects to regularize the covariances matrices)
%       'GLR_CSP': CSP with Generic Learning Regularization
%       'GL_TR_CSP': Tikhonov Regularized CSP with Generic Learning
%       'TR_CCSP1": Tikhonov Regularized Composite CSP method 1
%       'CCSP1': Composite CSP method 1 (de-emphasize matrices made from
%           fewer trials)
%       'CCSP2': Composite CSP method 2 (weigh matrices according to the
%           KL divergence measure)
%dataSet: a string identifying the data set used for evaluation
%       'BCI_III_DSIVa': BCI competition III (2005), data set IVa (foot and
%           right hand motor imagery) - 5 subjects - 118 electrodes\
%       'BCI_III_DSIIIa': BCI competition III (2005), data set IIIa (left and
%           right hand motor imagery only) - 3 subjects - 60 electrodes
%       'BCI_IV_DSIIa': BCI competition IV, data set IIa (left and right
%           hand motor imagery only) - 9 subjects - 22 electrodes
%printTopo: value = 1: the topography for the first and last 3 fitlers will
%   be saved into jpeg and eps figures
%           value = 0: the topography will not be saved
%
%output: accuracy: a vector v for which v(i) is the accuracy obtained by
%   the method on data from subject i
%
%by Fabien Lotte (fprlotte@i2r.a-star.edu.sg)
%created: 29/12/09
%last revised: 28/10/2010

%adding the path to use the necessary functions
% addpath '..\CSP_Algorithms\'
% addpath '..\CSP_Algorithms\CSP\';
% addpath '..\LDA_Classification';
% addpath '..\Utilities';
% addpath '..\EEG_Data\ReadingDataCode';

%loading EEG data from files
display('reading EEG data...');
if strcmp(dataSet,'BCI_III_DSIVa')
    if exist('trainingEEGSignals.mat','file')==0
        read_BCI_III_DSIVa;
    end
    disp('reading BCI competition III, data set IVa');
    load 'trainingEEGSignals.mat'
    load 'testingEEGSignals.mat'
    load 'elecCoord118.mat' %the electrodes 3D coordinates
    emap = load('EMapAll.mat');
% elseif strcmp(dataSet,'BCI_III_DSIIIa')
%     disp('reading BCI competition III, data set IIIa');
%     if exist('trainingEEGSignals.mat','file')==0
%         read_BCI_III_DSIIIa_LR;
%     end
%     load 'trainingEEGSignals.mat'
%     load 'testingEEGSignals.mat'
%     load 'elecCoord60.mat' %the electrodes 3D coordinates
%     emap = load('EMapAllEasyCap.mat');
% elseif strcmp(dataSet, 'BCI_IV_DSIIa')
%     disp('reading BCI competition IV, data set IIa');
%     if exist('..\EEG_Data\BCI_IV_DSIIa\trainingEEGSignals.mat','file')==0
%         read_BCI_IV_DSIIa([1 2])
%     end
%     load '..\EEG_Data\BCI_IV_DSIIa\trainingEEGSignals.mat'
%     load '..\EEG_Data\BCI_IV_DSIIa\testingEEGSignals.mat'  
%     load '..\EEG_Data\elecCoord22.mat' %the electrodes 3D coordinates
%     emap = load('..\EEG_Data\EMapAll.mat');
else
    disp('!! ERROR !! Incorrect data set! Possible sets are: BCI_III_DSIVa, BCI_III_DSIIIa and BCI_IV_DSIIa');    
    return
end       
emap = emap.EMapAll;

display('done !');

nbSubjects = length(trainingEEGSignals);
accuracy = zeros(1,nbSubjects);
nbFilterPairs = 3; %we use 3 pairs of filters

%parameters for band pass filtering the signals in the mu+beta band (8-30 Hz)
low = 8;
high = 30;
order = 100; %we use a 5th-order butterworth filter

%variables to compute the computation time
trainingTime = zeros(nbSubjects,1);

%hyperparameters of the RCSP approaches
alphaExp = 10:-1:1;
alphaList = 10.^-alphaExp; %regularization parameter of the objective function
rList = [0.01 0.05 0.1 0.5 0.8 1.0 1.2 1.5]; %parameter of the SR_CSP
gammaList = 0:0.1:0.9; %first regularization parameter for covariance matrix regularization
betaList = 0:0.1:0.9; %second regularization parameter for covariance matrix regularization

k = 10; %for hyperparameter selecting using k-fold cross validation

root = ['Pictures' algo '\' dataSet '\']; %to save the filter topography pictures

%filtering each trial in the mu+beta frequency band
%Note: Dieter Devlaminck made me realize (thank you Dieter!) that this
%filtering scheme was not optimal. It would be indeed better to filter the
%whole signal first, rather than each trial individually. I have added this
%possibility in the read_BCIXX_DSXX functions (see the function header). 
%The results presented in the IEEE-TBME paper are based on a filtering of each trial individually.
 for s=1:nbSubjects
     trainingEEGSignals{s} = eegButterFilter(trainingEEGSignals{s}, low, high, order);
     testingEEGSignals{s} = eegButterFilter(testingEEGSignals{s}, low, high, order);
 end

for s=1:nbSubjects    
    %Learning the Regularized CSP filters, depending on the algorithm chosen
    tic; %to evaluate the learning time
    
    otherSubjects = [1:(s-1) (s+1):nbSubjects];
    
    if strcmp(algo,'CSP')
        disp('Learning CSP filters assuming invertible covariance matrices');
        CSPMatrix = learnCSPLagrangian(trainingEEGSignals{s});
     elseif strcmp(algo,'CSPPhaseAmp')
         disp('Learning CSP filters with joint diagonalization');
         [M,N,P]=size(trainingEEGSignals{1,1}.x);
            for j= 1:P
            for i=1:N
                phase(:,i,j) = hilbert(trainingEEGSignals{1,1}.x(:,i));
            end
            end
            phase1= atan2(imag(phase), real(phase));
         amp_phase1= (trainingEEGSignals{1,1}.x).* exp(1i*phase1);
         CSPMatrix = learnCSPLagrangianphaamp(trainingEEGSignals{1,1},amp_phase1);
     elseif strcmp(algo, 'CSPPhase')
         disp('Learning CSP filters for phase');
         [M,N,P]=size(trainingEEGSignals{1,1}.x);
            for j= 1:P
            for i=1:N
                phase(:,i,j) = hilbert(trainingEEGSignals{1,1}.x(:,i));
            end
            end
            phase1= atan2(imag(phase), real(phase));
         CSPMatrix = learnCSPLagrangianp(trainingEEGSignals{1,1}, phase1);
%     elseif strcmp(algo, 'DL_CSP_diff')
%         disp('Learning CSP filters with Diagonal Loading and different regularization parameter values for each class');
%         bestGamma = DL_CSP_diffWithBestParams(trainingEEGSignals{s}, [gammaList;gammaList], k, nbFilterPairs);
%         disp('bestGamma: '); disp(bestGamma);
%         CSPMatrix = learn_DL_CSP_diff(trainingEEGSignals{s}, bestGamma);
%     elseif strcmp(algo, 'DL_CSP')
%         disp('Learning CSP filters with Diagonal Loading regularization (cross validation)');
%         bestGamma = DL_CSPWithBestParams(trainingEEGSignals{s}, gammaList, k, nbFilterPairs);
%         CSPMatrix = learn_DL_CSP(trainingEEGSignals{s}, bestGamma);
     elseif strcmp(algo,'TR_CSP')            
         disp('Learning Tikhonov Regularized CSP filters');    
         %step 1: find the optimal hyperparameters        
         bestAlpha = TR_CSPWithBestParams(trainingEEGSignals{1,1}, alphaList, k, nbFilterPairs); 
         %step 2: learn the TR CSP with the best hyperparameters
         CSPMatrix = learn_TR_CSP(trainingEEGSignals{1,1}, bestAlpha);
%     elseif strcmp(algo,'WTR_CSP')            
%         disp('Learning Weighted Tikhonov Regularized CSP filters');    
%         %learning the weight from other subjects        
%         weightVector = learnWTR_CSP_WeightVector(trainingEEGSignals(otherSubjects), nbFilterPairs);
%         %step 1: find the optimal hyperparameters        
%         bestAlpha = WTR_CSPWithBestParams(trainingEEGSignals{s}, weightVector, alphaList, k, nbFilterPairs);
%         disp(['BestAlpha: ' num2str(bestAlpha)]);
%         %step 2: learn the WTR CSP with the best hyperparameters
%         CSPMatrix = learn_WTR_CSP(trainingEEGSignals{s}, weightVector, bestAlpha);
%     elseif strcmp(algo,'DL_TR_CSP')
%         disp('Learning Tikhonov Regularized CSP filters with Diagonal Loading');
%         %step 1: find the optimal hyperparameters        
%         [bestAlpha bestGamma] = DL_TR_CSPWithBestParams(trainingEEGSignals{s}, alphaList, gammaList, k, nbFilterPairs);
%         disp(['BestAlpha: ' num2str(bestAlpha) ' - bestGamma: ' num2str(bestGamma)]);
%         %step 2: learn the DL TR CSP with the best hyperparameters
%         CSPMatrix = learn_DL_TR_CSP(trainingEEGSignals{s}, bestAlpha, bestGamma);
%     elseif strcmp(algo,'SR_CSP')        
%         disp('Learning Spatially Regularized CSP filters');    
%         %step 1: find the optimal hyperparameters        
%         [bestAlpha bestR] = SR_CSPWithBestParams(trainingEEGSignals{s}, elecCoord, alphaList, rList, k, nbFilterPairs);
%         disp(['Best r: ' num2str(bestR)]);
%         %step 2: learn the SR CSP with the best hyperparameters
%         CSPMatrix = learn_SR_CSP(trainingEEGSignals{s}, elecCoord, bestAlpha, bestR);
%     elseif strcmp(algo,'DL_SR_CSP')
%         disp('Learning Spatially Regularized CSP filters with Diagonal Loading');    
%         %step 1: find the optimal hyperparameters        
%         [bestAlpha bestGamma bestR] = DL_SR_CSPWithBestParams(trainingEEGSignals{s}, elecCoord, alphaList, gammaList, rList, k, nbFilterPairs);
%         disp(['Best r: ' num2str(bestR) ' - bestAlpha: ' num2str(bestAlpha) ' - bestGamma: ' num2str(bestGamma)]);
%         %step 2: learn the SR CSP with the best hyperparameters
%         CSPMatrix = learn_DL_SR_CSP(trainingEEGSignals{s}, elecCoord, bestAlpha, bestGamma, bestR);
%     elseif strcmp(algo,'SSR_CSP')
%         disp('Learning Subjects Subset Regularized CSP filters');         
%         otherTrainingEEGSignals = trainingEEGSignals(otherSubjects);        
%         %step 1: identifying the subset of subjects to use
%         bestSubset = SSR_CSP_SelectSubjects(trainingEEGSignals{s}, otherTrainingEEGSignals, nbFilterPairs);
%         %step 2: identifying the best regularization parameter
%         bestBeta = SSR_CSPWithBestParams(trainingEEGSignals{s}, otherTrainingEEGSignals(bestSubset), betaList, k, nbFilterPairs);
%         %step 3: learning the SSR CSP with the best hyperparameter
%         CSPMatrix = learn_SSR_CSP_fixedSelection(trainingEEGSignals{s}, otherTrainingEEGSignals(bestSubset), bestBeta);
%         clear otherTrainingEEGSignals;
%     elseif strcmp(algo,'GLR_CSP')
%         disp('Learning Generic Learning Regularized CSP filters');         
%         [bestGamma bestBeta] = GLR_CSPWithBestParams(trainingEEGSignals{s}, trainingEEGSignals(otherSubjects), gammaList, betaList, k, nbFilterPairs);
%         CSPMatrix = learn_GLR_CSP(trainingEEGSignals{s}, trainingEEGSignals(otherSubjects), bestGamma, bestBeta);
%     elseif strcmp(algo,'GL_TR_CSP')
%         disp('Learning Generic Learning Tikhonov Regularized CSP filters');         
%         [bestAlpha bestBeta] = GL_TR_CSPWithBestParams(trainingEEGSignals{s}, trainingEEGSignals(otherSubjects), alphaList, betaList, k, nbFilterPairs);
%         disp(['BestAlpha: ' num2str(bestAlpha) ' - bestBeta: ' num2str(bestBeta)]);
%         CSPMatrix = learn_GLR_CSP(trainingEEGSignals{s}, trainingEEGSignals(otherSubjects), bestAlpha, bestBeta);
%     elseif strcmp(algo,'TR_CCSP1')
%         disp('Learning Tikhonov Regularized Composite CSP method 1 filters');         
%         [bestAlpha bestBeta] = TR_CCSP1WithBestParams(trainingEEGSignals{s}, trainingEEGSignals(otherSubjects), alphaList, betaList, k, nbFilterPairs);
%         disp(['BestAlpha: ' num2str(bestAlpha) ' - bestBeta: ' num2str(bestBeta)]);
%         CSPMatrix = learn_TR_CCSP1(trainingEEGSignals{s}, trainingEEGSignals(otherSubjects), bestAlpha, bestBeta);
%     elseif strcmp(algo, 'CCSP1')
%         disp('Learning Composite CSP method 1');            
%         bestBeta = CCSP1WithBestParams(trainingEEGSignals{s}, trainingEEGSignals(otherSubjects), betaList, k, nbFilterPairs);
%         CSPMatrix = learn_CCSP1(trainingEEGSignals{s}, trainingEEGSignals(otherSubjects), bestBeta);
%     elseif strcmp(algo, 'CCSP2')
%         disp('Learning Composite CSP method 2');        
%         bestBeta = CCSP2WithBestParams(trainingEEGSignals{s}, trainingEEGSignals(otherSubjects), betaList, k, nbFilterPairs);
%         CSPMatrix = learn_CCSP2(trainingEEGSignals{s}, trainingEEGSignals(otherSubjects), bestBeta);
    else
        disp('!! ERROR !! Incorrect CSP algorithm !');
        disp('Possible algorithms are: CSP, CSPJointDiag, DL_CSP_auto, DL_CSP, TR_CSP, WTR_CSP, DL_TR_CSP, SR_CSP, DL_SR_CSP, SSR_CSP, GLR_CSP, GL_TR_CSP, TR_CCSP1, CCSP1, CCSP2');
        return;
    end
    trainingTime(s) = toc;
    
    %plotting the best CSP filter, if needed
    if printTopo==1
        %creating the map of channel location                
        elecIndex = ismember(emap(1,:),trainingEEGSignals{s}.c);    
        emap = emap(:,elecIndex);        
        for f=1:nbFilterPairs
            plotElecPotentials(emap,CSPMatrix(f,:));
%             filenameFig = ['' root 'S' num2str(s) '_FilterPair' num2str(f) '_C1.fig'];
%             filenameEPS = ['' root 'S' num2str(s) '_FilterPair' num2str(f) '_C1.eps'];
            filenameJPG = ['' root 'S' num2str(s) '_FilterPair' num2str(f) '_C1.jpeg'];            
            print(filenameFig);
            print('-depsc',filenameEPS);
            print('-djpeg',filenameJPG);
            close;

            plotElecPotentials(emap,CSPMatrix(end-f+1,:));
%             filenameFig = ['' root 'S' num2str(s) '_FilterPair' num2str(f) '_C2.fig'];
%             filenameEPS = ['' root 'S' num2str(s) '_FilterPair' num2str(f) '_C2.eps'];
            filenameJPG = ['' root 'S' num2str(s) '_FilterPair' num2str(f) '_C2.jpeg'];
            print(filenameFig);
            print('-depsc',filenameEPS);
            print('-djpeg',filenameJPG);
            close;
        end
    end
        
    %extracting CSP features from the training set
    trainFeatures = extractCSPFeatures(trainingEEGSignals{s}, CSPMatrix, nbFilterPairs);
    
    %training a LDA classifier on these features
    ldaParams = LDA_Train(trainFeatures);
    
    %extracting CSP features from the testing set
    testFeatures = extractCSPFeatures(testingEEGSignals{s}, CSPMatrix, nbFilterPairs);
    
    %classifying the test features with the learnt LDA
    result = LDA_Test(testFeatures, ldaParams);    
    accuracy(s) = result.accuracy;
    disp(['test set accuracy for subject' num2str(s) ' = ' num2str(accuracy(s)) ' %']);    
end

%computing average time needed for training
meanTrain = mean(trainingTime); stdTrain = std(trainingTime);
disp(['Training time: ' num2str(meanTrain) ' - std: ' num2str(stdTrain)]);
%accuracyFilename = ['..\GeneratedData\accuracy_' algo '-' dataSet '.mat'];
save(['accuracy_' algo '-' dataSet '.mat'],'accuracy');