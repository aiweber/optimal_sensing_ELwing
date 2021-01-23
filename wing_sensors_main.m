%%
% This script finds optimal sensor locations using the SSPOC method. Run
% this script with optimal_sensing_ELwing as working directory.
%
% Steps:
%   1) Generate simulated strain data using Euler-Lagrange method
%   2) Transform simulated data with neural sensors
%   3) Dimensionality reduction on sensor data (PCA)
%   4) Optimization to find sensor locations
%
% Makes use of the cvx package, available from: http://cvxr.com/cvx/


%% set load, save, plot preferences

clearvars

simIter = 1; % number of times to run entire simulation (if noiseless spiking, only do 1)
plotFlag = 1; % plot results of each iteration (0 or 1)
dropout = 0;
nDraws = 5;  % how many times to do dropout on given dataset (only used if dropout = 1)

loadPars = 0; % 1 to load saved paramters, 0 to enter manually using function makeParameterStruct()
savePars = 0;
parFileName = 'WeberDefaultPars.mat'; %include .mat in file name, used for both loading and saving

loadResults = 0; % 1 to add to an existing results file, 0 to create new one
saveResults = 1; % 1 to save results to file, 0 does not save
resultsFileName = 'test_run'; %include .mat in file name, used for both loading and saving


%% set up to save results

if loadResults == 1
    % check if file exists, load file
    if isfile(resultsFileName)
        disp(['loading ', resultsFile])
        load(resultsFileName)
    else
        disp(['existing file named ' , resultsFileName,'.mat not found, creating new results table'])
        sparseSensorResults = table();
    end
else
    sparseSensorResults = table();
    if saveResults == 1
        disp(['Warning: This will overwrite any existing results file named ' resultsFileName '.mat'])
    end
end

%% initialize parameters

dataDir = pwd;
if loadPars == 1
    % check if file exists, load file
    
    assert(isfile(parFileName),'existing file named %s not found', parFileName)
    load(parFileName)
else
    Pars = makeParameterStruct();
    if savePars
        save(parFileName, 'Pars')
    end
end
ParCombo = parCombo(Pars); % set up for parameter sweep
parComboTable = struct2table(ParCombo,'AsArray',1);

%% Run simulation
nSensors = Pars.chordElements*Pars.spanElements;

for iPar = 1:length(ParCombo)

    % 1) Generate simulated strain data using Euler-Lagrange method
    StrainSet = eulerLagrangeSimWrapper(ParCombo(iPar),dataDir);
    
    % 2) transform simulated data with neural sensors & split into train and
    % test sets
    [X_ne, G_ne] = neuralTransformationOfData(StrainSet,ParCombo(iPar));
    
    for iter = 1:simIter
        
        timePtsPerSpikeRep = size(X_ne,2)/(ParCombo(iPar).sampFreq / ParCombo(iPar).flapFrequency)-length(unique(G_ne)); % first spike of each condition is removed
        X = zeros(size(X_ne,1),timePtsPerSpikeRep*ParCombo(iPar).spikeReps);
        G = zeros(1,timePtsPerSpikeRep*ParCombo(iPar).spikeReps);
        
        for spRep = 1:ParCombo(iPar).spikeReps
            thisRepIdx = (spRep-1)*timePtsPerSpikeRep+1:spRep*timePtsPerSpikeRep;
            [X(:,thisRepIdx), G(thisRepIdx), ~] = convertProbFiringToSpikes(X_ne,G_ne,ParCombo(iPar));
        end
        
        [X, G, XTest, GTest] = trainTestSplit(X,G,ParCombo(iPar));
        
        % 3) Dimensionality reduction on sensor data (PCA)
        
        %%% standardize training data
        Xmean = mean(X,2);
        Xstd = std(X,[],2);
        Xstd(Xstd<1e-14) = 1; % avoid divide by zero error
        XNorm = (X-Xmean)./repmat(Xstd,1,size(X,2));   
        [w_t, Psi] = dimReductionForSspoc(XNorm, G, ParCombo(iPar));
        
        % 4) Optimization to find sensor locations
        
        [sensors, cutoffLim, s] = sspocOptim(w_t, Psi, ParCombo(iPar), length(unique(G)));
        [~, I_top] = sort( sum(abs(s),2),'descend');
        
        % calculate classification accuracy
        sensorsSort = I_top(1:Pars.rmodes);
        if length(sensors)<5
            sensors = sensorsSort(1:5);
        end
        
        sensors10 = sensorsSort(1:10);
        
        acc10 = classAccuracyLin(X, G, XTest, GTest, sensors10);
        accAll = classAccuracyLin(X, G, XTest, GTest, 1:size(X,1));
        acc = classAccuracyLin(X, G, XTest, GTest, sensors);
                 
        if dropout
            dropoutAcc = zeros(9,nDraws);
            sensorsKeep = zeros(9,nDraws,9);
            for nKeep = 1:9  % drop all but 1 sensor
                for drawNum = 1:nDraws
                    sensorsTemp = sensors10(randperm(10)); % randomly shuffle top 10 sensors
                    sensorsTemp = sensorsTemp(1:nKeep);
                    dropoutAcc(nKeep,drawNum) = classAccuracyLin(X, G, XTest, GTest, sensorsTemp);
                    sensorsKeep(nKeep,drawNum,1:nKeep) = sensorsTemp;
                end
            end            
        end
        
        % Format results for saving
        results = formatResults(sensors, acc, sensors10, acc10, accAll, ParCombo(iPar));
        if dropout
            dropoutResults = table(mat2cell(dropoutAcc,size(dropoutAcc,1),size(dropoutAcc,2)),mat2cell(sensorsKeep,size(sensorsKeep,1),size(sensorsKeep,2),size(sensorsKeep,3)));
            results = [results dropoutResults];
        end
        sparseSensorResults = [sparseSensorResults; results]; %add new results to existing table
        
        if plotFlag
            plotSensorLocation(sensors10, acc10, Pars);
        end
        
        disp(['parameter set ' num2str(iPar) '/' num2str(length(ParCombo)) ', simIter ' num2str(iter) '/' num2str(simIter) ' done'])
    end %end parameter sweep loop
end %end iteration loop

%% save results to file

if saveResults == 1
    save(resultsFileName, 'sparseSensorResults')
    disp(['saved updated results in ', resultsFileName])
end