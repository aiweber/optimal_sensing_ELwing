function [X, G, XTest, GTest] = trainTestSplit(X,G,SplitDataPrs)
% [X, G, XTest, GTest] = trainTestSplit(X,G,SplitDataPrs)
%
% Randomly splits data into train (X,G) and test (XTest,GTest) datasets,
% using SplitDataPrs.trainFraction to determine fraction of data that will
% be used for training

nDataPts = size(X,2);
idx = randperm(nDataPts);
trainTestCutoff = round(nDataPts*SplitDataPrs.trainFraction);

XTest = X(:,idx(trainTestCutoff+1:end));
GTest = G(idx(trainTestCutoff+1:end));
X(:,idx(trainTestCutoff+1:end)) = [];
G(idx(trainTestCutoff+1:end)) = [];