function [X_spt,G_spt,Xspikes] = convertProbFiringToSpikes(X,G,ParStruct)
% function [X_spt,G_spt,Xspikes] = convertProbFiringToSpikes(X,G,NeuralPrs)
%
% This function converts probability of firing data (output of
% "neuralTransformationOfData" function) into time-to-first-spike data.
% Note that this function creates a matrix of the time to first spike
% *within each wingstroke* at each sensor location.  The number of data
% points in the data matrix will be greatly reduced (e.g., if there are 40
% time points per wingstroke in the original X matrix, these 40 points will 
% be reduced to a single data point in the X_spt matrix).
%
% The function inputs are
%   X: original data matrix, sensor locations x time points
%   G: class vector for matrix X
%   ParStruct: parameter structure which contains
%       refPer: refractory period, ms
%       sampFreq: sampling frequency, points/sec
%       flapFrequency: frequency of flapping, flaps/sec
%
% The function outputs are
%   X_spt: new data matrix of times to first spike, sensor locations x wingstroke number
%   G_spt: new matrix of classes for X_spt
%   Xspikes: matrix of spiking data (0s and 1s) at original time resolution

refPerInTPts = round(ParStruct.refPer*ParStruct.sampFreq/1000);
maxT = size(X,2);

Xspikes = zeros(size(X));
Xtemp = X;  % running calculation of probability of spiking with refractory period incorporated
for locNum = 1:size(X,1)  % over all sensor locations
    for t = 1:size(X,2) % over time pts
        rateVal = Xtemp(locNum,t);  % probability of firing (0 to 1)
        randVal = rand;
        if randVal<rateVal
            Xspikes(locNum,t) = 1;  % neuron spiked
            Xtemp(locNum,t+1:min(t+refPerInTPts,maxT)) = 0;  % set P(firing) to 0 for refractory period
        end
    end
end

classes = unique(G);
nClasses = length(classes);

% convert to matrix of 0s and 1s to time of first spike within each wingbeat
ptsPerCondition = size(X,2)/nClasses; % assumes equal number of points in each class
ptsPerFlap = 1/ParStruct.flapFrequency * ParStruct.sampFreq;
flapsPerCondition = floor(ptsPerCondition/ptsPerFlap);
extraPtsPerCondition = rem(ptsPerCondition,ptsPerFlap);

for iClass = 1:nClasses
    % throw out first wingstroke of each condition (these can have artifacts from simulation)
    Xspikes(:,(flapsPerCondition-1)*ptsPerFlap*(iClass-1)+1:((flapsPerCondition-1)*ptsPerFlap*(iClass-1)+ptsPerFlap)) = [];
    G((flapsPerCondition-1)*ptsPerFlap*(iClass-1)+1:((flapsPerCondition-1)*ptsPerFlap*(iClass-1)+ptsPerFlap)) = [];
    % throw out extraPtsPerCondition
    Xspikes(:,(flapsPerCondition-1)*ptsPerFlap*iClass+1:((flapsPerCondition-1)*ptsPerFlap*iClass+extraPtsPerCondition)) = [];
    G((flapsPerCondition-1)*ptsPerFlap*iClass+1:((flapsPerCondition-1)*ptsPerFlap*iClass+extraPtsPerCondition)) = [];
end


XspikesReshape = reshape(Xspikes, size(X,1), ptsPerFlap, []); 
Greshape = reshape(G, 1, ptsPerFlap, []); 

X_spt = zeros(size(X,1),nClasses*(flapsPerCondition-1));
G_spt = zeros(1,nClasses*(flapsPerCondition-1));
for t = 1:size(XspikesReshape,3)  % for each wingstroke
    [~,tmp] = max(XspikesReshape(:,:,t),[],2);
    X_spt(:,t) = tmp;
    zeroIdx = sum(XspikesReshape(:,:,t),2)==0; % no spikes during this wingstroke
    X_spt(zeroIdx,t) = 0; % originally 0, also ptsPerFlap*10
    G_spt(t) = mode(Greshape(:,:,t),2);
end

