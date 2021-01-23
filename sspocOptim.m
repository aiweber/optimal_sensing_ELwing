function [sensors, cutoffLim, s] = sspocOptim(w_t, Psi, SspocPrs, numClasses)
% [sensors, cutoffLim, s] = sspocOptim(w_t, Psi, SspocPrs)
%
% This function performs optimization to find optimal sensor locations.
% Solves Psi'*s = w, using elastic net penalty on s
%
% Inputs:
%   w_t: decision vector to reconstruct
%   Psi: reduced data matrix; nSensorLocs x SspocPrs.rmodes
%   SspocPrs: structure containing parameters for optimization.  
%    Fields:
%     elasticNet: relative weight of L1 and L2 norms for optimization,
%                 between 0 and 1
%     rmodes: number of modes kept during dim reduction 
%     penaltyVec: penalty on values of s; used for implementing spatial 
%                 penalties.  Must be zero or a vector of length nSensors.
%     eps: tolerance for w vector reconstruction (only used if >2 classes)
%     numClasses: if numClasses=2, Psi'*s == w_t
%                 if numClasses>2, frob(Psi'*s-w_t) < eps
%     
% Outputs:
%   sensors: sparse matrix of optimal sensor locations
%   cutoffLim: cutoff value used to determine number of sensors kept
%   s: full vector of sensor weights found in optimization step (before 
%      cutoffLim applied)
%
% This function makes use of the cvx package (http://cvxr.com/cvx/)

[n, ~] = size(Psi);  % number of sensor locations

% checks to make sure penaltyVec is of correct size
if isempty(SspocPrs.penaltyVec)
    SspocPrs.penaltyVec = 0;
end
if length(SspocPrs.penaltyVec) ~= n && numel(SspocPrs.penaltyVec)~=1
    disp('error: Incorrect size penaltyVec. Must be zero or a vector of length nSensors.')
end
if isrow(SspocPrs.penaltyVec)
    SspocPrs.penaltyVec = SspocPrs.penaltyVec';
end

if length(SspocPrs.penaltyVec) > 1  % if a vector
    SspocPrs.penaltyVec = repmat(SspocPrs.penaltyVec,1,numClasses-1);
end


if numClasses == 2
    cvx_begin quiet
    variable s( n );
    minimize(   SspocPrs.elasticNet*norm(s,1) ...
        + (1-SspocPrs.elasticNet)*norm(s,2) ...
        + norm(s.*SspocPrs.penaltyVec,2)  );
    subject to
        Psi'*s == w_t; %#ok<EQEFF>
    cvx_end
else
    v = ones(numClasses-1,1);
    cvx_begin quiet
    variable s( n, numClasses-1 );
    minimize(   SspocPrs.elasticNet*( norm(s,1) + SspocPrs.columnCoupling*norm(s*v,1)) ...
        + (1-SspocPrs.elasticNet)*norm(s,2) ...
        + norm(s.*SspocPrs.penaltyVec ,2)  );
    subject to
        norm(Psi'*s-w_t,'fro')<=SspocPrs.eps; %#ok<VUNUS>
    cvx_end
end
    
[~, I_top] = sort( sum(abs(s),2),'descend');
sensorsSort = I_top(1:SspocPrs.rmodes);

cutoffLim = norm(s, 'fro')/SspocPrs.rmodes*sqrt(numClasses-1);
sensors = sensorsSort(  abs(s(sensorsSort))>= cutoffLim );
