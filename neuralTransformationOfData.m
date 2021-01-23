function [X, G] = neuralTransformationOfData(StrainSet, Pars)
% [X, G] = neuralTransformationOfData(StrainSet,Pars)
%
% Takes in strain data and converts it to neurally encoded strain via a
% linear-nonlinear model.
%
% Inputs:
%   Pars: structure with parameters for simulation.
%     Fields:
%     rots: vector of rotation rates simulated (rad/s)
%     staWidth: width of STA (ms)
%     staFreq: frequency of STA (Hz)
%     staDelay: time delay before oscillations of STA begin (ms)
%     nldShift: horizontal offset (threshold shift) of NLD
%     nldGrad: slope parameter of NLD
%     normalizeVal
%     subSamp
%
% Outputs:
%   X: data matrix; nSensorLocs x nDataPts
%   G: vector indicating classes of columns in X; 1 x nDataPts
%   SspocPrs: structure containing parameters for optimization.  Fields:
%     rmodes: number of modes to keep during dim reduction
%     wTrunc: desired number of sensors

fieldsStrainSet = fields(StrainSet);
nSensorLocs = size(StrainSet.(fieldsStrainSet{1}),1);

staFunc = @(t) cos(  Pars.staFreq*(t+Pars.staDelay)) .* exp(-(t+Pars.staDelay).^2 / Pars.staWidth.^2  );
nldFunc = @(s) (  1 ./ ( 1+exp(-Pars.nldGrad.*(s-Pars.nldShift)) ) - 0.5  ) + 0.5;

staT = -19:1000/Pars.sampFreq:0;

f = staFunc(staT);
f = f-mean(f);
if Pars.staFreq < .1
    f = ones(size(f));
end
k = sqrt(1/sum(f.^2)); 
staFilt = fliplr(k*f /0.2003 *1000/Pars.sampFreq);  


n_conv = ( Pars.simStartup  *Pars.sampFreq*Pars.subSamp +2 -length(staT) )...  % set time points of data that will be used for convolution
    : Pars.simEnd*Pars.sampFreq*Pars.subSamp;
n_out = round((Pars.simEnd-Pars.simStartup) * Pars.sampFreq*Pars.subSamp);

convMat = zeros(nSensorLocs,n_out*length(Pars.pitchRots));
G = zeros(1,n_out*length(Pars.pitchRots));
signTag = 'NP'; % negative or positive 

for iRot = 1:length(Pars.pitchRots)      % convolve data with neural filter
    thisRollRot = Pars.rollRots(iRot);
    thisPitchRot = Pars.pitchRots(iRot);
    thisYawRot = Pars.yawRots(iRot);
    
    rotSignIdxs = round((sign([thisRollRot thisPitchRot thisYawRot])+3)/2);
    rotString = ['roll' signTag(rotSignIdxs(1)) num2str(abs(thisRollRot)) '_pitch' signTag(rotSignIdxs(2)) num2str(abs(thisPitchRot)) '_yaw' signTag(rotSignIdxs(3)) num2str(abs(thisYawRot))] ;
    
    sSet = StrainSet.(['strain_' rotString]);
    strainConv = zeros(nSensorLocs,n_out);
    
    for locNum = 1:nSensorLocs
        strainConv(locNum,:) = conv(sSet(locNum,n_conv),staFilt,'valid');
    end
    convMat(:,(iRot-1)*n_out+1:(n_out*iRot)) = strainConv;
    G(:,(iRot-1)*n_out+1:(n_out*iRot)) = iRot;
end %end rotation loop

X = nldFunc( convMat/Pars.normalizeVal/Pars.subSamp );  % run filtered data through nonlinear decision function

end

