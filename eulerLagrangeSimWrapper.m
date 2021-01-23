function StrainSet = eulerLagrangeSimWrapper(Pars,dataDir)
% StrainSet = eulerLagrangeSimWrapper(ElPrs,dataDir)
%
% This function either loads or creates a structure of strain data.
%
% Inputs:
%   Pars: structure with parameters for Euler-Lagrange simulation
%   dataDir: directory to either look for already simulated data or to save
%     new data to
%
% Outputs:
%   StrainSet: structure of strain values for each parameter
%     combination
%
% This function relies on the eulerLagrange function to run new simulations.

StrainSet = struct();
signTag = 'NP';  % positive or negative
includeVec = [Pars.xInclude Pars.yInclude];
includeDirs = 'xy';
for iRot = 1:length(Pars.pitchRots)
    thisRollRot = Pars.rollRots(iRot);
    thisPitchRot = Pars.pitchRots(iRot);
    thisYawRot = Pars.yawRots(iRot);
    Efact = Pars.E/(3e7);
    Estring = ['_E' num2str(floor(Efact)) '-' num2str(floor(Efact*10)-floor(Efact)*10) num2str(floor(Efact*100)-floor(Efact*10)*10) num2str(floor(Efact*1000)-floor(Efact*100)*10)];
    
    rotSignIdxs = round((sign([thisRollRot thisPitchRot thisYawRot])+3)/2);
    rotString = ['roll' signTag(rotSignIdxs(1)) num2str(abs(thisRollRot)) '_pitch' signTag(rotSignIdxs(2)) num2str(abs(thisPitchRot)) '_yaw' signTag(rotSignIdxs(3)) num2str(abs(thisYawRot))] ;
    fieldName = ['strain_' rotString];
    StrainSet.(fieldName) = [];
    
    for includeIdx = 1:2
        Pars.xInclude = 0;
        Pars.yInclude = 0;
        if includeVec(includeIdx)
            
            strainDir = includeDirs(includeIdx);
            Pars.([strainDir 'Include']) = 1;
            fName = ['simStrain_' rotString '_' strainDir '_phi' num2str(round(Pars.phi_dist(iRot)*1000)) '_theta' num2str(round(Pars.theta_dist(iRot)*1000)) '_psi' num2str(round(Pars.psi_dist(iRot)*1000)) Estring '_sf' num2str(round(Pars.sampFreq)) '.mat'];
                        
            if (Pars.runELSim == 0) && exist([dataDir '/simulated_wing_data/' fName],'file')    % load data if it has already been run and saved
                loadedWorkspace = load([dataDir '/simulated_wing_data/' fName]);
                if ~isfield(loadedWorkspace,'simStrain')
                    disp(['Error: simStrain not saved for ' fName]);
                    return
                end
                if size(loadedWorkspace.simStrain,2)<Pars.sampFreq*Pars.simEnd
                    disp(['Error: Not enough time points in previously saved dataset.  Delete file ' fName ' and regenerate with additional time points.'])
                    return
                end

                StrainSet.(fieldName) = [StrainSet.(fieldName); loadedWorkspace.simStrain(:,1:Pars.sampFreq*Pars.simEnd)];
                clear loadedWorkspace
            else       % otherwise generate and save data
                
                if Pars.sampFreq>1e3  % if desired sampling frequency is greater than 1e3, simulate at 1e3Hz and then interpolate
                    desiredSampFreq = Pars.sampFreq;
                    Pars.sampFreq = 1e3;
                    desiredSimEnd = Pars.simEnd;
                    Pars.simEnd = desiredSimEnd+1/Pars.sampFreq;
                    disp(['simulating new data for rotation ' rotString])
                    strainDataTemp.simStrain = eulerLagrange(thisRollRot, thisPitchRot, thisYawRot, Pars.phi_dist(iRot), Pars.theta_dist(iRot), Pars.psi_dist(iRot), Pars);
                    % interpolate simulated data
                    tSim = 0: 1/Pars.sampFreq : (Pars.simEnd-1/Pars.sampFreq);
                    tDesired = 0: 1/desiredSampFreq : (desiredSimEnd-1/desiredSampFreq);
                    simStrain = interp1(tSim,strainDataTemp.simStrain',tDesired,'spline')';
                    Pars.sampFreq = desiredSampFreq; % set sampFreq in parameter structure back to original (desired) sampFreq
                    Pars.simEnd = desiredSimEnd; % set simEnd in parameter structure back to original (desired) simEnd
                else
                    disp(['simulating new data for rotation ' rotString])
                    simStrain = eulerLagrange(thisRollRot, thisPitchRot, thisYawRot, Pars.phi_dist(iRot), Pars.theta_dist(iRot), Pars.psi_dist(iRot), Pars);
                end
                
                StrainSet.(fieldName) = [StrainSet.(fieldName); simStrain];
                save([dataDir '/simulated_wing_data/' fName],'simStrain','Pars');
            end
        end
    end
    
end


