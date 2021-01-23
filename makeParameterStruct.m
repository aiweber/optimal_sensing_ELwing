function [Pars] = makeParameterStruct()
% [Pars] = makeParameterStruct()
%
% Function that sets parameter values for wing_sensors_main


% parameters for Euler-Lagrange simulation
Pars.runELSim      = 0;      % 0 to skip E-L simulation if a strain set already exists for these parameters
Pars.rollRots      = [0 0];  % phi; rotation rates to simulate (rad/s)
Pars.pitchRots     = [0 0];  % theta; rotation rates to simulate (rad/s)
Pars.yawRots       = [0 10]; % psi; rotation rates to simulate (rad/s); Thomas used 0,10
Pars.phi_dist      = 0.312 *[1 1];  % roll (flapping) axis
Pars.theta_dist    = 0.1 *[0 0];    % pitch axis
Pars.psi_dist      = 0.1 *[1 1];    % yaw axis
Pars.E             = 3e7*10^(-1/4);    % Young's modulus: [3e7*10^(-5/8),3e7*10^(9/16)]; stiffness factor 1 = 3e7

Pars.chordElements = 26;
Pars.spanElements  = 51;
Pars.simStartup    = 1;  % in seconds
Pars.simEnd        = 4;  % in seconds
Pars.flapFrequency = 25;
Pars.harmonic      = 0.2;
Pars.sampFreq      = 1e4;
Pars.xInclude      = 0;  % chordwise
Pars.yInclude      = 1;  % spanwise

% parameters that control neural encoding
%    go across columns for parameters that get passed to functions as vectors;
%    go down rows for parameters to be swept over    
Pars.staWidth      = 4; 
Pars.staFreq       = 1; % a/2/pi cycles/ms
Pars.staDelay      = 5;
Pars.nldShift      = .2;
Pars.nldGrad       = 50;

Pars.normalizeVal  = 3.7732e-4; 
Pars.subSamp       = 1;

Pars.refPer        = 15; % refractory period, ms (only used if Pars.spikesFlag = 1)
Pars.spikeReps     = 10;  % how many times to generate spikes from simulated data (needed because spiking data is downsampled to once/wingstroke)

% parameters for splitting into training and test sets
Pars.trainFraction  = 0.9; % fraction of data to use for training vs test

% parameters for SSPOC
Pars.rmodes             = 20;   % how many modes to keep during dim reduction
Pars.wTrunc             = 3;    % desired number of sensors; choose an integer between 1 and 30
Pars.elasticNet         = 0.9;
Pars.penaltyVec         = 0;    % must be row vector [1, Pars.chordElements x Pars.spanElements] or 0
Pars.eps                = 1e-6; % tolerance for w reconstruction (only used if numClasses>2)
Pars.columnCoupling     = 0;    % column coupling weight to promote nonzero rows of s (only used if numClasses>2)
