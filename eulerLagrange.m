function [strain, deform] = eulerLagrange(rotRoll, rotPitch, rotYaw, ph, th, ps, par)
% [strain, deform] = eulerLagrange(rotRoll, rotPitch, rotYaw, ph, th, ps, par)
%
% Function to solve ordinary differential equations for plate model of wing flapping with concurrent body rotation 
%
% Inputs
%   rotRoll = rotation rate in roll axis (x-axis) [rad/s]
%   rotPitch = rotation rate in pitch axis (y-axis) [rad/s]
%   rotYaw = rotation rate in yaw axis (z-axis) [rad/s]
%   ph = disturbance level around roll (flapping) axis [rad/s]
%   th = disturbance level around pitch axis [rad/s]
%   ps = disturbance level around yaw axis [rad/s]
%   par = structure with simulation parameters 
%
% Outputs: 
%   strain = matrix of strain at specified locations and time stamps
%   deform = matrix of deformation at specified locations and time stamps
%
% Created by Annika Eberle 
% September 20, 2013
%   Updated: 2017/07/03  (TLM)
%   Updated: 2020/06/25 (AIW)


%% Parameterize model 
%Geometry and nodes
    nodes   = 4;
    a       = 1.25;             %half chord in cm
    b       = 2.5;              %half span in cm
    xpos    = [-a a a -a];      %x position of plate nodes
    ypos    = [0 0 2*b 2*b];    %y position of plate nodes

%Material properties
    E       = par.E;            %Young's modulus (converting from kg/m/s2 to kg/cm/s2) - currently for moth (500 GPa), but for acrylic:3e9*10^-2
    nu      = 0.35;             %Poisson's ratio for the plate - currently for moth. for acrylic: 0.35 
    G       = E / (2*(1+nu));   %Shear modulus for isotropic material
    h       = 1.27e-2;            %plate height in cm -- currently for moth, but for acrylic:1.27e-2
    density = 1180 * (1e-2)^3;  %density of plate in kg/cm^3 (converting from m^3)
    
%Simulation parameters 
    dampingfactor = 63;    %multiplier for velocity-proportional damping via mass matrix 
    flapamp = 15; 
    par.freq0 = 1;
    par.freqEnd = 10;
    par.nFreq = 15;

%% Kinematics of input angles for flapping and rotation
% create symbolic variables
syms x y t

% Specify properties for sigmoidal startup
sigprop = [1;10;3];
sigd = sigprop(1);
sigc = sigprop(2);
sign = sigprop(3);
sigmoid = (sigd.*(2*pi*par.flapFrequency*t).^sign) ./ (sigc+sigd.*(2*pi*par.flapFrequency*t).^sign);

% Local movements: flapping + disturbances

% Specify local flapping function
phiFlap = deg2rad(flapamp) ...
    .*(  sin(2*pi*par.flapFrequency*t) ...
    + par.harmonic*sin(2*pi*2*par.flapFrequency*t) ) .* sigmoid;
% Generate flapping disturbance
phi_dot_disturbance = ph*whiteNoiseDisturbance(par);
phi_disturbance = int(phi_dot_disturbance);
% Generate rotating disturbance
theta_dot_disturbance = th*whiteNoiseDisturbance(par);
theta_disturbance = int(theta_dot_disturbance);
psi_dot_disturbance = ps*whiteNoiseDisturbance(par);
psi_disturbance = int(psi_dot_disturbance);

% Generate total local movements
phi = phiFlap + phi_disturbance;
theta = theta_disturbance;
psi = psi_disturbance; 

% Global movements: rotation
% Specify global rotation function
globalangle(1) = rotRoll*t.*sigmoid;  % x-axis; roll (flapping axis); phi
globalangle(2) = rotPitch*t.*sigmoid;  % y-axis; pitch; theta
globalangle(3) = rotYaw*t.*sigmoid;  % z-axis; yaw; psi

%Velocity and acceleration of the body (i.e. center base of plate)
    v0  = [0 0 0];
    dv0 = [0 0 0];

%Shape functions and their spatial derivatives 
    N = shapefunc2(a,b,nodes,xpos,ypos); %generate shape functions 
    N = [N(3,:).';N(4,:).']; %put into matrix form 
    for i = 1:6
        dxi(:,i) = [diff(N(i),x,2);diff(N(i),y,2);2*diff(diff(N(i),x),y)]; %compute second spatial derivative 
    end

%% Generate function with equations for ODEs
thisODEfile = ['PlateODE' num2str(round(rand*1000)) '.m'];
% ensure no previous version of this ODE file exists 
if exist(thisODEfile) == 2
    display('PlateODE exists; deleting now ')
end
%     delete('functions/PlateODE.m')
clear(thisODEfile)
    
[M K Ma Ia Q] = createODEfile_rotvect(a,b,E,G,nu,h,density,dampingfactor,phi,theta,psi,globalangle,N,dxi,thisODEfile);

% pause(4) %make sure file saves before solving the ODE 
% this section is problematic
iter =1; 
while exist(thisODEfile, 'file') ~= 2 && iter<5
    pause(2)
    iter = iter+1; 
end 

%% Solve ODE 
%Specify initial conditions and damping matrix
    init = zeros(2*6,1);
%Solve ODE
    disp('solving ode')
    options = odeset('RelTol',1e-5); % originally 1e-5
    teval = 0: 1/par.sampFreq : (par.simEnd-1/par.sampFreq); 
    [~,Y] = ode45(@(t,y) feval(thisODEfile(1:end-2),t,y,v0,dv0,M,K,Ma,Ia,Q),teval,init,options);
% Delete PlateODE to ensure next simulation can create it's own 
    delete(thisODEfile)
    
%% Postprocess results 

    %Specify spatial locations where the solution will be evaluated (26x51)
    disp('postprocessing')
    xeval = linspace(-a, a, par.chordElements);        % Previously -a:2*a/10:a
    yeval = linspace(0, 2*b, par.spanElements);
    [x,y]=meshgrid(xeval,yeval);
    
    % evaluate position
    for i = 1:6
        deform_N(i,:,:) = eval(N(i))';
    end
    for j = 1:length(Y(:,1))
        for i = 1:6
            deform_i(i,:,:)= deform_N(i,:,:)*squeeze(Y(j,i));
        end
        deform(j,:,:) = squeeze(sum(deform_i,1)) ;
    end
    
    %Evaluate shape functions and their derivatives for strains, and spatial derivatives of strain
    for i = 1:6
        strainxi(i,:,:) = eval(dxi(1,i))'; %for normal strain along x axis
        strainyi(i,:,:) = eval(dxi(2,i))'; %for normal strain along y axis
    end
    
    %Multiply by solution to ODE and sum over all components to solve for actual strains and displacements
    for j = 1:length(Y(:,1))
        %Multiply over all components
        for i = 1:6
            strainx1(i,:,:)= strainxi(i,:,:)*squeeze(Y(j,i));
            strainy1(i,:,:)= strainyi(i,:,:)*squeeze(Y(j,i));
        end
        %Sum over all components
        strainx(j,:,:) = squeeze(sum(strainx1,1))*-h/2;
        strainy(j,:,:) = squeeze(sum(strainy1,1))*-h/2;
    end
    
    % different strain components can be given as output, yInclude is standard
    if par.xInclude == 1 && par.yInclude == 0
        strain = [strainx(:,:)]';
    elseif par.xInclude == 0 && par.yInclude == 1
        strain = [strainy(:,:)]';
    elseif par.xInclude == 1 && par.yInclude == 1
        strain = [strainx(:,:), strainy(:,:)]';

    else
        error('Either par.xInclude or par.yInclude must be nonzero')
    end
end