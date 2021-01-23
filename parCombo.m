function ParStruct = parCombo(par)
% ParStruct = parCombo(par)
%
%create structure with all unique combinations of parameters 
%for running parameter sweeps

%input should be a structure, each field contains a parameter type
%parameters to be swept over should be in rows


ParStruct = struct([]);  %create empty structure, to fill in parameter combinations

parNames = fields(par); %names of all the parameters
nParams = length(parNames); %how many parameters
[dimRow,~] = structfun(@size, par, 'UniformOutput', false); %dimensions of each parameter list
dimRow = struct2array(dimRow); %dimRow is dimension to sweep over
dimRow(dimRow==0) = 1; %give empty fields dimension of 1
dimProduct =prod(dimRow); %number of parameter combinations

%create hypercube of combinations for each named parameter
oneBlock =  ones(dimRow);
for k = 1:nParams
    
    repVec = ones(nParams, 1)';
    repVec(k) = dimRow(k);
    
    %use parameter indices to create paramCube, in case parameter is vector
    repIndex = 1:dimRow(k);
    paramCube = bsxfun(@times, oneBlock,  reshape( repIndex, repVec)  );
    
    %assign actual parameter values
    paramVec = par.(parNames{k});
    for k2= 1:dimProduct
        if isempty(paramVec)
            ParStruct(k2).(parNames{k})=paramVec;
        else
            ParStruct(k2).(parNames{k})=rmmissing(paramVec(paramCube(k2),:));
        end
    end
    
end

end
