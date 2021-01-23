function results = formatResults(sensors, acc, sensors10, acc10, accAll, Pars)
% results = formatResults(sensors, acc, sensors10, acc10, accAll, Pars)
%
% Creates a single row table entry containing sensor locations, accuracy, and parameters 
% 
timestamp = datetime;
sensors = mat2cell(sensors,length(sensors),1);
sensors10 = mat2cell(sensors10,length(sensors10),1);
results = table(timestamp, sensors, acc, sensors10, acc10, accAll);
fields = fieldnames(Pars);
for i = 1:length(fields)
    Pars.(fields{i}) = num2cell(Pars.(fields{i}),[1 2]);
end
parTable = struct2table(Pars, 'AsArray', true);

results = [results parTable];
end

