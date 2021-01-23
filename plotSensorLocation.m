function [] = plotSensorLocation(sensors, acc, Pars)
% plotSensorLocation
% plots locations of sets of sensors on rectangular wing, 
% with sensors belonging to each "set" colored differently

%find number of sensor sets
numElements = Pars.chordElements*Pars.spanElements;
maxSensor = max(sensors);
numSets = idivide(maxSensor,int16(numElements),'ceil');

f1 = figure('DefaultAxesFontSize', 18);  
hold on;

for iSet = 1:numSets
    %identify which sensors belong in this set
    ind = (sensors <= numElements*iSet) & (sensors > numElements*(iSet -1));
    [sensorsRowIdx1,sensorsColIdx1] = ind2sub([Pars.chordElements,Pars.spanElements],(sensors(ind)-double((iSet-1)*numElements)));
    h1 = plot(sensorsColIdx1,sensorsRowIdx1,'o','MarkerSize',10);
    set(h1, 'markerfacecolor', get(h1, 'color'));
end

box on
xlabel('span distance'); ylabel('chord distance')
xlim([1 Pars.spanElements]); ylim([1 Pars.chordElements])
h = title(['accuracy = ' num2str(round(acc*1000)/10) '%']);
P = get(h,'Position');
set(h,'Position',[P(1) P(2)+0.3 P(3)]);