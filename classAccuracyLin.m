function [acc, classesOut] = classAccuracyLin(X, G, XTest, GTest, sensors)
% [acc, classesOut[ = classAccuracyLin(X, G, XTest, GTest, sensors)
%
% Calculate classification accuracy on test set.  First, use optimized
% sensor locations to find w, centroids, and decision
% boundary on training data. Then assess accuracy on test data.
%

nSensorLocs =  size(XTest,1);
classes = unique(GTest);
nClasses = numel(classes);
nSensors = length(sensors);

Phi = zeros(nSensors, nSensorLocs);    % construct the measurement matrix Phi
for iSensors = 1:nSensors
    Phi(iSensors, sensors(iSensors)) = 1;
end

% learn new classifier for sparsely measured data
w_sspoc= LDA_n(Phi * X, G);
XProj = w_sspoc' * (Phi * X);

% compute centroid of each class in classifier space
centroid = zeros(nClasses-1, nClasses);
s_dev = zeros(nClasses-1, nClasses);
for i = 1:nClasses
    centroid(:,i) = mean(XProj(:,G==classes(i)), 2);
    s_dev(:,i) = std(XProj(:,G==classes(i)),[],2);
end

XProjTest = w_sspoc' * (Phi * XTest);

classesOut = ones(1,size(XProjTest,2));

% use sparse sensors to classify X
if nClasses == 2
    
    sep_line = mean(centroid);
    
    if centroid(1) < sep_line
        classesOut(XProjTest>sep_line) = 2;
    else
        classesOut(XProjTest<=sep_line) = 2;
    end
    
       
else
    for i = 1:size(XProjTest,2)
        D = centroid - repmat(XProjTest(:,i), 1, size(centroid,2));
        D = sqrt(sum(D.^2, 1));
        
        [~, classesOut(i)] = min(D);
    end
end


acc = sum(classesOut==GTest)./length(GTest);