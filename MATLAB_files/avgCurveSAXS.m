function [ avgValuesMatrix ] = avgCurveSAXS (curves,...
    numExpTimes, lengthIq, numSimRuns )
% % Input: curveList
%m where curveList contains the simulated curves from different 
%m DESSA runs, all for the same parameter set (same grid, s).
%m curveList(:,i,j) = col Vector [I(q)]
%m i is simulation run, j is assembly time point.

%m Output: avgValuesMatrix
%m Each col of avgValuesMatrix is the I(q) curve associated with
%m a particular experimental time point. 
%m std_func_of_x - Matrix. 
%m Each row gives, for a particular #repeats used to generate rep SAXS experiments,
%m the std of the SAXS curves at successive time points.

%m ASSUMPTIONS: 
%m 1. each of the experimental I(q) curves has measurements 
%m taken at the same number of q points. The Crysol default is 50. 
 
display('In avyCurveSAXS.m.')

%numSimRuns = numel(curves);

%m curveList is cell with dimensions [numSimRuns by numExpTimes]
%m Each col of avgValuesMatrix is the decay curve associated with
%m a particular experimental time point. 
avgValuesMatrix = zeros(lengthIq,numExpTimes);


for j = 1:numExpTimes
 
    % Average the curves over every simulation run for exp time j.
    temp = sum(curves(:,:,j),2); % sum over each simulation run at time j.
    avgValuesMatrix(:,j) = temp/numSimRuns;
    
end

end

%{
        % Obtain vector of standard deviations (one at each exp time)
        std_C = std(representativeCurves,0,3);
        
        % Sum the standard deviations across exp times for this chunk.
        std_C = sum(std_C);
%}