function [ currCurve ] = ensureAdequateTimepoints(simInfo,currCurve)
% CCMV VERSION
% This function is called from compObj.m in the case where 
% DESSA outputs data only at the experimental time points given.
% However, in the event the simulation finishes, there may be no
% DESSA output at 1 or more late experimental time points.

% This file operates under the assumption that if a DESSA
% trajectory finishes at say, 9s, yet we have experimental data
% at 12s, then the state of the simulated system would not change
% between 9s and 12s, so we can use the output at 9s as a surrogate
% for 12s.
% Similarly, the last time point for which DESSA outputs data can
% be used as an unchanging surrogate for all remaining (experimental)
% time points at which we need data.

sim = currCurve{1};
vec = currCurve{2};

% display(['The size of sim is: ',num2str(size(sim,1))])
% display(['The size of vec is: ',num2str(size(vec,1))])

% These are the times used by genSyntheticDataSAXS2.m
eTimes = simInfo.dessaTimesUsed;%[.01:.2:16,16.4:.6:20]; 

% These are the times used in the actual synthetic data stored
% on lane cluster.
% % expCurve = simInfo.data; 
% % eTimes = expCurve{2};

numDessaTimes = size(sim,1);
% length(eTimes)

if numDessaTimes < length(eTimes)
    
    display('ensureAdequateTimepoints.m: if case reached.')
    
    % Duplicate last DESSA time point so we have enough of them.
    numDupRows = length(eTimes) - numDessaTimes;
    lastRow = sim(end,:);
    toAdd = repmat(lastRow,numDupRows,1);
    sim = [sim;toAdd];
%     size(sim)
    
    % For 90 subunits, last chunk rows are end-89:end
    % For 450 subunits (5 capsids), last chunk rows are end-449:end
%     lastChunk = vec(end-89:end,:);
    lastChunk = vec(end-449:end,:);
    toAddVec = repmat(lastChunk,numDupRows,1);
    vec = [vec;toAddVec];
    
%     size(lastChunk)
%     size(vec)
    
else
    display('ensureAdequateTimepoints.m, else case, assert next.')
    assert( numDessaTimes == length(eTimes) )
    display('ensureAdequateTimepoints.m: else case, assert passed.')
end

currCurve = {sim,vec};
display('Leaving ensureAdequateTimepoints.m')
end

