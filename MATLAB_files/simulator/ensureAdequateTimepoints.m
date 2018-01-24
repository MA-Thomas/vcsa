function [ e ] = ensureAdequateTimepoints(entries)
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

for i = 1:length(entries)

    
    
end


end

