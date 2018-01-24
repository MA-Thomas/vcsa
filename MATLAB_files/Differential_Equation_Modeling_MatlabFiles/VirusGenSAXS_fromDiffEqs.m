function [z_data] = VirusGenSAXS_fromDiffEqs(t,NUMmatrix,t_gt,b,...
    n_mers_cell,face_centers,crysolIq)
% Inputs: t,NUMmatrix are output of ode45, suitably altered from conc to
%         Number units.
%         t_gt are the ground truth time points. etc.
% Outputs: z_data (SAXS data)

qvec = [0:0.01:0.5]'; % Default q values crysol uses to calculate I(q). 
I_system = zeros(length(crysolIq),1);

%{
              NO LONGER NEEDED BC I FIGURED OUT HOW TO SPECIFY THE
              TIMEPOINTS RETURNED BY ode45

% % Select the appropriate time points to consider from t (i.e. t_diffeq)
% Initialize time_indices.
% For each time point t_gt(i), find the time point t(j) that is larger than
% t_gt(i) and s.t.  j is not already in 'time_indices'.
% % % % % time_indices = zeros(1,length(t_gt));
% % % % % for i = 1:length(t_gt)
% % % % %     list = t >= t_gt(i);
% % % % %     larger_times = find(list);
% % % % %     
% % % % %     if ~isempty(larger_times)
% % % % %         % update j.
% % % % %         j = larger_times(1); 
% % % % %     else
% % % % %         % don't update j. Add the same j to time_indices again.
% % % % %     end
% % % % %     
% % % % %     time_indices(i) = j;
% % % % %     
% % % % % end
% % % % % 
% % % % % t_diffeq = t(time_indices);
% % % % % NUMmatrix = NUMmatrix(time_indices,:);
%}
t_diffeq = t';


% % Construct simMatrix. First convert conc matrix to Number matrix.
% Cols of conc correspond to successive n-mers, from monomer to complete.
simMatrix = [t_diffeq, NUMmatrix];
z_data = zeros(length(qvec),1,length(t_diffeq));


% % From VirusGen_forClusterSAXS.m:
% CURRENT INTERPRETATION Dec.18 2016
% % I(q) = SUM_over_partialAssemblies( crysolIq*numSubs_i * S(q)_i ) where 
% % S(q)_i is structure factor partial assembly i. 
% % We still use the pairwise distances within each partial assembly
% % to calculate S(q)_i. 
% % There is no structure factor at the level of the entire system.
% % I.e., we take the entire system to be dilute (each PA far from others)
% % so the structure factor at the level of the entire system is 1.

% I am going to treat each pentamer as a subunit. 
% The ordering of columns in NUMmatrix is the same as the ordering 
% of b(:,[1,2]). 
% So we can use b(:,[1,2]) to choose which matrices in n_mers_cell to look
% at in order to determine the subunits present in an oligomer, and hence
% the pairwise distances present in that oligomer. 
[~,ib,iu] = unique( b(:,[1,2]),'rows' );
species = b(ib,[1,2]);

for it = 1:length(t_diffeq)
    for is = 1:size(species,1)
        
        % If there are oligomers of the current size at the current time,
        % and they are dimer or larger:
        if NUMmatrix(it,is) > 0           
            size_species = species(is,1);
            type_species = species(is,2);
            if size_species > 1
                % Dimer or larger case 
                list_of_faces_present = find(sum(n_mers_cell{size_species}{type_species},1));

                coordinates_of_faces = face_centers(list_of_faces_present,:);
                distanceMatrix = squareform(pdist(coordinates_of_faces,'euclidean'));

                % Calculate STRUCTURE FACTOR for current partial assembly.
                % FORM FACTOR for current partial assem is (N * crysolIq).
                for q = 1:length(crysolIq)
                    Sq_vector(q) = (1/size_species)* ...
                        sum( exp(-1i * qvec(q) * distanceMatrix(:) ) );          
                end
                
                % I_system(q) = SUM_over_i_partialAssemblies( crysolIq*numSubs_i * S(q)_i )
                % Add to total intensity, I_system, the contribution from
                % the current partial assembly. These contributions add up
                % unweighted (with weight 1) because we're treating the the
                % system as a dilute mixture of partial assemblies. This 
                % effectively means the sturcture factor at the level of the
                % collection of partial assemblies is 1.
                I_system = I_system + (size_species * crysolIq).*Sq_vector';            
            
            
            else
                % Monomer case
                I_system = I_system + (size_species * crysolIq);
            end
        
        end
        
        z_data(:,1,it) = abs(I_system);
    end
    
    % Reset I_system for the next snapshot
    I_system = zeros(length(crysolIq),1);
    
    
end
    
    




end

