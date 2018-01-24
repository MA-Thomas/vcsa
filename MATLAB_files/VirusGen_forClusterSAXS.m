function z_data = VirusGen_forClusterSAXS(matrix,mfmatrix,crysolIq_cell)

% FORMER INTERPRETATION
% % Our S(q) comes from the wiki structure factor article, eq.4.
% % Each partial assembly gets a single structure factor, S(q).
% % Assumes that within a partial assembly, all subunits are identical.
% % This assumption may approximately hold for the different subunit types
% % of CCMV. But for other viruses, we may need a better formula for S(q).

% FORMER (UPDATED) INTERPRETATION
% % I(q) = (crysolIq*numSubunits)*S(q) where S(q) is structure factor
% % for the entire system considered as a whole (not of a partial assem).
% % We still use the pairwise distances within each partial assembly
% % to calculate S(q), only there are also no the intra assebly pairwise
% % distances to consider also. Setting these distances to inf means


% CURRENT INTERPRETATION Dec.18 2016
% % I(q) = SUM_over_partialAssemblies( crysolIq*numSubs_i * S(q)_i ) where 
% % S(q)_i is structure factor partial assembly i. 
% % We still use the pairwise distances within each partial assembly
% % to calculate S(q)_i. 
% % There is no structure factor at the level of the entire system.
% % I.e., we take the entire system to be dilute (each PA far from others)
% % so the structure factor at the level of the entire system is 1.


% % TODO: not urgent. Sq below is defined to be a sum of complex exponentials.
% % Intensity, however, should be real valued (THE OUPUT IS REAL VALUED)

% crysolIq is the I(q) curve for a single cc subunit obtained w/ crysol.
% It is a column vector.
display('VirusGen_forClusterSAXS RUNNING')

crysolIq = crysolIq_cell{1};
% % % display(['VIrusGEn_forCLusterSAXS.m: size of crysolIq is: ',...
% % %     num2str(size(crysolIq))]);

% matrix = load(matrixv);
% mfmatrix = load(mfmatrixs);
% % % display(['VirusGen_forClusterSAXS.m: size of matrix is: ',...
% % %     num2str(size(matrix))]);
% % % display(['VirusGen_forClusterSAXS.m: size of mfmatrix is: ',...
% % %     num2str(size(mfmatrix))]);
numSubunits = mfmatrix(1,2);
% % % display(['VirusGen_forClusterSAXS.m: number of subunits is: ',...
% % %     num2str(numSubunits)]);




qvec = [0:0.01:0.5]'; % Default q values crysol uses to calculate I(q). 
assert(isequal(length(qvec),length(crysolIq)))

display('VirusGen_forClusterSAXS.m: made it past qvec,crysol assert');


I_system = zeros(length(crysolIq),1);


numRows = size(matrix,1);
numTimeSteps = numRows/numSubunits;
startIndex = 1; %for indexing rows of z_data
currentrow = 1;

z_data = zeros(length(qvec),1,numTimeSteps); 
timestep = 0;


% % % display('VirusGen_forClusterSAXS: about to enter loop over snapshots.')
while currentrow < numRows
    
    timestep = timestep+1;
    %display(['timestep is :',num2str(timestep)])
%     if currentrow < (numRows-10*numSubunits+1)         
%         currentrow = currentrow + numSubunits;
%         display('Skipped')
%         continue;        
%     end     
    
    currentTimeStep = matrix(currentrow:currentrow+numSubunits-1,:); 
    %m A chunk of the total matrix at a single time step in the simulation.
    
    currentTimeStep = sortrows(currentTimeStep,1);

    
    %subunitTypes = matrix(currentrow:currentrow+numSubunits-1,17);
  
    % Get unique assembly ids, C.
    [C,ia,ic] = unique(currentTimeStep(:,2));
    
   
    % Monomer partial assemblies will contribute S(q)=1 to their
    % Intensity Curves.
    partialAssemblySizes = zeros(length(C),1);
    
    % Iterate over all (partial) assemblies.
    I_partial_assemblies = zeros(length(qvec),length(C));
    
%     display(['About to iterate over the ',num2str(length(C)),...
%         ' partial assemblies at this round.'])
    for pa = 1: length(C)        
        
               
        % Get the rows of this partial assembly in currentTimeStep. 
        cTS_rows_pa = currentTimeStep(:,2)==C(pa);        
        
        % Subunit position coordinates for this partial assembly.
        cTS_pa = currentTimeStep(cTS_rows_pa,3:5);
        
        % Record partial assembly size.
        %partialAssemblySizes(pa) = size(cTS_pa,1);       
 
    
        % If partial assembly is dimer or larger, get all 
        % pairwise distances. 
        % see https://statinfer.wordpress.com/2011/11/14/efficient-matlab-i-pairwise-distances/
        % This function treats the location vectors as columns.
        % Apply transpose to ensure they are cols!!
        N = size(cTS_pa,1); % number of subunits in current partial assem.
        if N > 1
            distanceMatrix = sqrt(bsxfun(@plus,dot(cTS_pa',cTS_pa',1)',...
                dot(cTS_pa',cTS_pa',1))-2*(cTS_pa*cTS_pa'));
            
            % Calculate STRUCTURE FACTOR for current partial assembly.
            % FORM FACTOR for current partial assem is (N * crysolIq).
            for q = 1:length(crysolIq)
                Sq_vector(q) = (1/N)* sum( exp(-1i * qvec(q) * ...
                    distanceMatrix(:) ) );          
            end
            
            % I_system(q) = SUM_over_i_partialAssemblies( crysolIq*numSubs_i * S(q)_i )
            % Add to total intensity, I_system, the contribution from
            % the current partial assembly. These contributions add up
            % unweighted (with weight 1) because we're treating the the
            % system as a dilute mixture of partial assemblies. This 
            % effectively means the sturcture factor at the level of the
            % collection of partial assemblies is 1.
            I_system = I_system + (N * crysolIq).*Sq_vector';
            
            
        elseif N == 1
            I_system = I_system + (N * crysolIq);
           
        end
    end    

    
    % Add to z_data the total Intensity curve for this time step.
%     interval = startIndex:startIndex+length(qvec)-1;
   
%     z_data(interval,3) = abs(I_system); %sum(I_partial_assemblies,2);
%     z_data(interval,1) = qvec;
%     z_data(interval,2) = repmat(currentTimeStep(1,12),length(qvec),1);
    z_data(:,1,timestep) = abs(I_system);

%     startIndex = startIndex+length(qvec);
    currentrow = currentrow + numSubunits;

    
    % JULY 13 2017: I'M MAKING THIS ADDITION NOW. PROB SHOULD HAVE BEEN
    % HERE ALL ALONG. RESET I_system FOR THE NEXT SNAPSHOT.
    % FOR PARAMETER INFERENCE, GROUND TRUTH DATA WAS CREATED WITHOUT THIS
    % ADDITION, SO I'LL LEAVE THIS COMMENTED OUT. BUT WHEN USING REAL WORLD
    % (NO KNOWN GROUND TRUTH PARAMETERS) DATA, THIS SHOULD BE UNCOMMENTED.
% % %     I_system = zeros(length(crysolIq),1);
    % ----------------------------------------------------
    
end
format long

assert(isequal(numTimeSteps,timestep))
display('VirusGen_forClusterSAXS.m: Leaving now.')
end





