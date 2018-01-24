function [ simData,virustype ] = importSimSAXS ( simFile,virusXML )

    display('In importSimSAXS.m')
    tic;
    simData = [];

    try        
        fid = fopen(simFile,'r')        
    catch fopenErr        
        fprintf('Failed to open file: %s\n',simFile);
        throw(fopenErr);        
    end
    
    if fid < 0
        fprintf('Failed to open file: %s\n',simFile)
        return
    end
        
    % Look at the first line in the .dat file
    str = fgetl(fid);
    
    if isempty(str) || str(1) == -1        
        fprintf('Error opening file: %s\n',simFile);        
        return    
    end 
  
    fclose(fid);
    


% Indices start at zero (strange for matlab) when indexing dlmread.
  

  % -----------------------------------------------------------------------  
    line_one = dlmread(simFile,' ',[0 0 0 1]);
    numlines = line_one(1);
    divline = line_one(2); % ----------------

    if numlines == 314 && divline == 628 && divline2 == 999 && ...
            divline3 == 999 && divline4 == 999
        simData = [];
        trajectories = [];
        virustype = [];
        display('importSimSAXS: simFile still contains 314 628 999 999 999 as first line. Return.');
        
        return
    end

    sim_matrix = dlmread(simFile,' ',[1 0 divline-2 201]);
  
    vec_matrix = dlmread(simFile,' ',[divline 0 numlines-1 12]);
    
    
    C = strfind(virusXML,'CCMV');
    c = strfind(virusXML,'ccmv');
    HPV = strfind(virusXML,'HPV');
    hpv = strfind(virusXML,'hpv');
    HBV = strfind(virusXML,'HBV');
    hbv = strfind(virusXML,'hbv');
    numSubunits = sim_matrix(1,2);
%     numTimeSteps = size(traj_matrix,1);
    
%     if ~isequal(size(timeList,1),size(reactionList,1))
%         display('importSimSAXS: length of timeList, and reactionList not equal. Return.')
%         simData = 'failed'
%         trajectories = 'failed'
%         virustype = 'failed'
%         return
%     end
    
    if ~isempty(c) || ~isempty(C)
        virustype = 'ccmv';
    end
    if ~isempty(hpv) || ~isempty(HPV)
        virustype = 'hpv';
    end
    if ~isempty(hbv) || ~isempty(HBV)
        virustype = 'hbv';
    end
    

    simData = {sim_matrix,vec_matrix};    
    
    clear sim_matrix vec_matrix
    display('Leaving importSimSAXS.m')
end