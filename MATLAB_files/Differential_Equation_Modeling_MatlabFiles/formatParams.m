function k_matrix = formatParams(parameters,b)
% For each 12D parameter vector (each row) in 'parameters', 
% there will be a corresponding row in k_matrix.
% The rows of 'k_matrix' will have the same length as 'conc'.

[unique_mn_rows,ib,ic] = unique( b(:,[1,2]),'rows' );
conc = b(ib,[1,2]);

k_matrix = zeros( size(parameters,1), size(conc,1) );

% Iterate over all oligomers
for o = 1:size(conc,1)
    
    % oligomer_size takes values in {1,2,...,12}
    oligomer_size = conc(o,1);
    
    % OPTIONAL: 4 parameter search
    % Pretend there are only 4 meaningfully distinct oligomer
    % sizes.
    % assign oligomers of sizes 1,2 to the same parameter.
    % 3-6 to the same parameter; 7,10 to the same parameter; 11,12...
    % This implies there are 4 distinct parameters rather than 12.
    if size(parameters,2)==4
       if oligomer_size < 3
           oligomer_size = 1; 
       elseif oligomer_size < 7
           oligomer_size = 2;
       elseif oligomer_size < 11
           oligomer_size = 3;
       else
          oligomer_size = 4;
       end           
    end

    % OPTIONAL: 6 parameter search
    if size(parameters,2)==6
       if oligomer_size < 3
           oligomer_size = 1; 
       elseif oligomer_size < 5
           oligomer_size = 2;
       elseif oligomer_size < 7
           oligomer_size = 3;
       elseif oligomer_size < 9
           oligomer_size = 4;
       elseif oligomer_size < 11
           oligomer_size = 5;
       else
          oligomer_size = 6;
       end           
    end
    
    % OPTIONAL: 2 parameter search
    if size(parameters,2)==2
       if oligomer_size < 4
           oligomer_size = 1; 
       else
           oligomer_size = 2;
       end
           
    end
   
   % Iterate over all parameter vectors
   for p = 1:size(parameters,1)

       k_matrix(p,o) = parameters(p,oligomer_size);
       
   end
    
end


end

