function binary_not_close_vecs = returnBinaryNotCloseVecsARD(...
    binary_not_close_vecs,lengthscales,remaining_params,p)
% NOTE: this will add a row to <binary_not_close_vecs> for the current
% <lengthscale>.
% p is a row vector. lengthscale is a row VECTOR. binary_not_close_vecs is a matrix of row vectors. 

    p = repmat(p,size(remaining_params,1),1);
    diffMatrix = abs(remaining_params - p);
    
    lengthscales = repmat(lengthscales,size(remaining_params,1),1);
    
    % Binary matrix specifying, for each point (row), which dimensions are
    % NOT close to p, based on the characteristic lengthscale associated
    % with that dimension.
    notCloseMatrix = diffMatrix > lengthscales;
    
    % In order for a point to count as close, ALL dimensions must be close.
    % I.e. for point i to count as close, the sum across row i must be == 0
    %      for point i to count as far, the sum across row i must be > 0.
    colSumVec = sum(notCloseMatrix,2);
    binaryVec = colSumVec > 0;   
    
    if isempty(binary_not_close_vecs)
        binary_not_close_vecs = binaryVec;
    else
        binary_not_close_vecs = [binary_not_close_vecs, binaryVec];
    end
end

