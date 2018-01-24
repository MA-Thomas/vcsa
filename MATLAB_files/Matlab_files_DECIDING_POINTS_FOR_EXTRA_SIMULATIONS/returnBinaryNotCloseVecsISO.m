function binary_not_close_vecs = returnBinaryNotCloseVecsISO(...
    binary_not_close_vecs,lengthscale,remaining_params,p)
% NOTE: this will add a row to <binary_not_close_vecs> for the current
% <lengthscale>.
% p is a row vector. lengthscale is a SCALAR. 
% Each col in binary_not_close_vecs is a list specifying whether each point
% (i.e. row) of <remaining_params> is far from p, as determined by the
% current lengthscale. 

    % Examine distance btw <p> and each point in <remaining_params>.
    p = repmat(p,size(remaining_params,1),1);
    colOfDistances = sqrt(sum((remaining_params - p) .^ 2,2))   
    
    % binaryVec tells whether each point in <remaining_params> is 
    % far from <p>.
    binaryVec = colOfDistances > lengthscale;

    if isempty(binary_not_close_vecs)
        binary_not_close_vecs = binaryVec;
    else
        binary_not_close_vecs = [binary_not_close_vecs,binaryVec];
    end
end
