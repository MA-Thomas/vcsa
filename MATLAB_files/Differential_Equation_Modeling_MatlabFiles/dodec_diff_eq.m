function [ dcdt ] = dodec_diff_eq( t,conc,a,b,s,k )


O_monomer = 5;
%k_on = 10^2;

dcdt = zeros(length(conc),1);

[unique_mn_rows,ib,iu] = unique( b(:,[1,2]),'rows' );
% Thus, b( ib,[1,2] ) will return a list of all existing species.
% 'cn' below will iterate through each of these species.
assert( isequal( size( b( ib,[1,2] ),1 ), length(conc) ) )

% for cn = 1:length(conc)
for cn = 1:size( b( ib,[1,2] ),1 )   

    term1 = 0;
    term2 = 0;
    term3 = 0;
    term4 = 0;

    lenk = length(k);
    switch lenk
    case 1    
    % -----------------------------------------------------
    % Model of Misra & Schwartz (Eq.1)
    % k_on is independent of oligomer. I.e. it is a scalar
    % -----------------------------------------------------
    k_on = k;
    % We will sum over fixed (j,k), where j==m-1. (See Eq.1 Misra&Schwartz)
    % List the rows of 'b' whose [j,k] (i.e. cols [3,4]) values 
    % correspond to conc(cn) (i.e. to b(ib(cn), [1,2]) ).
    % We will then sum over the list's various (m,n) associated with 
    % the single (j,k).
    b_mat = repmat( b(ib(cn),[1,2]), size(b,1), 1 );
    list = sum( b(:,[3,4]) == b_mat, 2) == 2;

    % Term 1. "Raising [j,k] by breaking (m,n)"
    % Sum over the the (m,n) values corresponding to the current (j,k).
    term1 = sum( b(list,5).*s(list,5).*conc(iu(list)) );
   
    % Term 2. "Lowering [j,k] by forming (m,n)"
    term2 = sum( a(list,5) )*O_monomer*conc(cn)*conc(1);
  
    
    % We will sum over fixed (j,k), where j==p+1. (See Eq.1 Misra&Schwartz)
    % List the rows of 'b' who's [j,k] (i.e. cols [1,2]) values 
    % correspond to conc(cn) (i.e. to b(ib(cn), [1,2]) ).    
    b_mat = repmat( b(ib(cn),[1,2]), size(b,1), 1 );
    list = sum( b(:,[1,2]) == b_mat, 2) == 2;    
    
    % Term 3. "Lowering [j,k] by breaking (j,k) bonds"
    term3 = sum( b(list,5).*s(list,5) )*conc(cn); 
    
    % Term 4. "Raising [j,k] by bond formation"
    % How to compute [p,q]? 
    % We akready have the list of rows of 'a'. We need to know to which
    % elements of conc they correspond. Look at a(list,[1,2]).
    [C,ib_list,ia_list] = intersect(b( ib,[1,2] ),a(list,[1,2]),'rows');
    if ~isempty(C)
        term4 = sum( a(list,5).*conc(ib_list) )*O_monomer*conc(1);
    end
    total = k_on*(term1-term2) - k_on*(term3-term4);
    
    dcdt(cn) = k_on*(term1-term2) - k_on*(term3-term4);
    
    otherwise
    % ------------------------------------------------------------
    % My Model: k is oligomer dependent for the forward reactions
    % ------------------------------------------------------------ 
    k_on = 10^2; % ground truth value - to be used for backward reactions
    k_oligomer = k'; % length(k_oligomer) == length(conc)
    
    b_mat = repmat( b(ib(cn),[1,2]), size(b,1), 1 );
    list = sum( b(:,[3,4]) == b_mat, 2) == 2;

    % Term 1. "Raising [j,k] by breaking (m,n)"
    term1 = k_on*sum( b(list,5).*s(list,5).*conc(iu(list)) );
   
    % Term 2. "Lowering [j,k] by forming (m,n)"
    term2 = sum( a(list,5) )*O_monomer*conc(cn)*conc(1)*k_oligomer(cn);
  
    
    % We will sum over fixed (j,k), where j==p+1. (See Eq.1 Misra&Schwartz)
    % List the rows of 'b' who's [j,k] (i.e. cols [1,2]) values 
    % correspond to conc(cn) (i.e. to b(ib(cn), [1,2]) ).    
    b_mat = repmat( b(ib(cn),[1,2]), size(b,1), 1 );
    list = sum( b(:,[1,2]) == b_mat, 2) == 2;    
    
    % Term 3. "Lowering [j,k] by breaking (j,k) bonds"
    term3 = k_on*sum( b(list,5).*s(list,5) )*conc(cn); 
    
    % Term 4. "Raising [j,k] by bond formation"
    [C,ib_list,ia_list] = intersect(b( ib,[1,2] ),a(list,[1,2]),'rows'); 
    if ~isempty(C)
%         display('a,conc,k_olig below')
%         a(list,5)
%         conc(ib_list)
%         k_oligomer(ib_list)
        term4 = sum( a(list,5).*conc(ib_list).*k_oligomer(ib_list) )*O_monomer*conc(1);
    end
    
    total = (term1-term2) - (term3-term4);
    
    dcdt(cn) = (term1-term2) - (term3-term4);    
    end
    
end

end

