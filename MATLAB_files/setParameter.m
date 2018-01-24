function [ bs, cs ] = setParameter ( bsRef, csRef, paraList, x )
%{
EXAMPLE FOR 2 PARAMS (csRef set to [])
>> [ bs, cs ] = setParameter_testing ( optInfo.BS, csRef, ...
             optInfo.paraList, [.5 .5] )

bs = 

    'bst0a'    'bst1b'    [1x2 double]
    'bst0a'    'bst1d'    [1x2 double]
    'bst0b'    'bst1a'    [1x2 double]
    'bst0b'    'bst1c'    [1x2 double]
    'bst0c'    'bst0d'    [1x2 double]


cs =

     []

% AND
>> bs{:,3}

ans =

    0.0465   -2.9986


ans =

    0.0465   -2.9986


ans =

    0.0465   -2.9986


ans =

    0.0465   -2.9986


ans =

    0.0465   -2.9986

%}




    if exist('bsRef','var') && ~isempty(bsRef)
        
        %m bsTime is matrix of current parameters. 
        %m bsRef is the struct optInfo.BS 
        bsTime = cell2mat(bsRef(:,3));
        bs = bsRef;
        
    else bs = [];
        
    end
    
    if exist('csRef','var') && ~isempty(csRef)
        
        csTime = cell2mat(csRef(:,3));
        cs = csRef;
        
    else cs = [];
        
    end
    
    %m UPDATE PARAMETERS, bsTime.
    %m Params (elements of bsTime) are in 'log space', so
    %m we can ADD grid points, x(i), directly. These param
    %m values will be exponents (10^param) in 'real space' xml files.
    for i = 1 : size(paraList,2);
        
        if paraList{i}{1} == 'b' && ~isempty(bs)
            
            bsPara1 = paraList{i}{2}; %m e.g. [1 2 3 4 5] (which on rates)
            bsPara2 = paraList{i}{3}; %m e.g. [1] (which off rates)
            bsTime(bsPara1,bsPara2) = bsTime(bsPara1,bsPara2) + x(i);
            
            %mt: SO, THE ORDER OF TERMS IN OFFSET/GRID VECTOR x IS THE
            %mt: same as the ordering I specify in paraList in begin.m.
            
        end

        if paraList{i}{1} == 'c' && ~isempty(cs)
            
            csPara = paraList{i}{2};
            %csPara2 = paraList{i}{3};
            csTime(csPara) = csTime(csPara) + x(i);
            
        end

    end

    %m STORE UPDATED PARAMS IN STRUCT.
    if ~isempty(bs)

        for i = 1 : size(bsTime,1), bs{i,3} = bsTime(i,:); end
        
    end
    
    if ~isempty(cs)

        for i = 1 : size(csTime,1), cs{i,3} = csTime(i,:); end
        
    end
    
end