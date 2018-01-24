function [ scoreList ] = fitMultiCurveSAXS ( simInfo, ...
    curves, simNum )

% simNum = no. of grids

tic;
%m    INPUT: curves(:,idx,:) is a particular DESSA->VirusGen simulation. 
%m    FORMAT: curves(:,idx,:) <- z_data
%m            and z_data(:,1,timestep) = abs(I_system); 
%m            (see VirusGen_forClusterSAXS.m)

%m    OUTPUT: the RMSD of an averaged (over repeats) SAXS experiment.

    display('In fitMultiCurveSAXS.m')
    lengthIq = simInfo.lengthIq; % number of points in Crysol I(q) curve.
    concNum = numel(simInfo.conc);
    scoreTemp = zeros(concNum,1);
    scoreList = [];
    expCurve = simInfo.data; %m Experimental Curves
    %m FORMAT: expCurve{1} = SAXS output for all simulation times  
    %                      = [I(q1,t1)  I(q1,t2) ...
    %                         I(q2,t1)  I(q2,t2) ...
    %                         ...       ...   
    %                         I(qn,t1)  I(qn,t2) ... ]
    %          expCurve{2} = row vector of simulation times.
   

    %m numExpTimes is the number of experimental times.    
    numExpTimes = length(expCurve{2});
    
    %m Determine experimental time points list, eTimes
    eTimes = expCurve{2};
    
    
%     display(['fitMultiCurveSAXS.m: numExpTimes and size(eTimes) are: ',...
%         num2str(numExpTimes),' ',num2str(size(eTimes))]);
       
    if exist('curves','var') && ~isempty(curves) 
      
        %m  Smooth out stochastic noise with avgCurve.m
        %m  i is remnant of Lu's code. I'll keep for now.        

        
        % (Update: fitMultiCurveSAXS.m now called for
        % each grid point separately. Thus, simNum:=1)
        for s = 1 : simNum             
           
            % %for i = 1 : concNum
            i = 1; % For any code that uses the value of 'i' below.            
        
            %m Each col of av_matrix is the decay curve associated with
            %m a particular experimental time point.
%             display(['fitMultiCurveSAXS.m: size(curves) is: ',...
%                 num2str(size(curves))]);            
            [av_matrix] = avgCurveSAXS(curves,numExpTimes,lengthIq,...
                simInfo.repeat);   
            ida = size(av_matrix,2);                
            assert(size(av_matrix,2) == length(eTimes));

            % % MT: save all simulated curves if begin.m specifies to.
            if isfield(simInfo, 'SAXScurves')
                
                % NOTE: the first cell will be empty bc I do 'end+1'
                dessaTimesUsed = simInfo.dessaTimesUsed;
                simInfo.SAXScurves{end+1} = {av_matrix,dessaTimesUsed};                         
                display('genSyntheticDataSAXS (fitMultiCurve.m): add av_matrix to simInfo for this grid')
                simInfo_synthetic = simInfo;
                save(simInfo.SAXS_File,'simInfo_synthetic'); 
                display('genSyntheticDataSAXS  (fitMultiCurve.m): SYNTHETIC DATA SAVED')
            end

            %m Calc RMSD, step 1
            for w = 1:ida
                
                assert(isequal(size(expCurve{1}(:,w)),size(av_matrix(:,w))))
%                 display('fitMultiCurveSAXS.m: assert sizes equal passed.')
                diff = expCurve{1}(:,w) - av_matrix(:,w);
                assert(~isequal(diff,expCurve{1}(:,w)))
                scoreTemp(i,w) = diff'*diff;  

 
                assert(~isnan(sum(expCurve{1}(:,w))));
                assert(~isnan(sum(av_matrix(:,w))));             
                assert(~isnan(scoreTemp(i,w)));
                assert(scoreTemp(i,w) >= 0);
                
            end 
            fprintf('fitMultiCurveSAXS.m: made it past asserts.\n');
            % %end % This 'end' is for the i (concentration) loop.            
            
            %m Calc RMSD, step 2
            %scoreList(s)   
%             display('fitMultiCurveSAXS.m: ');
            scoreList(s) = (sum(sum(scoreTemp,2))/numel(expCurve{1}))^0.5; %concNum)^0.5;

        end

        
    else            
        display('Error: cannot find simulated data.');
        return        
    end   
    
    
    scoreList = scoreList'; %m must be col vector for getOffset() in
                            %m optimizer.m
    
    display('fitMultiCurveSAXS.m: now leaving fitMultiCurveSAXS.m.');
    display(['toc = ',num2str(toc)])
end