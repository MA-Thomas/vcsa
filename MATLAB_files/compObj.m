function [  ] = compObj(homeFolder)
% Loads from and Saves to /home/, not /scratch/
% Note: this function assumes all jobs have already been submitted and
% finished running. 
    display('In compObj.m') 
    
    % If we are at the point in the overall inference pipeline where
    % previously evaluated points are to be re-evaluated with more
    % simulations, run compObj_resampledPoints() instead of compObj()
    if exist([homeFolder,'resampled_points_info.mat'])==2
        load([homeFolder,'resampled_points_info.mat'])
        
        if run_compObj_resampledPoints == 1
            compObj_resampledPoints(homeFolder)
            return
        end
        % The value of run_compObj_resampledPoints will be updated to <0> 
        % in inLoop2_BayesianOpt_GP.m (if previously set to 1) after
        % optInfo.offset and optInfo.score are updated with the results of
        % the extra simulations.
    end

    
    
    load([homeFolder,'paramFile.mat'])
    display('paramFile.mat loaded from home')
    
    load([homeFolder,'jobNum.mat'],'jobNum');
    display('jobNum.mat loaded from home')

    load([homeFolder,'dataList.mat'],'dataList');
    display('dataList.mat loaded from home')
    % dataList.mat contains a cell array of strings.
    % These strings  are paths to dessa data on /home/ 
    % (-see prepareJob).    

    load([homeFolder,'loop_eval_criteria.mat'],'minScore','noImprove','count');
    display('loop_eval_criteria - count - loaded from home')
    % So I can use 'count' to name y.mat with the proper iter.
    
    load([homeFolder,'xFile.mat'],'x');    
    
    simNum = jobNum / simInfo.repeat / numel(simInfo.conc)

    lengthIq = simInfo.lengthIq
    expCurve = simInfo.data
    numExpTimes = length(simInfo.dessaTimesUsed)%length(expCurve{2});

    tic
    
    % m-th column of 'list' will contain the jobNums for the m-th param.
    numRepeats = simInfo.repeat
    list = reshape(1:jobNum,numRepeats,jobNum/numRepeats)
    y = zeros( size(list,2), 2 )
    
    % Iterate over the diff params
    for s = 1:size(list,2)
        display(['s is: ',num2str(s)])
        curves = zeros(lengthIq,numRepeats,numExpTimes);        
        RMSDs = zeros( numRepeats,1 );
        
        % Iterate in parallel over the repeats at the current param, s.
        % j denotes the simulation repeat within current param, s.
        % j_vec(j) denotes the current jobNum.
        j_vec = list(:,s);
        
        % THIS SHOULD BE parfor, BUT compObj DOESN'T RUN WHEN IT IS.
        for j = 1:length(j_vec)
        display(['j_vec(j) is: ',num2str(j_vec(j))])
        dataPath = strrep(dataList{j_vec(j)},fileInfo.workFolder,...
            fileInfo.dessaFolderOnHome) 
        
        if exist(dataPath,'file')
            %[simCell,divLine] = importSimSLS(dataPath) 
            [currCurve,virustype] = importSimSAXS(dataPath,...
                fileInfo.xmlRule);
            % NOTE: traj_matrix NO LONGER AN OUTPUT OF importSimSAXS.
            % 6/6/17
            % For CCMV, six cols for each timestep of traj_matrix are:
            % #monomerType0 #monomerType1 #boundAtA #boundAtB #boundAtC #boundAtD
            
            % In case something failed in the import, move on to nxt file.
            % The RMSD of this skipped file will be set to inf.            
            if isequal(currCurve,'failed') && ...
                    isequal(virustype,'failed')
                
                RMSDs(j) = inf;
                display('failed')
                continue
            end
            
            display('compObj.m: importSim called')
            display(' ')         
            
            % % --- SLS Eq Here ---------------
            % - see bottom of this file for the SLS code.
            % % --- SAXS Eq Here --------------     
            if  ~isempty(currCurve) && ...                    
                isreal(currCurve{1}) && isreal(currCurve{2})                   
                display('compObj.m: about to call VirusGen')
                
                % Used if DESSA configured to ouput ONLY at exp time points
                currCurve = ensureAdequateTimepoints(simInfo,currCurve);
                
                curves(:,j,:) = VirusGen_forClusterSAXS(currCurve{2},...
                    currCurve{1},simInfo.crysolIq);
                display('compObj.m: VirusGen_forClusterSAXS called')
                display(' ')
                
                % Append RMSD of j-th repeated SAXS experiment, i.e.
                % curves(:,idx,:), to list of RMSDs used to calc std of
                % RMSDs (noise at particular grid point).                
                RMSDs(j) = calcRMSDforNoise(curves(:,j,:),expCurve{1});

               display(['compObj.m: Latest RMSD is: ',num2str(RMSDs(j))])

                   
            else
                display('compObj.m: BAD currCurve. QUIT')
                quit;    
            end
        end   
        end
        % Call fitMultiCurveSAXS for each grid point, not all grids.
        % This means we don't have to hold data from all grid pnts
        % in memory at once. (7-18-2016)

        % <score> is the RMSD of an averaged (over repeats) SAXS experiment.
        % calcRMSDforNoise() returns the RMSD of a single repeat.
        [score] = fitMultiCurveSAXS(simInfo,curves,1); %SAXS%
        Rs = RMSDs;
        y(s,:) = [score,...
            sqrt(mean( ((Rs(Rs<inf)) - score).*((Rs(Rs<inf)) - score) ))];
        y(s,:)

        % Save must be done by a function when using parfor.
        parsave_y([homeFolder,'y_iter',num2str(count),'.mat'],'y',y);

        % % MT: save all simulated curves.
        SAXScurve_file = [homeFolder,'SAXSfile_iter_',num2str(count),...
            '_grid_',num2str(s),'.mat'];
        parsave_curves(SAXScurve_file,'curves',curves);

        % % MT: save trajectory and RMSD data.
        Lmat_name = [homeFolder,'Lfile_iter_',num2str(count),...
            '_grid_',num2str(s),'.mat']
        parsave_RMSDs(Lmat_name,'RMSDs',RMSDs)%,'traj_cell')         

        
    end   
    
    %%}
    
    % Calculate Covariance (Kernel) Matrix. 
    % Return Kernel as optInfo.Kernel
    % ------- temp block begin --------
    optInfo.Lmat_names = []
    virustype = 'ccmv'
    for c = 1:count
        
        optInfo.Lmat_names{c} = [];
        numf = length(dir([homeFolder,'Lfile_iter_',num2str(count),...
                '_grid_*.mat']));
        for ss = 1:numf
            Lmat_name = [homeFolder,'Lfile_iter_',num2str(count),...
                '_grid_',num2str(ss),'.mat'];
            optInfo.Lmat_names{c}{end+1} = Lmat_name;
        end
    
    end
    % ------- temp block end --------  
        
    yF = ['y_iter',num2str(count),'.mat']; 
    save([homeFolder,'storage/',yF],'y');
    
    display(' ')
    display('compObj.m: All jobs finished.\n All Objectives Saved...');
    
    % Reset # repeats.
%     simInfo.repeat = 300;
    save([homeFolder,'paramFile.mat'],'fileInfo','qInfo','simInfo','optInfo');
    
    display('compObj.m: paramFile.mat saved to homeFolder.')
    display(' ')   
    
    toc
    %}
    display('Post compObj.m')
end


function r = calcRMSDforNoise(simSAXSexp,realSAXSexp)
% input: simSAXSexp is curves(:,idx,:) - each col vector is Iq curve.
% realSAXSexp is expCurve{1}(:,:) - eac col vector is Iq curve.
assert(isequal( size(simSAXSexp,3),size(realSAXSexp,2) ));
SE = zeros(1,size(simSAXSexp,3));

for u = 1:size(simSAXSexp,3)
    diff = simSAXSexp(:,1,u) - realSAXSexp(:,u);
    %SE(u) = diff'*diff/length(diff); %MT: This division is part of the "mean" in RMSD. IT MUST COME AFTER both sums which comprise the mean, not btw them.
    SE(u) = diff'*diff;
end

r = (sum(SE)/numel(realSAXSexp))^0.5;
end

function parsave_y( string1, string2, y )
save(string1,string2)
end

function parsave_curves( string1, string2, curves )
save(string1,string2)
end

function parsave_RMSDs( string1, string2, RMSDs )
save(string1,string2)
end

function compObj_resampledPoints(homeFolder)
    load([homeFolder,'resampled_points_info.mat'])
    display('In compObj.m (compObj_resampledPoints)')
    display('resampled_points are: ',resampled_points)
    % NOTE: I will manually update xFile.mat when it is time to resample
    % points. The pipeline will then be started from inLoop1.m
    % I will also ensure that resampled_points_info.mat is in the
    % same folder - for use in compObj.m and inLoop2_BayesianOpt_GP.m.
    % 'resampled_points_info.mat' IS CREATED IN
    % 'decide_points_to_resample.m'.
    % NOTE: Contents of resampled_points_info.mat:
    %       run_compObj_resampledPoints ~ {1,0}
    %       resampled_points ~ [ round, point#_in_round ]
    %       close_points_to_remove ~ [ round, point#_in_round ]
    
    % The order (rows) of resampled_points corresponds to the order (rows)
    % in xFile.mat.  That is, x(1,:) is the same point that was previously
    % evaluated as the <resampled_points(1,2)> point in round
    % <resampled_points(1,1)>.
    % <close_points_to_remove> points will be removed from optInfo.offset
    % and optInfo.score in inLoop2_BayesianOpt_GP.m
    % run_compObj_resampledPoints==1 tells compObj.m above to use
    % compObj_resampledPoints() here.
    
    
    load([homeFolder,'paramFile.mat'])
    display('paramFile.mat loaded from home')
    
    load([homeFolder,'jobNum.mat'],'jobNum');
    display('jobNum.mat loaded from home')

    load([homeFolder,'dataList.mat'],'dataList');
    display('dataList.mat loaded from home')
    % dataList.mat contains a cell array of strings.
    % These strings  are paths to dessa data on /home/ 
    % (-see prepareJob).    

    load([homeFolder,'loop_eval_criteria.mat'],'minScore','noImprove','count');
    display('loop_eval_criteria - count - loaded from home')
    % So I can use 'count' to name y.mat with the proper iter.
    
    load([homeFolder,'xFile.mat'],'x');    
    
    simNum = jobNum / simInfo.repeat / numel(simInfo.conc)

    lengthIq = simInfo.lengthIq
    expCurve = simInfo.data
    numExpTimes = length(simInfo.dessaTimesUsed)%length(expCurve{2});

    tic
    
    % m-th column of 'list' will contain the jobNums for the m-th param.
    numRepeats = simInfo.repeat
    list = reshape(1:jobNum,numRepeats,jobNum/numRepeats)
    y = zeros( size(list,2), 2 )
    
    % Iterate over the diff params
    for s = 1:size(list,2)
        display(['s is: ',num2str(s)])
        curves_new = zeros(lengthIq,numRepeats,numExpTimes);        
        RMSDs_new = zeros( numRepeats,1 );
        
        % Iterate in parallel over the repeats at the current param, s.
        % j denotes the simulation repeat within current param, s.
        % j_vec(j) denotes the current jobNum.
        j_vec = list(:,s);
        
        % THIS SHOULD BE parfor, BUT compObj DOESN'T RUN WHEN IT IS.
        for j = 1:length(j_vec)
        display(['j_vec(j) is: ',num2str(j_vec(j))])
        dataPath = strrep(dataList{j_vec(j)},fileInfo.workFolder,...
            fileInfo.dessaFolderOnHome) 
        
        if exist(dataPath,'file')
            [currCurve,virustype] = importSimSAXS(dataPath,...
                fileInfo.xmlRule);

            if isequal(currCurve,'failed') && ...
                    isequal(virustype,'failed')
                
                RMSDs_new(j) = inf;
                display('failed')
                continue
            end
            
            display('compObj.m: importSim called')
            display(' ')         
            
            % % --- SLS Eq Here ---------------
            % - see bottom of this file for the SLS code.
            % % --- SAXS Eq Here --------------     
            if  ~isempty(currCurve) && ...                    
                isreal(currCurve{1}) && isreal(currCurve{2})                   
                display('compObj.m: about to call VirusGen')
                
                % Used if DESSA configured to ouput ONLY at exp time points
                currCurve = ensureAdequateTimepoints(simInfo,currCurve);
                
                curves_new(:,j,:) = VirusGen_forClusterSAXS(currCurve{2},...
                    currCurve{1},simInfo.crysolIq);
                display('compObj.m: VirusGen_forClusterSAXS called')
                display(' ')
                
                % Append RMSD of j-th repeated SAXS experiment, i.e.
                % curves(:,idx,:), to list of RMSDs used to calc std of
                % RMSDs (noise at particular grid point).                
                RMSDs_new(j) = calcRMSDforNoise(curves_new(:,j,:),expCurve{1});

               display(['compObj.m: Latest RMSD is: ',num2str(RMSDs_new(j))])

                   
            else
                display('compObj.m: BAD currCurve. QUIT')
                quit;    
            end
        end   
        end
        
        % -----------------------------------------------------------------
        % -----------------------------------------------------------------        
        % Now load the SAXS curves (at the current point) evaluated in a
        % previous round.
        load([homeFolder,'SAXSfile_iter_',num2str(resampled_points(s,1)),...
            '_grid_',num2str(resampled_points(s,2)),'.mat']);
        curves_old = curves;
        clear curves        
        curves = [curves_old, curves_new];
        numRepeats_old = size(curves_old,2);
        numRepeats_total = numRepeats+numRepeats_old;
        
        % Now load the RMSDs (at the current point) evaluated in a previous 
        % round.
        load([homeFolder,'Lfile_iter_',num2str(resampled_points(s,1)),...
            '_grid_',num2str(resampled_points(s,2)),'.mat']);
        RMSDs_old = RMSDs;
        clear RMSDs
        RMSDs = [RMSDs_old;RMSDs_new];
        % -----------------------------------------------------------------
        % -----------------------------------------------------------------
        
        % Call fitMultiCurveSAXS for each grid point, not all grids.
        % This means we don't have to hold data from all grid pnts
        % in memory at once. (7-18-2016)

        % <score> is the RMSD of an averaged (over repeats) SAXS experiment.
        % calcRMSDforNoise() returns the RMSD of a single repeat.
        simInfo_resampled = simInfo;
        simInfo_resampled.repeat = numRepeats_total;
        [score] = fitMultiCurveSAXS(simInfo_resampled,curves,1); %SAXS%
        Rs = RMSDs;
        y(s,:) = [score,...
            sqrt(mean( ((Rs(Rs<inf)) - score).*((Rs(Rs<inf)) - score) ))];
        y(s,:)

        % Save must be done by a function when using parfor.
        parsave_y([homeFolder,'y_iter',num2str(count),'.mat'],'y',y);

        % % MT: save all simulated curves.
        SAXScurve_file = [homeFolder,'SAXSfile_iter_',num2str(count),...
            '_grid_',num2str(s),'.mat'];
        parsave_curves(SAXScurve_file,'curves',curves);

        % % MT: save trajectory and RMSD data.
        Lmat_name = [homeFolder,'Lfile_iter_',num2str(count),...
            '_grid_',num2str(s),'.mat']
        parsave_RMSDs(Lmat_name,'RMSDs',RMSDs)%,'traj_cell')         

        
    end   
    
    %%}
    

    % ------- temp block begin --------
    optInfo.Lmat_names = []
    virustype = 'ccmv'
    for c = 1:count
        
        optInfo.Lmat_names{c} = [];
        numf = length(dir([homeFolder,'Lfile_iter_',num2str(count),...
                '_grid_*.mat']));
        for ss = 1:numf
            Lmat_name = [homeFolder,'Lfile_iter_',num2str(count),...
                '_grid_',num2str(ss),'.mat'];
            optInfo.Lmat_names{c}{end+1} = Lmat_name;
        end
    
    end
    % ------- temp block end --------
   
        
    yF = ['y_iter',num2str(count),'.mat']; 
    save([homeFolder,'storage/',yF],'y');
    
    display(' ')
    display('compObj.m: All jobs finished.\n All Objectives Saved...');
    
    save([homeFolder,'paramFile.mat'],'fileInfo','qInfo','simInfo','optInfo');
    
    display('compObj.m: paramFile.mat saved to homeFolder.')
    display(' ')   
    
    toc
    %}
    display('Post compObj.m')

end

