function [ ] = inLoop2(homeFolder)
% Loads from and Saves to /home/, not /scratch/
    display('In inLoop2.m')
    load([homeFolder,'paramFile.mat']);
%     load([homeFolder,'iterPrefix.mat']);
    load([homeFolder,'xFile.mat']); 
    load([homeFolder,'loop_eval_criteria.mat'],'count','minScore','noImprove'); 
    display(['inLoop2.m: paramFile.mat,iterPrefix.mat,xFile.mat,',...
        'loop_eval_criteria.mat exist and are loaded.'])

    %load([homeFolder,'y.mat']);
    load([homeFolder,'y_iter',num2str(count),'.mat'],'y');
    display(['inLoop2.m: y_iter',num2str(count),'.mat exists and is loaded.'])
    
    

    nPara = size(optInfo.paraList,2)
    s = 3%0.05; %Search small box enclosing sphere the points were originally sampled from. %optInfo.initGridSize;
    
    y
    tempFile = [fileInfo.homeWorkFolder,fileInfo.prefix,'.mat']
%     if exist(tempFile,'file'), system(['rm -f ',tempFile]); end
    
    optInfo.score
    optInfo.offset
    optInfo.score = [optInfo.score;y(:,1)]
    optInfo.offset = [optInfo.offset;x]

    save([homeFolder,'paramFile.mat'],'fileInfo','qInfo','simInfo','optInfo',...
        'count','minScore','noImprove');

    display('InLoop2.m: Making prediction...');

    lb = ones(nPara,1) * optInfo.lb * s;
    ub = ones(nPara,1) * optInfo.ub * s;
    dx = (ub - lb) * optInfo.eps;
%     params = struct('bounds',{lb,ub},'nreq',xSize,'p',0.5);
    % MT: In a comparison with my Gaussian Process approach using 7 kernels and
    % each with 3 exploration/exploitation parameter values (21 total new
    % points), I will request that SNOBFIT return 21 new points. SNOBFIT
    % will determine the distribution of these, from more local to global.
    params = struct('bounds',{lb,ub},'nreq',21,'p',0.5);
    
    if count == 1	
        display('inLoop2.m: count==1 case.')
        try
            [newPoints,xBest,yBest] = snobfit(tempFile,x,y,params,dx);
        catch ME
        end
        if ~exist('newPoints','var')
           ME.identifier
           ME.message
           display('QUITING inLoop2.m NOW')
           quit
        end
    else
        display('inLoop2.m: count~=1 case.')
        try
            [newPoints,xBest,yBest] = snobfit(tempFile,x,y,params);
        catch ME
        end
        if ~exist('newPoints','var')
           ME.identifier
           ME.message
           display('QUITING inLoop2.m NOW')
           quit
        end

    end

    fprintf('\nNew requested offsets (i.e., grid, x):\n');
    x = newPoints(:,1:nPara)
    xSize = size(x,1)
    %MT_Jan25 cleanup(fileInfo,iterPrefix);
    display('inLoop2.m: Sample scores:\n');%m I uncommented this.
    fprintf('%e\n',[minScore;y(:,1)]);%m I uncommented this.
    fprintf('inLoop2.m: Best offset of this iteration:\n');
    xBest
    fprintf('\n inLoop2.m: Best score of this iteration: %e\n',yBest);

    
    if yBest < minScore

        [bs,cs] = setParameter(optInfo.BS,optInfo.CS,...
            optInfo.paraList,xBest);
        predictPara.BS = bs;
        predictPara.CS = cs;

        if minScore == inf 

            display('Initialization finished.');
            fprintf('Initial parameters:\n');

        else

            display('Improved!');
            fprintf('Improved! New parameters:\n');

        end

        for i = 1 : size(bs,1)

            fprintf('%s %s %.15e %.15e\n',bs{i,1},bs{i,2},...
                bs{i,3}(1),bs{i,3}(2));

        end

        for i = 1 : size(cs,1)

            fprintf('%s %s %.15e\n',cs{i,1},cs{1,2},cs{i,3});
            display('inLoop2.m: fprintf(cs stuff) done')
        end

        minScore = yBest
        noImprove = 0

    else

        display('No improvement. Decrease grid size...');

        noImprove = noImprove + 1;

    end

    %mt: displayed above.   
% % %     fprintf('\nNew requested offsets:\n');
% % % 
% % %     for i = 1 : xSize
% % % 
% % %         fprintf('%.15e ',x(i,:));
% % %         fprintf('\n');
% % % 
% % %     end

    fprintf('Iteration %d finished.\n\n',count);
    fprintf('Iteration %d finished.\n\n',count);

    save([homeFolder,'xFile.mat'],'x','xSize');
    save([homeFolder,'loop_eval_criteria.mat'],'minScore','noImprove','count')
    save([homeFolder,'storage/newPoints',num2str(count),'_incl_SNOBFIT_type.mat'],'newPoints')
    display('first save done (xFile and loop...mat).')
    % Save another copy of loop_eval to be used in PI_script.sh (source)
    system(['[ -f ',homeFolder,'loop_eval_criteria.txt ] && rm -r ',...
        homeFolder,'loop_eval_criteria.txt']);
    formatspec = '%s \n';
    fid_loop = fopen([homeFolder,'loop_eval_criteria.txt'],'w');
    fprintf(fid_loop,formatspec,['minScore=',num2str(minScore)]);
    fprintf(fid_loop,formatspec,['noImprove=',num2str(noImprove)]);
    fprintf(fid_loop,formatspec,['count=',num2str(count)]);
    fclose(fid_loop);
    
    save([homeFolder,'newTime.mat'],'predictPara');
    save([homeFolder,'paramFile.mat'],'fileInfo','qInfo','simInfo','optInfo');
    display('Post inLoop2.m')
end

