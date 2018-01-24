function [ ] = inLoop1(homeFolder)
    display('In inLoop1.m')
    load([homeFolder,'paramFile.mat']);
    load([homeFolder,'loop_eval_criteria.mat'],'minScore','noImprove','count'); 
    
    
    nPara = size(optInfo.paraList,2); %Dimension of parameter space
    
    % If this is first time in inLoop1.m, prepare parameter matrix, x.
    % Otherwise, load the new points generated in inLoop2.m of prev round.
    if exist([homeFolder,'xFile.mat'],'file')
        load([homeFolder,'xFile.mat']); 
        display('inLoop1.m: xFile.mat loaded')
    else
%         s = optInfo.initGridSize;

         m = 100 % number of points (xSize)
         r = 3; %9  % radius of offset/parameter hypersphere.
         x = randsphere(m,nPara,r); % Sample interior points from nPara-dim sphere of radius r
         xSize = size(x,1);
    end    
    
    count = count + 1  
    
    jobList = [];
    dataList = [];
    iterPrefix = [fileInfo.prefix,'_iter',num2str(count)]
 
    display(xSize)
    display(x)
    display('')
    for i = 1 : xSize
        
        fprintf('\rPreparing sample points: %d out of %d',i,xSize);
        [bs,cs] = setParameter(optInfo.BS,optInfo.CS,...
            optInfo.paraList,x(i,:));
        jobPrefix = [iterPrefix,'/grid',num2str(i)];
        
        % Update [jobList,dataList] with entries for current param set.
        [jobList,dataList] = prepareJob(jobPrefix,jobList,...
            dataList,bs,cs,fileInfo,simInfo,iterPrefix,i);

    end
   
    fid_j = fopen([homeFolder,'job_list.txt'],'w')

    formatspec = '%s \n';
    for i=1:length(jobList)
        fprintf(fid_j,formatspec,jobList{i});
    end
    fclose(fid_j)
    
    jobNum = length(jobList)
    save([homeFolder,'dataList.mat'],'dataList');
    save([homeFolder,'jobNum.mat'],'jobNum');
    save([homeFolder,'iterPrefix.mat'],'iterPrefix');
    save([homeFolder,'loop_eval_criteria.mat'],'minScore','noImprove','count'); 
    save([homeFolder,'xFile.mat'],'x','xSize');
    
    % Save another copy of loop_eval to be used in PI_script.sh ('source' bash command)
    fid_loop = fopen([homeFolder,'loop_eval_criteria.txt'],'w');
    fprintf(fid_loop,formatspec,['minScore=',num2str(minScore)]);
    fprintf(fid_loop,formatspec,['noImprove=',num2str(noImprove)]);
    fprintf(fid_loop,formatspec,['count=',num2str(count)]);
    fclose(fid_loop);
    
    % Save variables at current round to storage/.
    dL = ['dataList',num2str(count),'.mat'];
    jN = ['jobNum',num2str(count),'.mat'];
    iP = ['iterPrefix',num2str(count),'.mat'];
    xF = ['xFile',num2str(count),'.mat'];    
    pF = ['paramFile',num2str(count),'.mat'];
    leC = ['loop_eval_criteria',num2str(count),'.mat']
    save([homeFolder,'storage/',dL],'dataList');
    save([homeFolder,'storage/',jN],'jobNum');
    save([homeFolder,'storage/',iP],'iterPrefix'); 
    save([homeFolder,'storage/',xF],'x','xSize');
    save([homeFolder,'storage/',pF],'fileInfo','optInfo','qInfo','simInfo');
    save([homeFolder,'storage/',leC],'count');
    
    display('Post inLoop1.m')    
    
end