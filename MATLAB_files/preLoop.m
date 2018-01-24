function [ ] = preLoop(homeFolder)
% Input: homeFolder is a string,e.g. '/home/marcust/HPV_v6_8param_snobfit/'
    display('In preLoop.m')    
    load([homeFolder,'paramFile.mat']);
    display(homeFolder)
    
    if ~isfield(optInfo,'eps'), optInfo.eps = 1e-3; end
    
    if ~isfield(optInfo,'BS'), optInfo.BS = []; end

    if ~isfield(optInfo,'CS'), optInfo.CS = []; end
    
    if ~isfield(optInfo,'initScore'), optInfo.initScore = inf; end
    
    if ~isfield(optInfo,'offset'), optInfo.offset = []; end
    
    if ~isfield(optInfo,'score'), optInfo.score = []; end
    
    if isfield(fileInfo,'log') && ~isempty(fileInfo.log)
        
        logFile = fopen(fileInfo.log,'a');
        
        if logFile < 0, logFile = fopen('/dev/null'); end
        
    else
        
        logFile = fopen('/dev/null');
    
    end

    fprintf('Begin search process for %s...\n\n',fileInfo.prefix);
    fprintf(logFile,'Begin search process for %s...\n\n',fileInfo.prefix);

    fprintf('Initial score: %e\n\n',optInfo.initScore);
    fprintf(logFile,'Initial score: %e\n\n',optInfo.initScore);

    predictPara.BS = optInfo.BS;
    predictPara.CS = optInfo.CS;
    minScore = optInfo.initScore;
    count = 0;
    noImprove = 0;
 
    save([homeFolder,'loop_eval_criteria.mat'],'minScore','noImprove','count')
    save([homeFolder,'paramFile.mat']);
   
    display('Post preLoop.m')
end