function [ jobList, dataList ] = prepareJob ( jobPrefix, jobList, dataList,...
    bindSite, confSwitch, fileInfo, simInfo, iter_jobPrefix,i_jobPrefix)
% NOTE: confSwitch not used anymore (confSwitch may be [])
  
    display('In prepareJob.m')
    if isfield(simInfo,'maximer'), maximer = simInfo.maximer;
        
    else maximer = 200;
        
    end
     
    interpreterPath = '#!/bin/bash';
    xmlTemplate = fileInfo.xmlRule;
    concNum = length(simInfo.conc);
    jobNum = simInfo.repeat*concNum;
    concFormat = ['%0',num2str(ceil(log(concNum+1)/log(10))),'d'];
    jobFormat = ['%0',num2str(ceil(log(simInfo.repeat+1)/log(10))),'d'];
    currJob = cell(1,jobNum);
    currData = cell(1,jobNum);

    xmlFile = cell(1,concNum);
    
    bsConc = bindSite;
    bsNum = size(bindSite,1);    

    
    %MT: modifyParameter will create the necessary file structure on
    %    homeWorkFolder so that fopen functions correctly below.
    xmlFile{i} = [fileInfo.homeWorkFolder,jobPrefix,'/',num2str(i,...
        concFormat),'.xml'];
    modifyParameter(xmlTemplate,xmlFile{i},bsConc,confSwitch);         

    for j = 1 : simInfo.repeat

         currPrefix = [jobPrefix,'/',num2str(j,jobFormat),'_conc',...
             num2str(i,concFormat)];
         currJob{(i-1)*simInfo.repeat+j} = [fileInfo.homeWorkFolder,...
             currPrefix,'.sh'];
         currData{(i-1)*simInfo.repeat+j} = [fileInfo.workFolder,...
             currPrefix,'.dat'];
         currDataOnHome{(i-1)*simInfo.repeat+j} = [fileInfo.dessaFolderOnHome,...
            strrep(jobPrefix,'/','_'),'_',num2str(j,jobFormat),...
            '_conc',num2str(i,concFormat),'.dat'];              
          

        if isfield(fileInfo,'dessaFolderOnHome')

        %m Print shell script commands, as strings, into currJob file.                 

             fid = fopen(currJob{(i-1)*simInfo.repeat+j},'w');                 

             fprintf(fid,'%s \n',interpreterPath);%MT Jan.28 2016 addition

             fprintf(fid,'%s \n','[ -d /scratch/marcust/ ] || mkdir /scratch/marcust/ ');
             fprintf(fid,'%s \n',['[ -d ',fileInfo.workFolder,' ] || mkdir ',fileInfo.workFolder]); 

             fprintf(fid,'%s \n',['[ -d ',fileInfo.workFolder,...
                 iter_jobPrefix,' ] || mkdir ',fileInfo.workFolder,...
                 iter_jobPrefix]); 

             fprintf(fid,'%s \n',['[ -d ',fileInfo.workFolder,...
                 iter_jobPrefix,'/grid',num2str(i_jobPrefix),...
                 ' ] || mkdir ',fileInfo.workFolder,iter_jobPrefix,...
                 '/grid',num2str(i_jobPrefix)]); 

             fprintf(fid,'%s %s %d %f > %s\n',fileInfo.command,...
                 xmlFile{i},simInfo.interval,simInfo.time,...
                 currData{(i-1)*simInfo.repeat+j});

             %MT: If Dessa returns exit status (-1), as is the case 
             %    when addPartner() of Subunit.java fails, rerun the
             %    simulation. Dessa will choose a new random seed which
             %    should solve the problem.
             %            Try new simulation a max of 3 times.
             %            (-1) mod 256 is 255, check for this exitcode.
             fprintf(fid, 'exitStatus=$? \n');  

             fprintf(fid, 'b=0 \n');
             fprintf(fid, 'while [[ $b -lt 4 && $exitStatus == 255 ]]; do\n' );
             fprintf(fid, '> %s\n ',currData{(i-1)*simInfo.repeat+j});
             fprintf(fid, 'set -o pipefail \n');
             fprintf(fid, '%s %s %d %f > %s\n',fileInfo.command,...
                 xmlFile{i},simInfo.interval,simInfo.time,...
                 currData{(i-1)*simInfo.repeat+j});                 
             fprintf(fid, 'exitStatus=$? \n');   
             fprintf(fid, 'let b=b+1; \n');
             fprintf(fid, 'done \n'); 
             % Note: "> filename" 
             % The > truncates file "filename" to zero length.
             % If file not present, creates zero-length file (same effect as 'touch').                  



            %%%% -New cluster setup. currData{..} now on /scratch.
            %%%% -dataTemp not necessary.

            % See shell_testing.m to test this section.
            % Writing Unix Bash Shell script to handle DESSA output.
            % 'divline' is the ln of dataTemp which 
            % separates sim from vec DESSA output.
            % 'numlines' is the total number of lines 
            % in dataTemp (incl newline).
            % 'newline' is [numlines divline].
            altstr1 = ['numlines=$(grep -c "." ' ...
                currData{(i-1)*simInfo.repeat+j} ')' ''];
% %                 altstr1 = ['numlines=$(grep -c "." ' ...
% %                     dataTemp ')' ''];

            %m -n flag returns line number of matched line.
            %m '|cut -f1 -d:' removes all but the line number itself.
            altstr2 = ['divline=$(grep -n -- ''----------------'' ' ...
                currData{(i-1)*simInfo.repeat+j} ' |cut -f1 -d:)'];

            %m Convert to strings
            str4 = 'strdivline=$( printf ''%.0f'' $divline)';
            str5 = 'strnumlines=$(printf ''%.0f'' $numlines)';

            % Concatenate strings
            str6 = 'newline="$strnumlines $strdivline "'; 

            % Now replace first (dummy) line of sim_output.txt with newline
            altstr7 = ['sed -i "1s/.*/$newline/" ' ...
                currData{(i-1)*simInfo.repeat+j} ];


            fprintf(fid,'%s \n %s \n %s \n %s \n %s \n %s \n', ...
            altstr1,altstr2,str3,str4,str5,str6,altstr7);          
            fclose(fid);           



% --------- Begin Section: COPY from /scratch/ (dataTemp) to home folder --

%                 dataOnHome = [fileInfo.dessaFolderOnHome,...
%                     strrep(jobPrefix,'/','_'),'_',num2str(j,jobFormat),...
%                     '_conc',num2str(i,concFormat),'.dat'];
            fidcj = fopen(currJob{(i-1)*simInfo.repeat+j},'a');
            fprintf(fidcj,'cp -bf %s %s\n',...
                currData{(i-1)*simInfo.repeat+j},...
                currDataOnHome{(i-1)*simInfo.repeat+j}); 


            fprintf(fidcj,'cmp %s %s\nexport COPY_CMP=$?\n\n',...
                currData{(i-1)*simInfo.repeat+j},...
                currDataOnHome{(i-1)*simInfo.repeat+j});

            %m LOOP BEGIN 
            fprintf(fidcj,'counter=0\n');
            fprintf(fidcj,'while [ $COPY_CMP -ne 0 ]\n do\n');             
            fprintf(fidcj,'cp -bf %s %s\n',...
                currData{(i-1)*simInfo.repeat+j},...
                currDataOnHome{(i-1)*simInfo.repeat+j});               
            fprintf(fidcj,'cmp %s %s\nexport COPY_CMP=$?\n\n',...
                currData{(i-1)*simInfo.repeat+j},...
                currDataOnHome{(i-1)*simInfo.repeat+j});
            fprintf(fidcj,'((counter++))\n');
            fprintf(fidcj,'if [ $counter -ge 4 ]; then\n');
            fprintf(fidcj,'COPY_CMP=0\n');
            fprintf(fidcj,'fi\n');
            fprintf(fidcj,'done\n');


            %m LOOP END

% %                %m Explanation: 'cmp x y' perfoms a byte wise comparison.
% %                %m Let's say there are 0 differences. Unix reserves $? as a
% %                %m variable name for storing the exit status (of 0). We
% %                %m save this status to a variable we call COPY_CMP.
% %                %m typing 'echo $COPY_CMP' on the commandline would display
% %                %m 0.
% %                %m LOOP END

            fprintf(fidcj,'rm -f %s\n',currData{(i-1)*simInfo.repeat+j});                
            fclose(fidcj);
% -------------- End Section -------------------



        else % UPDATE ELSE CASE TO MATCH IF CASE('newline' stuff)
             display('prepareJob.m - [isfield(fileInfo.dessaFolderOnHome) == false] case');
             fid = fopen(currJob{(i-1)*simInfo.repeat+j},'w');

             fprintf(fid,'%s %s %d %f %d %d > %s\n',fileInfo.command,...
                 xmlFile{i},simInfo.interval,simInfo.time,...
                 randList((i-1)*simInfo.repeat+j),...
                 currData{(i-1)*simInfo.repeat+j});

             fclose(fid);
             assert(false == true)
        end

    end
        

     display('Updating jobList and dataList')
     jobList = [jobList,currJob];
     dataList = [dataList,currDataOnHome];

     display('end of prepareJob.m')

end

function [ xmlOut ] = modifyParameter ( xmlInFile, xmlOutFile, bindSiteList, confSwitchList )

    %XML_READ reads xml files and converts them into Matlab's struct tree.
    [xmlIn, xmlRoot] = xml_read(xmlInFile);
    pref.StructItem = false;
    
    if nargin < 3
        
        disp('Not enough parameters; nothing modified.');
        xml_write(xmlOutFile,xmlIn,xmlRoot,pref);
        
    else
    
        xmlOut = modifyBSTime(xmlIn,bindSiteList);
    
        if nargin >= 4
        
            xmlOut = modifyCSTime(xmlOut,confSwitchList);
        
        end

        xml_write(xmlOutFile,xmlOut,xmlRoot,pref);
        
    end
    
end

function [ xmlOut ] = modifyBSTime ( xmlIn, bindSiteList )

    xmlOut = xmlIn;
    bindingSiteType = xmlIn.BindingSiteTypes.BindingSiteType;
    numBindSite = size(bindingSiteType,1);
    numBSModify = size(bindSiteList,1);
    countModify = 0;
    
    for bsm = 1 : numBSModify
        
        bs1 = bindSiteList{bsm,1}; %m i.e. 'bst0a'
        bs2 = bindSiteList{bsm,2}; %m i.e. 'bst0b'
        timeModify = bindSiteList{bsm,3};
        numItemModify = length(timeModify);
        countModify = countModify + numItemModify;
        
        for bst = 1 : numBindSite
            
            currBS = bindingSiteType(bst).ATTRIBUTE.name;
            ptnNum = length(bindingSiteType(bst).Partner);
            
            for ptn = 1 : ptnNum
            
                partner = bindingSiteType(bst).Partner(ptn).ATTRIBUTE.name;
            
                if (strcmpi(currBS,bs1) && strcmpi(partner,bs2)) || ...
                        (strcmpi(currBS,bs2) && strcmpi(partner,bs1))
                
                    if numItemModify >= 1
                    
                        xmlOut.BindingSiteTypes.BindingSiteType(bst). ...
                            Partner(ptn).ATTRIBUTE.bindTime = 10^timeModify(1);
                    
                    end
                
                    if numItemModify >= 2
                    
                        xmlOut.BindingSiteTypes.BindingSiteType(bst). ...
                            Partner(ptn).ATTRIBUTE.breakTime = 10^timeModify(2);
                    
                    end
                
                    if numItemModify >= 3
                    
                        xmlOut.BindingSiteTypes.BindingSiteType(bst). ...
                            Partner(ptn).ATTRIBUTE.fastBindTime = 10^timeModify(3);
                    
                    end
                
                end
                
            end
            
        end
        
    end
%{
    if countModify == 0
            
        disp('No modification made on Bind/break/fastbind time.');
            
    else
            
        disp('Bind/break/fastbind time modified.');    
            
    end
%}    
end

function [ xmlOut ] = modifyCSTime ( xmlIn, confSwitchList )
    
    xmlOut = xmlIn;
    numConformation = size(xmlIn.ConformationalSwitch.ConformationTime,1);
    numCSModify = size(confSwitchList,1);
    
    for csm = 1 : numCSModify
        
        csFrom = confSwitchList{csm,1};
        csTo = confSwitchList{csm,2};
        timeModify = confSwitchList{csm,3};

        for c = 1 : numConformation
            
            currFrom = xmlIn.ConformationalSwitch.ConformationTime(c).ATTRIBUTE.name;
            
            if strcmpi(currFrom,csFrom) && length(timeModify) >= 1
                
                for cst = 1 : size(xmlIn.ConformationalSwitch.ConformationTime(c).List,1)
                    
                    currTo = xmlIn.ConformationalSwitch.ConformationTime(c).List(cst).ATTRIBUTE.name;
                    
                    if strcmpi(currTo,csTo)
                        
                        xmlOut.ConformationalSwitch.ConformationTime(c).List(cst).ATTRIBUTE.time ...
                            = 10^timeModify(1);
                        
                    end
                    
                end
                
            end
            
            if strcmpi(currFrom,csTo) && length(timeModify) == 2
                
                for cst = 1 : size(xmlIn.ConformationalSwitch.ConformationTime(c).List,1)
                    
                    currTo = xmlIn.ConformationalSwitch.ConformationTime(c).List(cst).ATTRIBUTE.name;
                    
                    if strcmpi(currTo,csFrom)
                        
                        xmlOut.ConformationalSwitch.ConformationTime(c).List(cst).ATTRIBUTE.time ...
                            = 10^timeModify(2);
                        
                    end
                    
                end
                
            end
            
        end
        
    end
%{
    if numCSModify == 0
        
        disp('No modification made on Conformational Switch time.');
        
    else
        
        disp('Conformational Switch time modified.');
        
    end
%}    
end

% ------------------ JUNK CODE ---------------

% %             %m Marcus Addition: replaces prev lines by Lu. In fact, I
% %             %m alter much of Lu's code in this file to handle my DLS output
% %             %m instead of DESSA's sim output.
% %             currDLSJob{(i-1)*simInfo.repeat+j} = [fileInfo.workFolder,...
% %                 jobPrefix,'/',num2str(j,jobFormat),'.sh'];               
% %             currDLSData{(i-1)*simInfo.repeat+j} = [fileInfo.workFolder,...
% %                 jobPrefix,'/',num2str(j,jobFormat),'.dat'];            
            %m Marcus Explains:
            %m (i-1)*simInfo.repeat+j takes on the values 
            %m 0,1,2,...,simInfo.repeat,...,
            %m simInfo.repeat+1,...,.....,(concNum*simInfo.repeat + simInfo.repeat)




% %                 dataTemp = [fileInfo.tempFolder,...
% %                     strrep(jobPrefix,'/','_'),'_',num2str(j,jobFormat),...
% %                    '.dat'];


% %                fid = fopen(currDLSJob{(i-1)*simInfo.repeat+j},'w');



% %                 fprintf(fid,'%s %s %d %f %d %d > %s\n',fileInfo.command,...
% %                     xmlFile{i},simInfo.interval,simInfo.time,...
% %                     randList((i-1)*simInfo.repeat+j),maximer,dataTemp);
                %m This runs the simulator. Pipe DESSA+ console output to
                %m dataTemp. 
                %m My DESSA+ (diff from Lu's) has an args (command line)
                %m length and ordering which seems to be incompatible with
                %m Lu's. I'll re-write it below.
                
% %                 % eventsPerPrint ~ simInfo.interval is args[1] in DESSA+
% %                 % maxSimulationTime ~ simInfo.time is args[2] in DESSA+
% %                 k = 2.56e-7; %m args[3] in DESSA+
% %                 c = 145e-6;  %m args[4] in DESSA+
% %                 m = 250e3;   %m args[5] in DESSA+
% %                 % seedUsed ~ randList(..) is args[6] in DESSA+
% %                 % maxOutputSize ~ maximer is args[7] in DESSA+
% %                 fprintf(fid,'%s %s %d %f %d %d %d %d %d > %s\n',fileInfo.command,...
% %                     xmlFile{i},simInfo.interval,simInfo.time,...
% %                     k,c,m,randList((i-1)*simInfo.repeat+j),maximer,dataTemp);                


% %                 fprintf(fid,'cp -bf %s %s\n',dataTemp,...
% %                     currDLSData{(i-1)*simInfo.repeat+j});
                %m Explanation: copy byte-wise the contents of dataTemp to 
                %m currDLSData{..}.


% %                 fid = fopen(currDLSJob{(i-1)*simInfo.repeat+j},'w');


% %                 fprintf(fid,'cmp %s %s\nexport COPY_CMP=$?\n\n',...
% %                     dataTemp,currDLSData{(i-1)*simInfo.repeat+j});


% %                 fprintf(fid,'cp -bf %s %s\n',dataTemp,...
% %                     currDLSData{(i-1)*simInfo.repeat+j});


% %                 fprintf(fid,'cmp %s %s \n export COPY_CMP=$? \n done \n\n',...
% %                     dataTemp,currDLSData{(i-1)*simInfo.repeat+j});



% %                 fprintf(fid,'%s %s %d %f %d %d > %s\n',fileInfo.command,...
% %                     xmlFile{i},simInfo.interval,simInfo.time,...
% %                     randList((i-1)*simInfo.repeat+j),...
% %                     currDLSData{(i-1)*simInfo.repeat+j});
