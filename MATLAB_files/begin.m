function [] = begin(homeF, scratchF, prefx)
display('in Begin.m')
% % Input arguments --Altered by MARCUS
%optional items can be left blank or commented.

timerVal = tic;
%fileInfo.xmlRule = '/home/marcust/simulator/sim_hpv_init2_grid1_1__.xml';
%fileInfo.xmlRule = '/home/marcust/simulator/CCMV1_156_00.xml';    %MT:single capsid
fileInfo.xmlRule = '/home/marcust/simulator/CCMV_156_00_5capsids.xml';

fileInfo.command = ['java -Xmx5120M -jar '...
    '/home/marcust/simulator/DESSA_2017_SAXS_ReducedTimes_noLikelihoodPrintouts.jar'];

fileInfo.workFolder = scratchF

fileInfo.homeWorkFolder = homeF


fileInfo.dessaFolderOnHome = [homeF,'dessaFolderOnHome/']
A = exist(fileInfo.dessaFolderOnHome,'dir');
display('dessaFolderOnHome/ should exist. I.e., A=7. What is A? ')
display(A);
assert(A==7);

fileInfo.errFolder = [fileInfo.homeWorkFolder,'err/'];
%folder for error information
E = exist(fileInfo.errFolder,'dir');
assert(E==7);


fileInfo.tempFolder = fileInfo.workFolder;

fileInfo.log = [fileInfo.homeWorkFolder,'ccmv_PI.log'];
%log file. Remnant of prev version of code. 


fileInfo.debugFile1 = [fileInfo.homeWorkFolder,'debug1'];
fileInfo.debugFile2 = [fileInfo.homeWorkFolder,'debug2'];
fileInfo.debugFile3 = [fileInfo.homeWorkFolder,'debug3'];
fileInfo.debugFile4 = [fileInfo.homeWorkFolder,'debug4'];
%debug files for different level of debugging (remnant of prev version of code).


noise = -1;%0.025; % real 3-SLS, 50 repeat
%background noise level for score calculation.
%mt: noise <0 tells snobfit to estimate the actual noise using machine prec

% Remnant from earlier version of pipeline.
qInfo = []; 

simInfo.time = 20;
%simulation time as an argument to DESSA.

simInfo.dessaTimesUsed = [.01:.2:16,16.4:.6:20];%,21:300]; 
% When Dessa set to output only specific time points. 
% THIS MUST MATCH THE TIME POINTS DEFINED IN DESSA (file: src/selfassembly/Test.java)

simInfo.repeat = 300;%100;    %For creating ground truth experiment: 2500 or 5000;
%how many trajectories at each parameter point.

simInfo.interval = 1;
%events interval for DESSA output. This should always be 1 with current
%version of DESSA.

simInfo.maximer = 100;


simInfo.SAXScurves = cell(1,1);
simInfo.SAXS_File = [fileInfo.homeWorkFolder,'SynSAXS_CCMV_156_00_',...
    num2str(simInfo.repeat),'_repeats_5capsidXML.mat'];
%m For saving SLS curves which need to be averaged for visualization. 
%m OPTIONAL. Uncomment here if not needed. Check fitMultiCurve.m or
% fitMultiCurveSAXS.m for use.

load /home/marcust/simulator/crysol_Iq
simInfo.crysolIq = {crysol_Iq_cc}; %{crysol_Iq_cc,crysol_Iq_ab}
simInfo.lengthIq = length(crysol_Iq_cc); 
% number of points in Crysol I(q) curve. (q = 0:.01:.5)
% for use in fitMultiCurveSAXS.m
% option: crysol_Iq.mat also contains crysol_Iq_ab. 
% Use in VirusGen_forClusterSAXS.m if needed.

% % % load /home/marcust/simulator/hpv
% % % simInfo.data{1} = hpv145;
% % % simInfo.data{2} = hpv198;
% % % simInfo.data{3} = hpv220;

load /home/marcust/simulator/SynSAXS_CCMV_156_00_2500_repeats_5capsidXML.mat
simInfo.data = simInfo_synthetic.SAXScurves{2};
%MT: This is the ground truth experimental data created by averaging the
%    SAXS curves of 2500 repeated simulations at the GT parameter vector.

simInfo.conc = 1;%[0.53,0.72,0.8]; % 145,198,220 equivalent, diff [units].
%1; %[5.4; 8.2; 10.8];
%MT: CURRENT pipeline does not consider multiple concentrations.
%concentrations of coat proteins in experiment. The absolute unit was not
%important; but the unit must be consistent. Must be column vector.

% simInfo.K = 7.04e-8;%1;%[];
%scaling factor of Rayleigh ratio combined with concentration and molecular
%weight

%%%% ----------------------------------------------------------------------
%%%% -------------- GROUND TRUTH PARAMETERS DEFINED BELOW -----------------


%{
optInfo.BS = {'bst0a' 'bst0b' [1.963204407833257,-4.882236855678068]
              'bst0c' 'bst1a' [1.963204407833257,-4.952929196273922]
              'bst0d' 'bst1b' [1.963204407833257,-4.948210637772643]
              'bst1c' 'bst1d' [1.963204407833257,-5.519189689266716]};
%}
% % THIS IS FOR HBV
% % optInfo.BS = {'bst0a' 'bst0b' [log(91.8764758528154)/log(10),log(7.35886889953384e-06)/log(10)]
% %               'bst0c' 'bst1a' [log(91.8764758528154)/log(10),log(1.171150717264e-05)/log(10)]
% %               'bst0d' 'bst1b' [log(91.8764758528154)/log(10),log(1.171150717264e-05)/log(10)]
% %               'bst1c' 'bst1d' [log(91.8764758528154)/log(10),log(9.58960364746733e-06)/log(10)]};

% % THIS IS FOR CCMV
%{          
optInfo.BS = {'bst0a' 'bst1b' [5.2,3.17226e-05]
              'bst0a' 'bst1d' [5.2,3.17226e-05]
              'bst0b' 'bst1a' [5.2,3.17226e-05]
              'bst0b' 'bst1c' [5.2,3.17226e-05]
              'bst0c' 'bst0d' [5.2,3.17226e-05]};
%}

% % optInfo.BS = {'bst0a' 'bst1b' [log(158489)/log(10),log(1.00007)/log(10)]
% %               'bst0a' 'bst1d' [log(158489)/log(10),log(1.00007)/log(10)]
% %               'bst0b' 'bst1a' [log(158489)/log(10),log(1.00007)/log(10)]
% %               'bst0b' 'bst1c' [log(158489)/log(10),log(1.00007)/log(10)]
% %               'bst0c' 'bst0d' [log(158489)/log(10),log(1.00007)/log(10)]};
optInfo.BS = {'bst0a' 'bst1b' [log(5.1248)/log(10),log(2.87864e-05)/log(10)]
              'bst0a' 'bst1d' [log(5.1248)/log(10),log(2.87864e-05)/log(10)]
              'bst0b' 'bst1a' [log(5.1248)/log(10),log(3.174776e-05)/log(10)]
              'bst0b' 'bst1c' [log(5.1248)/log(10),log(3.174776e-05)/log(10)]
              'bst0c' 'bst0d' [log(5.1248)/log(10),log(3.463432e-05)/log(10)]};
% % optInfo.BS = {'bst0a' 'bst1b' [log(5.2)/log(10),log(3.17226e-05)/log(10)]
% %               'bst0a' 'bst1d' [log(5.0)/log(10),log(1.17226e-05)/log(10)]
% %               'bst0b' 'bst1a' [log(4)/log(10),log(6.17226e-05)/log(10)]
% %               'bst0b' 'bst1c' [log(6)/log(10),log(3.17226e-05)/log(10)]
% %               'bst0c' 'bst0d' [log(5.2)/log(10),log(9.17226e-05)/log(10)]};
% % optInfo.BS = {'bst0a' 'bst1b' [log(15.2)/log(10),log(3.17226e-03)/log(10)]
% %               'bst0a' 'bst1d' [log(25.0)/log(10),log(1.17226e-03)/log(10)]
% %               'bst0b' 'bst1a' [log(1)/log(10),log(3.17226e-04)/log(10)]
% %               'bst0b' 'bst1c' [log(6)/log(10),log(3.17226e-04)/log(10)]
% %               'bst0c' 'bst0d' [log(5.2)/log(10),log(3.17226e-04)/log(10)]};
% % optInfo.BS = {'bst0a' 'bst1b' [log(152)/log(10),log(.003478)/log(10)]
% %               'bst0a' 'bst1d' [log(131.2)/log(10),log(.0011246)/log(10)]
% %               'bst0b' 'bst1a' [log(10)/log(10),log(.0006172)/log(10)]
% %               'bst0b' 'bst1c' [log(15.49)/log(10),log(.0003709)/log(10)]
% %               'bst0c' 'bst0d' [log(52)/log(10),log(.0002063)/log(10)]};
% % optInfo.BS = {'bst0a' 'bst1b' [log(34.2)/log(10),log(1.17226e-06)/log(10)]
% %               'bst0a' 'bst1d' [log(34.2)/log(10),log(1.17226e-06)/log(10)]
% %               'bst0b' 'bst1a' [log(34.2)/log(10),log(1.17226e-06)/log(10)]
% %               'bst0b' 'bst1c' [log(34.2)/log(10),log(1.17226e-06)/log(10)]
% %               'bst0c' 'bst0d' [log(34.2)/log(10),log(1.17226e-06)/log(10)]};
% % optInfo.BS = {'bst0a' 'bst1b' [log(10)/log(10),log(3.17226e-04)/log(10)]
% %               'bst0a' 'bst1d' [log(10)/log(10),log(3.17226e-04)/log(10)]
% %               'bst0b' 'bst1a' [log(10)/log(10),log(3.17226e-04)/log(10)]
% %               'bst0b' 'bst1c' [log(10)/log(10),log(3.17226e-04)/log(10)]
% %               'bst0c' 'bst0d' [log(10)/log(10),log(3.17226e-04)/log(10)]};
% % optInfo.BS = {'bst0a' 'bst1b' [log(.352)/log(10),log(3.17226e-04)/log(10)]
% %               'bst0a' 'bst1d' [log(.352)/log(10),log(3.17226e-04)/log(10)]
% %               'bst0b' 'bst1a' [log(.352)/log(10),log(3.17226e-04)/log(10)]
% %               'bst0b' 'bst1c' [log(.352)/log(10),log(3.17226e-04)/log(10)]
% %               'bst0c' 'bst0d' [log(.352)/log(10),log(3.17226e-04)/log(10)]};


% %m This is for HPV. 
%MT: First optInfo.BS Lu's nucLim params. 
%    Second are July25 start. Third are July16 start. 
%    Fourth are July16 found.
%    Fifth are params on p1549 of Lu's SLS paper which I calculated
%    as wait times.
%    Sixth are Sept28 start
%    Seventh are the params from the xml's Lu believes he used in paper.

% % optInfo.BS = {'bst0a' 'bst0b' [log(79.4328)/log(10),log(1.25892541179417e-05)/log(10)]
% %               'bst0c' 'bst0d' [log(125.892)/log(10),log(7.94328234724282e-06)/log(10)]
% %               'bst0e' 'bst0e' [log(199.526)/log(10),log(5.01187233627272e-06)/log(10)]
% %               'bst0f' 'bst1a' [log(316.227)/log(10),log(3.16227766016838e-06)/log(10)]
% %               'bst0f' 'bst1b' [log(316.227)/log(10),log(3.16227766016838e-06)/log(10)]
% %               'bst0f' 'bst1c' [log(316.227)/log(10),log(3.16227766016838e-06)/log(10)]
% %               'bst0f' 'bst1d' [log(316.227)/log(10),log(3.16227766016838e-06)/log(10)]
% %               'bst0f' 'bst1e' [log(316.227)/log(10),log(3.16227766016838e-06)/log(10)]};
% optInfo.BS = {'bst0a' 'bst0b' [log(60)/log(10),log(1.25e-04)/log(10)]
%               'bst0c' 'bst0d' [log(110)/log(10),log(7.943e-05)/log(10)]
%               'bst0e' 'bst0e' [log(90)/log(10),log(5.011e-05)/log(10)]
%               'bst0f' 'bst1a' [log(310)/log(10),log(18.16e-06)/log(10)]
%               'bst0f' 'bst1b' [log(310)/log(10),log(18.16e-06)/log(10)]
%               'bst0f' 'bst1c' [log(310)/log(10),log(18.16e-06)/log(10)]
%               'bst0f' 'bst1d' [log(310)/log(10),log(18.16e-06)/log(10)]
%               'bst0f' 'bst1e' [log(310)/log(10),log(18.16e-06)/log(10)]};
%{
% % optInfo.BS = {'bst0a' 'bst0b' [log(69)/log(10),log(1.25e-04)/log(10)]
% %               'bst0c' 'bst0d' [log(130)/log(10),log(7.943e-05)/log(10)]
% %               'bst0e' 'bst0e' [log(100)/log(10),log(5.011e-05)/log(10)]
% %               'bst0f' 'bst1a' [log(310)/log(10),log(12.16e-06)/log(10)]
% %               'bst0f' 'bst1b' [log(310)/log(10),log(12.16e-06)/log(10)]
% %               'bst0f' 'bst1c' [log(310)/log(10),log(12.16e-06)/log(10)]
% %               'bst0f' 'bst1d' [log(310)/log(10),log(12.16e-06)/log(10)]
% %               'bst0f' 'bst1e' [log(310)/log(10),log(12.16e-06)/log(10)]};
% % optInfo.BS = {'bst0a' 'bst0b' [log(633.8697)/log(10),log(7.9433e-06)/log(10)]
% %               'bst0c' 'bst0d' [log(690.2398)/log(10),log(1.5524e-05)/log(10)]
% %               'bst0e' 'bst0e' [log(916.2205)/log(10),log(5.0119e-05)/log(10)]
% %               'bst0f' 'bst1a' [log(530.884)/log(10),log(1.2162e-06)/log(10)]
% %               'bst0f' 'bst1b' [log(530.884)/log(10),log(1.2162e-06)/log(10)]
% %               'bst0f' 'bst1c' [log(530.884)/log(10),log(1.2162e-06)/log(10)]
% %               'bst0f' 'bst1d' [log(530.884)/log(10),log(1.2162e-06)/log(10)]
% %               'bst0f' 'bst1e' [log(530.884)/log(10),log(1.2162e-06)/log(10)]};
% % optInfo.BS = {'bst0a' 'bst0b' [log(2.91105e6)/log(10),log(1/0.11)/log(10)]
% %               'bst0c' 'bst0d' [log(2.91105e6)/log(10),log(1/0.12)/log(10)]
% %               'bst0e' 'bst0e' [log(2.91105e6)/log(10),log(1/0.13)/log(10)]
% %               'bst0f' 'bst1a' [log(2.91105e6)/log(10),log(1/0.12)/log(10)]
% %               'bst0f' 'bst1b' [log(2.91105e6)/log(10),log(1/0.12)/log(10)]
% %               'bst0f' 'bst1c' [log(2.91105e6)/log(10),log(1/0.12)/log(10)]
% %               'bst0f' 'bst1d' [log(2.91105e6)/log(10),log(1/0.12)/log(10)]
% %               'bst0f' 'bst1e' [log(2.91105e6)/log(10),log(1/0.12)/log(10)]};
% % optInfo.BS = {'bst0a' 'bst0b' [log(9.50660e2)/log(10),log(9.2339e-2)/log(10)]
% %               'bst0c' 'bst0d' [log(11.5e2)/log(10),log(8.2669e-3)/log(10)]
% %               'bst0e' 'bst0e' [log(20.50660e4)/log(10),log(7.9161e-3)/log(10)]
% %               'bst0f' 'bst1a' [log(3.50660e4)/log(10),log(8.0191e-2)/log(10)]
% %               'bst0f' 'bst1b' [log(3.50660e4)/log(10),log(8.0191e-2)/log(10)]
% %               'bst0f' 'bst1c' [log(3.50660e4)/log(10),log(8.0191e-2)/log(10)]
% %               'bst0f' 'bst1d' [log(3.50660e4)/log(10),log(8.0191e-2)/log(10)]
% %               'bst0f' 'bst1e' [log(3.50660e4)/log(10),log(8.0191e-2)/log(10)]};
          
% % optInfo.BS = {'bst0a' 'bst0b' [log(9.50660e5)/log(10),log(9.23391337)/log(10)]
% %               'bst0c' 'bst0d' [log(9.50660e5)/log(10),log(8.26699933)/log(10)]
% %               'bst0e' 'bst0e' [log(9.50660e5)/log(10),log(7.91616814)/log(10)]
% %               'bst0f' 'bst1a' [log(9.50660e5)/log(10),log(8.01919196)/log(10)]
% %               'bst0f' 'bst1b' [log(9.50660e5)/log(10),log(8.01919196)/log(10)]
% %               'bst0f' 'bst1c' [log(9.50660e5)/log(10),log(8.01919196)/log(10)]
% %               'bst0f' 'bst1d' [log(9.50660e5)/log(10),log(8.01919196)/log(10)]
% %               'bst0f' 'bst1e' [log(9.50660e5)/log(10),log(8.01919196)/log(10)]};
%}
%optInfo.CS = {'bs0' 'bs1' 4
%              'bs1' 'bs0' 4};
%conformation switch time for optimization.
%OPTIONAL. You don't have to write this term if you don't want to play with
%conformation switch.

optInfo.initGridSize = 2;%1;
%initial grid size for generating the first parameter candidates.

optInfo.eps = 1e-3;
%threshold for improvement.
%OPTIONAL, default = 1e-3;


optInfo.lb = -1;
optInfo.ub = 1;
%Used with SNOBFIT as the global solver in inLoop2.m


% OBSOLETE: RandStream.setDefaultStream(RandStream('mt19937ar','seed',647456));
RandStream.setGlobalStream(RandStream('mt19937ar','seed',647456));
%set the random seed; OPTIONAL

fileInfo.paramFile = [fileInfo.homeWorkFolder,'paramFile.mat'];

format long

fileInfo.prefix = prefx


% % NOTE: for optInfo.paramList: {'b',1,1} means:
% %           column 1 - type is 'b' for binding not conformation switch.
% %           column 2 - which row(s) of optInfo.BS to consider. row 1
% %           column 3 - which column of optInfo.BS to consider. col 1

% % MT: for HBV
% % optInfo.paraList{1} = {'b',1,1};
% % optInfo.paraList{2} = {'b',2,1};
% % optInfo.paraList{3} = {'b',3,1};
% % optInfo.paraList{4} = {'b',4,1};
% % optInfo.paraList{5} = {'b',1,2};
% % optInfo.paraList{6} = {'b',2,2};
% % optInfo.paraList{7} = {'b',3,2};
% % optInfo.paraList{8} = {'b',4,2};

% % for HPV
% % optInfo.paraList{1} = {'b',1,1};
% % optInfo.paraList{2} = {'b',2,1};
% % optInfo.paraList{3} = {'b',3,1};
% % optInfo.paraList{4} = {'b',(4:8),1};
% % optInfo.paraList{5} = {'b',1,2};
% % optInfo.paraList{6} = {'b',2,2};
% % optInfo.paraList{7} = {'b',3,2};
% % optInfo.paraList{8} = {'b',(4:8),2};
%MT: So, we effectively have 8 parameters since binding sites 4:8
%    are all bonds of type A+ ~ A-. A single on-rate and a single 
% off-rate will be applicable to them.

% % optInfo.paraList{1} = {'b',(1:3),1}; %group first 3 on rates
% % optInfo.paraList{2} = {'b',(4:8),1}; %group last 5 on rates
% % optInfo.paraList{3} = {'b',(1:3),2}; %group first 3 off rates
% % optInfo.paraList{4} = {'b',(4:8),2}; %group last 5 off rates
%MT Jan.29 2016 - only 4 effective parameters.

% CCMV
% optInfo.paraList{1} = {'b',(1:5),1}; %group all on rates
% optInfo.paraList{2} = {'b',(1:5),2}; %group all off rates

optInfo.paraList{1} = {'b',1,1}; 
optInfo.paraList{2} = {'b',2,1}; 
optInfo.paraList{3} = {'b',3,1};
optInfo.paraList{4} = {'b',4,1};
optInfo.paraList{5} = {'b',5,1};
optInfo.paraList{6} = {'b',1,2};
optInfo.paraList{7} = {'b',2,2};
optInfo.paraList{8} = {'b',3,2};
optInfo.paraList{9} = {'b',4,2};
optInfo.paraList{10} = {'b',5,2};


%%%% -------- These are remnants from earlier version of pipeline ---------

% The order of cells of optInfo.likelihoods and optInfo.traj_matrices 
% corresponds to the ordering of optInfo.score and optInfo.offset 
% in inLoop2.m.  E.g. the first cell of optInfo.likelihoods contains the
% the trajectory likelihoods for all repeats corresponding to the first
% entry of optInfo.offset.
optInfo.likelihoods = {};
optInfo.traj_matrices = {};
optInfo.RMSDs = {};

optInfo.Kernel = [];
optInfo.K_blocks = [];
%%%% ----------------------------------------------------------------------



save(fileInfo.paramFile,'fileInfo','qInfo','simInfo','optInfo');

display('post save, begin.m')
end
