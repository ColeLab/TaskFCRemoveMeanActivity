function [output] = Preproc_HCPData_PostMinPreproc_TaskRest_BasisHRFModel(SUBJECTLIST)

% This function runs preprocessing and GLM analysis, after the HCP minimal preproccessing pipeline has been run.
% It is designed to parcellate a dataset into a set of regions, which are then preprocessed.
%
% Preprocessing steps included:
% 1) Parcellate a dense CIFTI file into a set of N time series, where N is the number of parcels/regions
% 2) Prepare nuisance regressors for removal. This includes spatial mask definition (e.g., white matter)
%    and extraction from volume fMRI data, as well as processing of motion parameters.
% 3) Preparation of task regressors, when task runs are present. This includes (with custom scripts)
%    conversion of behavioral timing data to a common format, convolving with a hemodynamic response function, 
%    and conversion to a GLM design matrix.
% 4) Rest fMRI nuisance regression for functional connectivity analyses, if rest data are present.
% 5) Task fMRI GLM along with nuisance regression (if task data are present), for functional connectivity and/or task activation analyses.
% 6) Temporal filtering of time series (optional).
%
% Note: Frequently customized variables are IN CAPS throughout the script

%Script author:
%Michael W. Cole, mwcole@mwcole.net, http://www.colelab.org
%
%Script version: 1.5
%Date: October 25, 2016
%
%Version history:
%1.1: Fixed a bug in which temporal filtering would run even if flagged not
%   to. Also changed standard TR duration.
%1.2: Added check for GLM design matrix deficiency
%1.3: NPROC now works properly. Added information on deleting intermediate 
%   temporary files to save disk space. Changed which directory files are 
%   saved to (now analysis_[ANALYSISNAME]). Now saves out 
%   output_GLM.taskdesignmat_hrf_tmasked, which indicates task time points 
%   after scrubbing.
%1.4: Fixed bug in which temporal filtering was only run on the last 
%   subject (instead of all subjects). Changed final output file to v7.3 to
%   accomodate large file sizes.
%1.5: Fixed a bug in which the hemispheres might have been flipped. Now each
%   hemisphere is loaded separately to ensure they are not flipped.

%% Parameters to customize for your analysis
addpath('/projects/AnalysisTools/')
addpath('/projects/AnalysisTools/gifti-1.6/')

%%Basic processing parameters
ANALYSISNAME='TaskFCMethods_HCPData_HRFBasis';
mfilenameRep='Preproc_PostMinPreproc_TaskAndRest_BasisHRFModel';
%subjList = list of subject numbers as a list of strings (used if not set by function call)
if isempty(SUBJECTLIST)
    SUBJECTLIST = {'100307','100408','101107','101309','101915','103111','103414','103818','105014','105115','106016','108828','110411','111312','111716','113619','113922','114419','115320','116524','117122','118528','118730','118932','120111','122317','122620','123117','123925','124422','125525','126325','127630','127933','128127','128632','129028','130013','130316','131217','131722','133019','133928','135225','135932','136833','138534','139637','140925','144832','146432','147737','148335','148840','149337','149539','149741','151223','151526','151627','153025','154734','156637','159340','160123','161731','162733','163129','176542','178950','188347','189450','190031','192540','196750','198451','199655','201111','208226','211417','211720','212318','214423','221319','239944','245333','280739','298051','366446','397760','414229','499566','654754','672756','751348','756055','792564','856766','857263','899885'};
    %Full subject list: {'100307','100408','101107','101309','101915','103111','103414','103818','105014','105115','106016','108828','110411','111312','111716','113619','113922','114419','115320','116524','117122','118528','118730','118932','120111','122317','122620','123117','123925','124422','125525','126325','127630','127933','128127','128632','129028','130013','130316','131217','131722','133019','133928','135225','135932','136833','138534','139637','140925','144832','146432','147737','148335','148840','149337','149539','149741','151223','151526','151627','153025','154734','156637','159340','160123','161731','162733','163129','176542','178950','188347','189450','190031','192540','196750','198451','199655','201111','208226','211417','211720','212318','214423','221319','239944','245333','280739','298051','366446','397760','414229','499566','654754','672756','751348','756055','792564','856766','857263','899885'};
end
numSubjs=length(SUBJECTLIST);
TR_INSECONDS=0.720;
FRAMESTOSKIP=5;

%Basic data parameters
RUNNAMES = {'rfMRI_REST1_RL', 'rfMRI_REST1_LR', 'rfMRI_REST2_RL', 'rfMRI_REST2_LR', 'tfMRI_EMOTION_RL','tfMRI_EMOTION_LR','tfMRI_GAMBLING_RL','tfMRI_GAMBLING_LR','tfMRI_LANGUAGE_RL','tfMRI_LANGUAGE_LR','tfMRI_MOTOR_RL','tfMRI_MOTOR_LR','tfMRI_RELATIONAL_RL','tfMRI_RELATIONAL_LR','tfMRI_SOCIAL_RL','tfMRI_SOCIAL_LR','tfMRI_WM_RL','tfMRI_WM_LR'};
%Rest runs: 'rfMRI_REST1_RL', 'rfMRI_REST1_LR', 'rfMRI_REST2_RL', 'rfMRI_REST2_LR'
%Task runs: 'tfMRI_EMOTION_RL','tfMRI_EMOTION_LR','tfMRI_GAMBLING_RL','tfMRI_GAMBLING_LR','tfMRI_LANGUAGE_RL','tfMRI_LANGUAGE_LR','tfMRI_MOTOR_RL','tfMRI_MOTOR_LR','tfMRI_RELATIONAL_RL','tfMRI_RELATIONAL_LR','tfMRI_SOCIAL_RL','tfMRI_SOCIAL_LR','tfMRI_WM_RL','tfMRI_WM_LR'
numRuns=length(RUNNAMES);
RESTRUNS=1:4;
TASKRUNS=5:18;
RUNLENGTHS = [1200, 1200, 1200, 1200, 176,176,253,253,316,316,284,284,232,232,274,274,405,405];
%Num TRs for rest runs: 1200, 1200, 1200, 1200
%Num TRs for task runs: 176,176,253,253,316,316,284,284,232,232,274,274,405,405
L_parcelCIFTIFile='/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii';
R_parcelCIFTIFile='/projects/AnalysisTools/ParcelsGlasser2016/Q1-Q6_RelatedParcellation210.R.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii';
parcellationName='Glasser2016';
NUMPARCELS=360;

%Data processing flags
GSR=0;      %GSR = 0 if no GSR, 1 if you want to include GSR
NPROC=4;    %Number of processors to use when computing GLMs
IMPLEMENT_MOTIONSCRUBBING=0;     %Set to 1 for yes, 0 for no
FDTHRESH=0.25;  %Framewise displacement; The threshold (in millimeters) for flagging in-scanner movement
TEMPORALFILTER=0;   %Flag indicating if data should be temporally filtered
HIGHPASS_HZ=0.008;
LOWPASS_HZ=0.09;

%Directories
BASEDIR='/projects/ExternalDatasets/HCPData/';
datadir=[BASEDIR '/data/minimalpreproc/'];
outputdatadir='/projects3/TaskFCMethods/data/';
timingfileDir=['/projects3/TaskFCMethods/data/' '/timingfiles/HCPTiming/'];
if ~exist(timingfileDir, 'dir'); mkdir(timingfileDir); end


%% Iterate through subjects

output=[];
output.SUBJECTLIST=SUBJECTLIST;
output.RUNNAMES=RUNNAMES;
output.RUNLENGTHS=RUNLENGTHS;
output.RESTRUNS=RESTRUNS;
output.TASKRUNS=TASKRUNS;
output.parcellationName=parcellationName;

for subjIndex=1:numSubjs
    subjNum = SUBJECTLIST{subjIndex};
    disp(['Processing subject ' subjNum]);
    
    tseriesMatSubj=zeros(NUMPARCELS,max(RUNLENGTHS),numRuns);
    
    %Directories
    subjDir=[datadir '/' subjNum '/'];
    SUBJDIROUTPUT=['/projects3/TaskFCMethods/data/' subjNum '/'];   %Typically set to be same as subjDir
    if ~exist(SUBJDIROUTPUT, 'dir'); mkdir(SUBJDIROUTPUT); end
    subjTemporaryAnalysisDir=[SUBJDIROUTPUT '/analysis_' ANALYSISNAME '/'];
    if ~exist(subjTemporaryAnalysisDir, 'dir'); mkdir(subjTemporaryAnalysisDir); end
    
    
    %% Downsampling grayordinate data to parcels
    
    disp('Downsampling grayordinate data to parcels')
    
    %Set to 1 if you want to run this procedure again for subjects that already had it run before (otherwise it will load from previously saved files)
    RERUN_EXTRACTPARCELS=0;
    
    runCount=1;
    for runName=RUNNAMES
        thisRunName=runName{1};
        
        L_parcelTSFilename=[subjTemporaryAnalysisDir '/' thisRunName '_Atlas.L.' parcellationName 'Parcels.32k_fs_LR.ptseries.nii'];
        R_parcelTSFilename=[subjTemporaryAnalysisDir '/' thisRunName '_Atlas.R.' parcellationName 'Parcels.32k_fs_LR.ptseries.nii'];
        
        if or(RERUN_EXTRACTPARCELS == 1, ~exist(L_parcelTSFilename, 'file'))
            
            subjRunDir=[subjDir '/MNINonLinear/Results/' thisRunName '/'];
            
            inputFile=[subjRunDir thisRunName '_Atlas.dtseries.nii'];
            
            eval(['!wb_command -cifti-parcellate ' inputFile ' ' L_parcelCIFTIFile ' COLUMN ' L_parcelTSFilename ' -method MEAN'])
            eval(['!wb_command -cifti-parcellate ' inputFile ' ' R_parcelCIFTIFile ' COLUMN ' R_parcelTSFilename ' -method MEAN'])
            
        end
        
        %Load parcellated data
        L_dat = ciftiopen(L_parcelTSFilename,'wb_command');
        R_dat = ciftiopen(R_parcelTSFilename,'wb_command');
        if size(L_dat.cdata,2)>RUNLENGTHS(runCount)
            disp(['WARNING: More TRs for this run than expected. Subject: ' subjNum ', Run: ' num2str(runCount)])
        end
        tseriesMatSubj(1:180,1:RUNLENGTHS(runCount),runCount)=L_dat.cdata(:,1:RUNLENGTHS(runCount));
        tseriesMatSubj(181:end,1:RUNLENGTHS(runCount),runCount)=R_dat.cdata(:,1:RUNLENGTHS(runCount));
        
        runCount=runCount+1;
    end
    

    
    %% Prepare nuisance regressors
    disp('Preparing nuisance regressors')
    %Using Freesurfer aparc+aseg masks
    
    %Set to 1 if you want to run this procedure again for subjects that already had it run before (otherwise it will load from previously saved files)
    RERUN_PREPNUISANCEREG=0;
    
    savedfile=[subjTemporaryAnalysisDir mfilenameRep '_' ANALYSISNAME '_nuisanceTSVars.mat'];
    
    if or(RERUN_PREPNUISANCEREG == 1, ~exist(savedfile, 'file'))
        
        subjMaskDir=[subjDir '/masks/'];
        if ~exist(subjMaskDir, 'dir'); mkdir(subjMaskDir); end
        
        %Resample Freesurfer segmented mask into functional space (using nearest neighbor interpolation); uses AFNI
        subjMNINonLinearDir=[subjDir '/MNINonLinear/'];
        exampleFunctionalVolFile=[subjDir '/MNINonLinear/Results/' RUNNAMES{1} '/' RUNNAMES{1} '.nii.gz'];
        eval(['!3dresample -overwrite -rmode NN -master ' exampleFunctionalVolFile ' -inset ' subjMNINonLinearDir 'aparc+aseg.nii.gz -prefix ' SUBJDIROUTPUT 'aparc+aseg_resampFunc.nii.gz']);
        
        %Load Freesurfer segmented mask
        aparc_aseg=load_nifti([SUBJDIROUTPUT '/aparc+aseg_resampFunc.nii.gz']);
        
        %Create gray matter mask
        maskValSet_graymatter=[8 9 10 11 12 13 16 17 18 19 20 26 27 28 47 48 49 50 51 52 53 54 55 56 58 59 60 96 97 1000:1035 2000:2035];
        grayMatterMask=ismember(aparc_aseg.vol,maskValSet_graymatter);
        
        %Create white matter mask
        maskValSet_whitematter=[2 7 41 46];
        whiteMatterMask=ismember(aparc_aseg.vol,maskValSet_whitematter);
        %Erode white matter mask by 2 voxels
        whiteMatterMask_eroded=imerode(whiteMatterMask,strel(ones(2,2,2)));
        
        %Create ventricle mask
        maskValSet_ventricles=[4 43 14 15];
        ventricleMask=ismember(aparc_aseg.vol,maskValSet_ventricles);
        %Erode ventricle matter mask by 2 voxels
        ventricleMask_eroded=imerode(ventricleMask,strel(ones(2,2,2)));
        
        %Create whole brain mask
        wholebrainMask=aparc_aseg.vol>0;
        
        %Load in nuisance time series for each run
        nuisanceTS_whitematter=zeros(max(RUNLENGTHS),numRuns);
        nuisanceTS_ventricles=zeros(max(RUNLENGTHS),numRuns);
        nuisanceTS_wholebrain=zeros(max(RUNLENGTHS),numRuns);
        nuisanceTS_motion=zeros(12,max(RUNLENGTHS),numRuns);
        FD_motion=zeros(max(RUNLENGTHS),numRuns);
        temporalMask=ones(max(RUNLENGTHS),numRuns);
        numFramesCensored=zeros(numRuns,1);
        
        runCount=1;
        for runName=RUNNAMES
            thisRunName=runName{1};
            
            subjRunDir=[subjDir '/MNINonLinear/Results/' thisRunName '/'];
            
            inputFile=[subjRunDir thisRunName '.nii.gz'];
            
            %Load data
            runData=load_nifti(inputFile);
            runData2D=reshape(runData.vol,size(runData.vol,1)*size(runData.vol,2)*size(runData.vol,3),size(runData.vol,4));
            if size(runData2D,2)>RUNLENGTHS(runCount)
                disp(['WARNING: More TRs for this run than expected. Subject: ' subjNum ', Run: ' thisRunName ', Run number: ' num2str(runCount)])
                runData2D=runData2D(:,1:RUNLENGTHS(runCount));
            end
            
            whiteMatterMask_eroded_1D=reshape(whiteMatterMask_eroded,size(whiteMatterMask_eroded,1)*size(whiteMatterMask_eroded,2)*size(whiteMatterMask_eroded,3),1);
            nuisanceTS_whitematter(1:RUNLENGTHS(runCount),runCount)=mean(runData2D(whiteMatterMask_eroded_1D,:),1);
            
            ventricleMask_eroded_1D=reshape(ventricleMask_eroded,size(ventricleMask_eroded,1)*size(ventricleMask_eroded,2)*size(ventricleMask_eroded,3),1);
            nuisanceTS_ventricles(1:RUNLENGTHS(runCount),runCount)=mean(runData2D(ventricleMask_eroded_1D,:),1);
            
            wholebrainMask_1D=reshape(wholebrainMask,size(wholebrainMask,1)*size(wholebrainMask,2)*size(wholebrainMask,3),1);
            nuisanceTS_wholebrain(1:RUNLENGTHS(runCount),runCount)=mean(runData2D(wholebrainMask_1D,:),1);
            
            %Note: derivatives are already included in motion time series
            motionvals=importdata([subjRunDir 'Movement_Regressors.txt'])';
            nuisanceTS_motion(:,1:RUNLENGTHS(runCount),runCount)=motionvals(:,1:RUNLENGTHS(runCount));
            
            %Skip first FRAMESTOSKIP frames
            temporalMask(1:FRAMESTOSKIP,runCount)=0;
            
            %Calculate framewise displacement (FD) according to Power et al. (2012)
            %Briefly: The sum of the absolute values of the translational and rotational displacements over all frames (in mm)
            %Note: HCP's minimal preprocessing pipeline uses the following ordering of the motion parameters (see https://github.com/Washington-University/Pipelines/blob/master/global/scripts/mcflirt_acc.sh):
            %trans x, trans y, trans z, rot x, rot y, rot z [rotations in degrees], then derivatives of those 6 (for 12 total)
            motionTS_dt=[zeros(size(nuisanceTS_motion,1),1) diff(squeeze(nuisanceTS_motion(:,1:RUNLENGTHS(runCount),runCount))')'];
            assumedRadius=50;
            rot_x=(2*assumedRadius*pi/360)*motionTS_dt(4,:);
            rot_y=(2*assumedRadius*pi/360)*motionTS_dt(5,:);
            rot_z=(2*assumedRadius*pi/360)*motionTS_dt(6,:);
            FD_motion(1:RUNLENGTHS(runCount),runCount)=abs(motionTS_dt(1,:))+abs(motionTS_dt(2,:))+abs(motionTS_dt(3,:))+abs(rot_x)+abs(rot_y)+abs(rot_z);
            %Apply temporal filtering to FD, to reduce the effect of respiration on FD measure. Based on Siegel et al. (2016) [Siegel JS, Mitra A, Laumann TO, Seitzman BA, Raichle M, Corbetta M, Snyder AZ (2016) ?Data Quality Influences Observed Links Between Functional Connectivity and Behavior?. Cereb Cortex. 1?11.http://doi.org/10.1093/cercor/bhw253]
            %Create temporal filter
            lopasscutoff=0.3/(0.5/TR_INSECONDS); % lowpass filter of 0.3 Hz
            %Using filter order of 1
            filtorder=1;
            [butta, buttb]=butter(filtorder,lopasscutoff);
            %Apply temporal filter
            filteredFD=filtfilt(butta,buttb,FD_motion(1:RUNLENGTHS(runCount),runCount));
            FD_motion(1:RUNLENGTHS(runCount),runCount)=filteredFD;
            
            %Implement motion scrubbing/censoring in the temporal mask
            if IMPLEMENT_MOTIONSCRUBBING
                threshedMotion=FD_motion(1:RUNLENGTHS(runCount),runCount)<FDTHRESH;
                %Mark one frame before high motion (consistent with Power et al., 2012)
                oneFrameBefore=find(~threshedMotion)-1;
                oneFrameBefore=oneFrameBefore(oneFrameBefore>0);
                threshedMotion(oneFrameBefore)=0;
                %Mark two frame after high motion (consistent with Power et al., 2012)
                twoFramesAfter=[find(~threshedMotion)+1 find(~threshedMotion)+2];
                twoFramesAfter=twoFramesAfter(twoFramesAfter<RUNLENGTHS(runCount));
                threshedMotion(twoFramesAfter)=0;
                %Add high motion points to temporal mask
                temporalMask(1:RUNLENGTHS(runCount),runCount)=temporalMask(1:RUNLENGTHS(runCount),runCount).*threshedMotion;
                percentTimePointsCensored=100*sum(~threshedMotion)/RUNLENGTHS(runCount);
                disp(['Marking high motion time points. ' num2str(percentTimePointsCensored) '% of time points marked for censoring/scrubbing for this run.'])
                numFramesCensored(runCount)=sum(~threshedMotion);
            end
            
            runCount=runCount+1;
        end
        
        if IMPLEMENT_MOTIONSCRUBBING
            nuisanceTSVars.numFramesScrubbedByRun=numFramesCensored;
            nuisanceTSVars.percentFramesScrubbed=100*sum(numFramesCensored)/sum(RUNLENGTHS);
        end
        
        %Organize nuisance regressors
        nuisanceTSVars.nuisanceTS_whitematter=nuisanceTS_whitematter;
        nuisanceTSVars.nuisanceTS_ventricles=nuisanceTS_ventricles;
        nuisanceTSVars.nuisanceTS_wholebrain=nuisanceTS_wholebrain;
        nuisanceTSVars.nuisanceTS_motion=nuisanceTS_motion;
        nuisanceTSVars.FD_motion=FD_motion;
        nuisanceTSVars.temporalMask=temporalMask;
        
        %Save nuisance time series to file (for more efficient processing in future when EXECUTE_PREPNUISANCEREG=0)
        disp(['Saving results to: ' savedfile])
        save(savedfile, 'nuisanceTSVars')
        
    else
        
        %Load nuisance time series from file (for more efficient processing when EXECUTE_PREPNUISANCEREG=0)
        disp(['Loading results from: ' savedfile])
        load(savedfile);        
        
    end
    
    output.nuisanceTSVars{subjIndex}=nuisanceTSVars;
    
    
    
    %% Prepare task regressors
    
    %Set to 1 if you want to run this procedure again for subjects that already had it run before (otherwise it will load from previously saved files)
    RERUN_TASKREGRESSORPREP=1;
    
    savedfile=[subjTemporaryAnalysisDir mfilenameRep '_' ANALYSISNAME '_TaskRegressorVars.mat'];
    
    %Only execute this if EXECUTE_TASKREGRESSORPREP==1 or the savedfile doesn't exist for this subject, but skip if no TASKRUNS (i.e., only rest data are included)
    if and(~isempty(TASKRUNS), or(RERUN_TASKREGRESSORPREP == 1, ~exist(savedfile, 'file')))
        
        %% PROJECT-SPECIFIC CUSTOM SCRIPTS THAT IMPORT BEHAVIORAL DATA INFO INTO MATLAB
        % The end result should be a "reference function"; a list of strings labeling each time point (TR) with behavioral info
        
        disp('Obtaining regressors for each task event');
        NUMCONDS=24;
        regressorMatrix=zeros(sum(RUNLENGTHS(TASKRUNS)),NUMCONDS);
        regressorMatrix(:,1)=importdata([timingfileDir num2str(subjNum) '_stimfile_HCPBasicBlockFactor_EV1_EMOTION_fear.1D']);
        regressorMatrix(:,2)=importdata([timingfileDir num2str(subjNum) '_stimfile_HCPBasicBlockFactor_EV2_EMOTION_neut.1D']);
        regressorMatrix(:,3)=importdata([timingfileDir num2str(subjNum) '_stimfile_HCPBasicBlockFactor_EV3_GAMBLINGwin.1D']);
        regressorMatrix(:,4)=importdata([timingfileDir num2str(subjNum) '_stimfile_HCPBasicBlockFactor_EV4_GAMBLINGloss.1D']);
        regressorMatrix(:,5)=importdata([timingfileDir num2str(subjNum) '_stimfile_HCPBasicBlockFactor_EV5_LANGUAGEstory.1D']);
        regressorMatrix(:,6)=importdata([timingfileDir num2str(subjNum) '_stimfile_HCPBasicBlockFactor_EV6_LANGUAGEmath.1D']);
        regressorMatrix(:,7)=importdata([timingfileDir num2str(subjNum) '_stimfile_HCPBasicBlockFactor_EV7_MOTORcue.1D']);
        regressorMatrix(:,8)=importdata([timingfileDir num2str(subjNum) '_stimfile_HCPBasicBlockFactor_EV8_MOTORlf.1D']);
        regressorMatrix(:,9)=importdata([timingfileDir num2str(subjNum) '_stimfile_HCPBasicBlockFactor_EV9_MOTORrf.1D']);
        regressorMatrix(:,10)=importdata([timingfileDir num2str(subjNum) '_stimfile_HCPBasicBlockFactor_EV10_MOTORlh.1D']);
        regressorMatrix(:,11)=importdata([timingfileDir num2str(subjNum) '_stimfile_HCPBasicBlockFactor_EV11_MOTORrh.1D']);
        regressorMatrix(:,12)=importdata([timingfileDir num2str(subjNum) '_stimfile_HCPBasicBlockFactor_EV12_MOTORt.1D']);
        regressorMatrix(:,13)=importdata([timingfileDir num2str(subjNum) '_stimfile_HCPBasicBlockFactor_EV13_RELATIONALrelation.1D']);
        regressorMatrix(:,14)=importdata([timingfileDir num2str(subjNum) '_stimfile_HCPBasicBlockFactor_EV14_RELATIONALmatch.1D']);
        regressorMatrix(:,15)=importdata([timingfileDir num2str(subjNum) '_stimfile_HCPBasicBlockFactor_EV15_SOCIALmental.1D']);
        regressorMatrix(:,16)=importdata([timingfileDir num2str(subjNum) '_stimfile_HCPBasicBlockFactor_EV16_SOCIALrnd.1D']);
        regressorMatrix(:,17)=importdata([timingfileDir num2str(subjNum) '_stimfile_HCPBasicBlockFactor_EV17_WM0bk_body.1D']);
        regressorMatrix(:,18)=importdata([timingfileDir num2str(subjNum) '_stimfile_HCPBasicBlockFactor_EV18_WM0bk_faces.1D']);
        regressorMatrix(:,19)=importdata([timingfileDir num2str(subjNum) '_stimfile_HCPBasicBlockFactor_EV19_WM0bk_places.1D']);
        regressorMatrix(:,20)=importdata([timingfileDir num2str(subjNum) '_stimfile_HCPBasicBlockFactor_EV20_WM0bk_tools.1D']);
        regressorMatrix(:,21)=importdata([timingfileDir num2str(subjNum) '_stimfile_HCPBasicBlockFactor_EV21_WM2bk_body.1D']);
        regressorMatrix(:,22)=importdata([timingfileDir num2str(subjNum) '_stimfile_HCPBasicBlockFactor_EV22_WM2bk_faces.1D']);
        regressorMatrix(:,23)=importdata([timingfileDir num2str(subjNum) '_stimfile_HCPBasicBlockFactor_EV23_WM2bk_places.1D']);
        regressorMatrix(:,24)=importdata([timingfileDir num2str(subjNum) '_stimfile_HCPBasicBlockFactor_EV24_WM2bk_tools.1D']);
        taskdesignmat=regressorMatrix;
        
        %Collapse block regressors by task
        numTasks=7;
        regressorMatrix_ByTask=zeros(size(regressorMatrix,1),numTasks);
        %Emotion task
        regressorMatrix_ByTask(:,1)=regressorMatrix(:,1)+regressorMatrix(:,2);
        %Gambling task
        regressorMatrix_ByTask(:,2)=regressorMatrix(:,3)+regressorMatrix(:,4);
        %Language task
        regressorMatrix_ByTask(:,3)=regressorMatrix(:,5)+regressorMatrix(:,6);
        %Motor task
        regressorMatrix_ByTask(:,4)=sum(regressorMatrix(:,7:12),2);
        %Relational task
        regressorMatrix_ByTask(:,5)=sum(regressorMatrix(:,13:14),2);
        %Social task
        regressorMatrix_ByTask(:,6)=sum(regressorMatrix(:,15:16),2);
        %WM task
        regressorMatrix_ByTask(:,7)=sum(regressorMatrix(:,17:24),2);
        
        %Move resulting timing files to timing files directory (if saving to AFNI .1D files)
        %eval(['!mv *stimfile*.1D ' timingfileDir]);
        %eval(['!mv *ConcatRef*.1D ' timingfileDir]);
        
        disp('Preparing HRF basis functions for the GLM')
        %Define HRF basis functions (obtained from FSL's FLOBS); sampled at 50ms resolution (20 Hz)
        numBasisFuncs=5;
        numBasisFuncTimepoints=414;
        samplingRateOfHRF=20;   %In Hz
        HRFBasisFuncs=zeros(numBasisFuncTimepoints,numBasisFuncs);
        HRFBasisFuncs(:,1)=[1.51E-07 1.47E-06 5.68E-06 1.47E-05 3.04E-05 5.39E-05 8.65E-05 0.000129654 0.000184896 0.000253506 0.000336978 0.000436366 0.000552952 0.00068782 0.000842129 0.001016785 0.001212671 0.001430859 0.00167222 0.001937615 0.002227943 0.002544423 0.002888172 0.003260602 0.003662795 0.00409579 0.004560572 0.005058326 0.005589911 0.006156055 0.006757744 0.007395931 0.008071811 0.008786051 0.009539181 0.010331838 0.011164482 0.012037122 0.01295003 0.013903195 0.014896861 0.015930278 0.017002178 0.018111268 0.019256204 0.020435608 0.021648066 0.022892116 0.024166277 0.025469035 0.026798845 0.028154146 0.029533331 0.030934794 0.032356904 0.033798011 0.035256463 0.036730573 0.038218667 0.039719058 0.041230055 0.042749973 0.044277111 0.045809791 0.047346439 0.048885683 0.050426045 0.051966021 0.053504083 0.055038629 0.056568271 0.058091623 0.059607314 0.061114132 0.062610739 0.064095806 0.065568138 0.067026487 0.068469631 0.069896513 0.071305915 0.072696579 0.07406733 0.075417067 0.076744856 0.078049651 0.07933052 0.080586689 0.081817308 0.083021504 0.084198778 0.085348518 0.086470186 0.087562965 0.088626072 0.089659002 0.090661287 0.091632317 0.092571638 0.093478638 0.094352985 0.095194439 0.09600285 0.096778009 0.097519613 0.098227298 0.098900881 0.099540201 0.100144816 0.100714481 0.101248979 0.101748117 0.102211867 0.102640272 0.103033262 0.103390897 0.103713192 0.10400003 0.104251473 0.104467461 0.104648206 0.104793714 0.104904073 0.104979456 0.105020009 0.105025835 0.104996872 0.104933144 0.104834824 0.104702098 0.10453514 0.10433423 0.104099645 0.103831591 0.103530268 0.103196039 0.102829239 0.102430269 0.101999415 0.1015371 0.101043718 0.100519859 0.099966037 0.099382724 0.098770453 0.098129814 0.097461401 0.096765824 0.096043648 0.095295489 0.094521667 0.093722743 0.092899168 0.092051434 0.091180007 0.090285603 0.089369042 0.088431015 0.087471827 0.086491523 0.085490637 0.084469538 0.083428878 0.082369134 0.081290657 0.080194098 0.079079967 0.077948838 0.076801221 0.07563781 0.074459474 0.073266929 0.072060853 0.070841882 0.06961069 0.068367811 0.067114127 0.065850626 0.064578103 0.063297275 0.062008828 0.06071344 0.059411909 0.058104959 0.056793414 0.055478126 0.054159689 0.052838674 0.051515884 0.050192201 0.048868532 0.047545509 0.046223748 0.044904078 0.043587189 0.042274048 0.040965595 0.039662732 0.038366378 0.037077107 0.035795621 0.03452296 0.033259852 0.032007103 0.030765494 0.029535843 0.028318917 0.02711556 0.025926542 0.024752484 0.023594027 0.022451727 0.021326027 0.020217455 0.019126531 0.018053825 0.016999952 0.015965327 0.014950239 0.013955074 0.012980092 0.012025305 0.011091182 0.01017807 0.009286251 0.008416244 0.007568324 0.006742544 0.005939034 0.005157898 0.004399306 0.003663209 0.002949466 0.002257843 0.001588295 0.000940844 0.000315325 -0.000288351 -0.000870334 -0.001430846 -0.00197016 -0.002488742 -0.002987122 -0.003465868 -0.003925743 -0.004367422 -0.004791148 -0.005197311 -0.005585976 -0.005957368 -0.006311633 -0.006648885 -0.006969275 -0.007272826 -0.00755965 -0.007829946 -0.008084014 -0.008322498 -0.008545754 -0.008754086 -0.008947826 -0.009127193 -0.009292597 -0.009444414 -0.009583004 -0.009708583 -0.009821543 -0.009922462 -0.01001171 -0.010089565 -0.010156344 -0.010212492 -0.010258358 -0.01029434 -0.010320867 -0.010338341 -0.010347164 -0.010347758 -0.010340369 -0.010325308 -0.010302812 -0.010273058 -0.010236306 -0.010192829 -0.010143002 -0.01008711 -0.010025503 -0.009958502 -0.009886326 -0.009809149 -0.009727064 -0.009640297 -0.009549202 -0.009454095 -0.00935522 -0.009252802 -0.00914711 -0.009038308 -0.00892658 -0.008812078 -0.008694947 -0.008575406 -0.00845366 -0.008329941 -0.008204561 -0.00807767 -0.007949426 -0.00781995 -0.007689369 -0.007557906 -0.007425679 -0.007292765 -0.007159278 -0.007025348 -0.006891177 -0.006756921 -0.006622711 -0.006488701 -0.006355065 -0.006221913 -0.006089349 -0.005957451 -0.005826277 -0.00569586 -0.005566331 -0.005437733 -0.005310147 -0.005183648 -0.005058319 -0.00493425 -0.004811534 -0.004690253 -0.004570392 -0.004452001 -0.004335164 -0.004219971 -0.004106494 -0.003994777 -0.003884849 -0.003776724 -0.003670441 -0.003566049 -0.003463594 -0.0033631 -0.003264597 -0.003168105 -0.003073644 -0.002981174 -0.002890688 -0.002802163 -0.00271561 -0.002631046 -0.002548492 -0.002467899 -0.002389193 -0.002312347 -0.002237352 -0.002164167 -0.002092782 -0.002023185 -0.001955351 -0.001889254 -0.001824827 -0.001762048 -0.001700905 -0.001641371 -0.001583475 -0.001527198 -0.001472474 -0.001419277 -0.001367569 -0.001317309 -0.001268452 -0.001220997 -0.001174933 -0.001130246 -0.00108693 -0.001044982 -0.001004374 -0.000965062 -0.000926975 -0.000890091 -0.000854378 -0.000819824 -0.000786415 -0.000754113 -0.000722864 -0.000692646 -0.000663442 -0.000635239 -0.000608024 -0.000581778 -0.000556465 -0.000532057 -0.000508559 -0.000485965 -0.000464273 -0.000443446 -0.000423455 -0.000404292 -0.000385931 -0.000368352 -0.000351513 -0.000335338 -0.0003198 -0.00030488 -0.000290559 -0.000276816 -0.000263622 -0.000250945 -0.00023877 -0.000227089 -0.000215896 -0.000205174];
        HRFBasisFuncs(:,2)=[-6.16E-07 -5.85E-06 -2.24E-05 -5.75E-05 -0.000117278 -0.000205007 -0.000324387 -0.000479146 -0.000672987 -0.000909615 -0.001192492 -0.00152425 -0.001908093 -0.002346263 -0.002840749 -0.003391801 -0.003999951 -0.0046659 -0.00538953 -0.006171052 -0.007010603 -0.007908905 -0.008866742 -0.009884974 -0.010963759 -0.012102212 -0.013299353 -0.014554769 -0.015867747 -0.017236773 -0.018660477 -0.0201368 -0.021664087 -0.023239698 -0.024860138 -0.026523383 -0.028226344 -0.029964839 -0.031734235 -0.033530489 -0.035350419 -0.037190539 -0.039047121 -0.040916412 -0.042794603 -0.044677851 -0.046562292 -0.048444011 -0.050319111 -0.052183675 -0.054033794 -0.055865571 -0.057675108 -0.059458553 -0.061212082 -0.06293191 -0.064614311 -0.066255589 -0.067852138 -0.069400412 -0.070896945 -0.072338367 -0.073721415 -0.07504293 -0.076300583 -0.077493288 -0.078619325 -0.079676761 -0.080663609 -0.081577478 -0.082417135 -0.083181498 -0.083869327 -0.084479921 -0.085012271 -0.085465613 -0.085839741 -0.086134121 -0.086347995 -0.086481425 -0.086533874 -0.086504728 -0.086393478 -0.08619999 -0.085924684 -0.085567336 -0.085128259 -0.084608065 -0.084007322 -0.083326447 -0.082566553 -0.081728743 -0.080814003 -0.079822898 -0.078756222 -0.077615288 -0.076401239 -0.075114779 -0.073757316 -0.072329713 -0.070833252 -0.069269671 -0.067640746 -0.065948273 -0.064193616 -0.062377923 -0.060502506 -0.05856875 -0.056577673 -0.054530769 -0.052429659 -0.050275863 -0.04807092 -0.04581687 -0.043515534 -0.041168786 -0.038778636 -0.036346477 -0.033873962 -0.031362879 -0.02881539 -0.026233322 -0.023618789 -0.020974301 -0.018302338 -0.01560518 -0.012885005 -0.010144203 -0.007385363 -0.004611043 -0.001823566 0.000974443 0.003780435 0.006591778 0.009405999 0.012220633 0.015033174 0.017840922 0.020641205 0.023431049 0.026207527 0.028967487 0.031707879 0.034426053 0.037119116 0.039784082 0.042417944 0.045018558 0.047583784 0.050111647 0.052599518 0.055044344 0.057443973 0.059796693 0.062099903 0.064350679 0.066546722 0.06868602 0.070766408 0.072786391 0.074744448 0.076637561 0.078463418 0.080220311 0.081907555 0.083524286 0.085070124 0.086543547 0.087943539 0.089269088 0.090518947 0.091692332 0.092788065 0.09380523 0.094743089 0.095602033 0.096381676 0.097081555 0.097701398 0.098241224 0.098701757 0.099084199 0.099389153 0.099615845 0.099764359 0.099834487 0.099827522 0.09974557 0.099589411 0.099359325 0.09905651 0.098683295 0.098241582 0.097732194 0.09715658 0.096516394 0.095813294 0.09504849 0.094223332 0.093340257 0.092401899 0.091409714 0.090365853 0.089272538 0.088131804 0.08694584 0.085717008 0.084447642 0.083140255 0.081797207 0.080420967 0.07901421 0.077579397 0.076118861 0.074634725 0.073129315 0.071604811 0.070063526 0.068507606 0.066939125 0.065360181 0.063772783 0.062178894 0.060580645 0.058979722 0.057378053 0.055777476 0.054179854 0.052586991 0.051000276 0.049421404 0.047852062 0.046293726 0.044747596 0.043215099 0.041697738 0.04019693 0.038713902 0.037249779 0.035805537 0.034381783 0.03297955 0.031599556 0.030242161 0.028906907 0.027593079 0.026301257 0.025031619 0.023785064 0.022562106 0.021363257 0.020189007 0.019039712 0.017916065 0.016818559 0.015747675 0.014703955 0.013687268 0.012697792 0.011735591 0.010800651 0.00989299 0.009012102 0.008157413 0.007328587 0.006525826 0.005748811 0.004996786 0.004269652 0.003567509 0.002890144 0.002236948 0.001607494 0.001001399 0.000418013 -0.000143159 -0.00068283 -0.001201537 -0.001699461 -0.002177283 -0.002635555 -0.003074487 -0.003494235 -0.00389502 -0.004277305 -0.004641275 -0.004987443 -0.005316346 -0.005628199 -0.005923015 -0.006200837 -0.006461939 -0.006706783 -0.006935907 -0.007149862 -0.007349052 -0.007533952 -0.007704848 -0.007861985 -0.008005542 -0.008135754 -0.008253079 -0.008357914 -0.008450712 -0.008532133 -0.008602486 -0.008661919 -0.008710586 -0.008748834 -0.008777175 -0.008795884 -0.008805364 -0.008806124 -0.008798606 -0.008783249 -0.008760394 -0.008730291 -0.00869335 -0.008650202 -0.008601198 -0.008546614 -0.008486672 -0.008421556 -0.008351462 -0.008276719 -0.008197513 -0.008114052 -0.008026543 -0.007935248 -0.007840441 -0.007742431 -0.007641523 -0.007537975 -0.007432055 -0.007323975 -0.00721398 -0.007102243 -0.006988879 -0.006874012 -0.006757737 -0.006640267 -0.006521791 -0.006402471 -0.006282473 -0.006162033 -0.006041376 -0.005920685 -0.005800094 -0.005679704 -0.005559546 -0.005439749 -0.005320481 -0.005201952 -0.005084224 -0.004967221 -0.00485099 -0.004735627 -0.004621227 -0.004507911 -0.004395832 -0.004285074 -0.004175694 -0.004067706 -0.003961166 -0.003856118 -0.003752649 -0.003650864 -0.003550772 -0.003452311 -0.003355503 -0.003260384 -0.003166977 -0.00307527 -0.002985295 -0.002897053 -0.002810517 -0.002725699 -0.002642596 -0.00256121 -0.002481527 -0.002403415 -0.002326852 -0.002251822 -0.002178344 -0.002106444 -0.002036156 -0.001967413 -0.001900167 -0.001834424 -0.001770215 -0.001707555 -0.001646451 -0.001586875 -0.001528828 -0.001472348 -0.001417462 -0.001364198 -0.001312525 -0.001262425 -0.001213908 -0.001166935 -0.001121488 -0.001077518 -0.001034878 -0.000993511 -0.000953371 -0.000914463 -0.000876762 -0.00084023 -0.000804835 -0.000770548 -0.000737355 -0.000705263 -0.000674253];
        HRFBasisFuncs(:,3)=[1.82E-06 1.52E-05 5.47E-05 0.000136622 0.000273526 0.000471283 0.000736563 0.001074855 0.001491855 0.001995458 0.002591106 0.003280847 0.004066564 0.00495131 0.005938251 0.007027233 0.008218646 0.009510582 0.010900576 0.012387891 0.013970199 0.015644605 0.017408639 0.019259984 0.021194135 0.023203793 0.025282401 0.027422997 0.029618813 0.031863747 0.034150005 0.036469365 0.038812692 0.041170703 0.04353589 0.045899308 0.048251947 0.050586782 0.052895955 0.055169882 0.057398168 0.059574056 0.061692205 0.06374739 0.065734468 0.067648404 0.069484291 0.071237323 0.072902869 0.074476433 0.075953693 0.077330501 0.078602879 0.079767062 0.080819484 0.081756792 0.08257586 0.083273782 0.0838479 0.084295802 0.084615325 0.08480459 0.084862172 0.084786793 0.084579428 0.084244849 0.083785943 0.08320438 0.082501708 0.081678465 0.080737034 0.079680363 0.078511311 0.077233588 0.075850088 0.074364637 0.072782281 0.071106707 0.069341191 0.067490149 0.06555648 0.063543136 0.061453236 0.059290703 0.057059753 0.054763582 0.052406247 0.049991957 0.047524408 0.045006521 0.042441696 0.03983399 0.037186678 0.034503262 0.031787278 0.029042423 0.026272512 0.023480668 0.020670424 0.017845234 0.015007606 0.012160911 0.009308749 0.006454302 0.00360053 0.000750839 -0.002092129 -0.00492558 -0.007745877 -0.010549172 -0.013331762 -0.01609012 -0.018820983 -0.021520697 -0.024185165 -0.026810832 -0.029393669 -0.031929544 -0.034414362 -0.036843603 -0.039213676 -0.041520306 -0.043759241 -0.04592691 -0.048019361 -0.050032055 -0.051959266 -0.053796277 -0.055538562 -0.057181475 -0.058720139 -0.060152312 -0.061475092 -0.062683458 -0.063772505 -0.064740207 -0.065584359 -0.066301642 -0.066887815 -0.067339784 -0.067654405 -0.067829579 -0.067863546 -0.06775598 -0.067505589 -0.067111298 -0.066572649 -0.065892429 -0.065073595 -0.064120089 -0.063032015 -0.061809204 -0.060454256 -0.058970463 -0.057358823 -0.055620747 -0.05375991 -0.051779852 -0.049681623 -0.047468186 -0.045144117 -0.042709808 -0.040168497 -0.037524218 -0.034782529 -0.031949894 -0.029033663 -0.026039466 -0.022973608 -0.01984235 -0.016651067 -0.013405458 -0.010111243 -0.006774272 -0.00340086 1.32E-06 0.003425863 0.006866322 0.010315953 0.013768558 0.017217781 0.020658001 0.02408341 0.027487527 0.030864353 0.034207118 0.037509785 0.040766829 0.043972333 0.047120476 0.050205985 0.053225257 0.056176025 0.05905357 0.061853455 0.064572466 0.067207077 0.069752181 0.072203888 0.074561845 0.076824971 0.078990497 0.081056696 0.083020995 0.084881447 0.086636887 0.088285882 0.089825765 0.091254794 0.092572673 0.093779568 0.094874962 0.095859771 0.096734867 0.09750083 0.098157868 0.098705572 0.099144438 0.099477256 0.099705193 0.099830539 0.099857021 0.099785363 0.099616742 0.099353127 0.098995253 0.09854476 0.098004404 0.097377243 0.096666665 0.095875471 0.095006714 0.094063854 0.093051327 0.091972035 0.090828721 0.089624738 0.088363592 0.087048889 0.085684416 0.08427401 0.082822074 0.081332989 0.079811457 0.078261696 0.076687724 0.075092375 0.07347826 0.071847305 0.070202479 0.068546969 0.066883292 0.065213201 0.063538619 0.061861183 0.060182527 0.058504809 0.056830495 0.05516166 0.053500007 0.051847449 0.050205517 0.048575509 0.046958413 0.045355078 0.043766677 0.042193995 0.04063801 0.039100148 0.037581457 0.03608289 0.034604981 0.033148428 0.031714053 0.03030222 0.028913392 0.027548023 0.026206405 0.024889181 0.023596314 0.022328017 0.021084688 0.019867012 0.018675593 0.017510319 0.016371534 0.015259262 0.014173421 0.013114555 0.012083194 0.011079822 0.010104486 0.009156973 0.008236852 0.007343737 0.006477419 0.005637581 0.004824097 0.004037016 0.003276544 0.002542638 0.001835079 0.001153942 0.000499209 -0.000129847 -0.000733376 -0.001311269 -0.001863413 -0.002390181 -0.002892508 -0.003370819 -0.003825777 -0.004258255 -0.004668946 -0.0050585 -0.005427223 -0.005775204 -0.006102965 -0.006411499 -0.006701325 -0.006972759 -0.007226197 -0.007461977 -0.007680241 -0.007881534 -0.008066471 -0.00823516 -0.008388037 -0.00852582 -0.008648734 -0.008757374 -0.00885227 -0.008934052 -0.00900341 -0.009060785 -0.009106613 -0.009141101 -0.009164341 -0.009176506 -0.009177977 -0.009169326 -0.009150858 -0.009122882 -0.00908603 -0.009041171 -0.00898917 -0.008930316 -0.008865013 -0.00879362 -0.008716483 -0.008633814 -0.008545976 -0.008453522 -0.008356776 -0.008255772 -0.008150789 -0.008042292 -0.007930768 -0.007816569 -0.007699948 -0.007581132 -0.007460294 -0.007337591 -0.007213183 -0.007087259 -0.00696015 -0.006832074 -0.006703167 -0.006573488 -0.0064433 -0.006312865 -0.006182357 -0.006051962 -0.005921941 -0.005792464 -0.005663567 -0.005535316 -0.005407808 -0.005281287 -0.005155914 -0.005031506 -0.00490804 -0.004785616 -0.004664382 -0.004544464 -0.004425976 -0.00430878 -0.004192796 -0.004078194 -0.003965096 -0.003853552 -0.003743606 -0.003635224 -0.003528405 -0.00342322 -0.003319775 -0.003218198 -0.003118495 -0.003020697 -0.002924864 -0.002831009 -0.002739188 -0.00264942 -0.002561438 -0.002475164 -0.002390532 -0.00230759 -0.00222631 -0.002146707 -0.002068791 -0.00199258 -0.001918111 -0.001845426 -0.001774544];
        HRFBasisFuncs(:,4)=[-2.59E-06 -2.15E-05 -7.77E-05 -0.000194547 -0.000389572 -0.000672975 -0.001055539 -0.001542332 -0.002138632 -0.002853077 -0.003692304 -0.004659239 -0.005758562 -0.006992551 -0.00836293 -0.009867055 -0.01150335 -0.013269035 -0.015160338 -0.017169972 -0.019291107 -0.021514884 -0.023828125 -0.026219399 -0.028681444 -0.03120815 -0.033789922 -0.036414235 -0.039068306 -0.041740446 -0.044419152 -0.047091252 -0.049741186 -0.052354149 -0.054915122 -0.057406303 -0.059813095 -0.062125211 -0.064330845 -0.066418486 -0.068374122 -0.070189845 -0.071861253 -0.073384236 -0.074754953 -0.075969851 -0.077025679 -0.077919475 -0.07864861 -0.079210764 -0.079603957 -0.079826534 -0.079877184 -0.079754944 -0.079459198 -0.078989683 -0.078346491 -0.07753008 -0.076541257 -0.075381194 -0.074051421 -0.072553861 -0.070891179 -0.069066104 -0.067084255 -0.06495664 -0.062691702 -0.060296419 -0.057777263 -0.055138983 -0.052390664 -0.049541624 -0.04660016 -0.043576631 -0.040480062 -0.037319376 -0.034103362 -0.030839672 -0.027534683 -0.024196341 -0.020832053 -0.017448566 -0.014052713 -0.01065146 -0.007252836 -0.003863491 -0.000489842 0.002862851 0.006189257 0.009483264 0.012739228 0.015952313 0.019117851 0.022230407 0.025286082 0.028281151 0.031210221 0.034068581 0.036852226 0.039556174 0.042175378 0.044706297 0.047146875 0.049494282 0.051745913 0.053898243 0.055946436 0.057885087 0.059707919 0.061410976 0.062990957 0.064443511 0.065764917 0.06695257 0.068003352 0.068913354 0.069678784 0.070293431 0.070752412 0.071051755 0.071190012 0.071164434 0.070973464 0.070616231 0.07009174 0.069398022 0.068531076 0.067489454 0.066273376 0.064883877 0.063320199 0.061587235 0.059689332 0.057627155 0.055400946 0.053017801 0.050486053 0.047810806 0.044995478 0.04204525 0.038965696 0.035763738 0.032447787 0.029028772 0.025515772 0.021916361 0.01824 0.014498015 0.01070179 0.006864208 0.002996719 -0.000891011 -0.004786558 -0.008677379 -0.012551616 -0.016399676 -0.02021189 -0.023977952 -0.027684046 -0.031312478 -0.03485147 -0.038286643 -0.041607358 -0.044799767 -0.047851871 -0.050752186 -0.053493077 -0.056065923 -0.058461511 -0.060674089 -0.062696508 -0.06452404 -0.066149255 -0.067564019 -0.068761334 -0.069738843 -0.070493512 -0.071022164 -0.071325264 -0.071404688 -0.07126197 -0.070899283 -0.070320961 -0.069524063 -0.068508279 -0.067272364 -0.065822171 -0.064168575 -0.062317888 -0.060274989 -0.058046278 -0.055643988 -0.053081209 -0.050365073 -0.047504141 -0.044507698 -0.041384786 -0.038142069 -0.034786644 -0.031330365 -0.02778559 -0.024160805 -0.020467759 -0.016716942 -0.012918261 -0.009081895 -0.005216181 -0.001328688 0.002572033 0.006476676 0.010374782 0.01425708 0.01811465 0.021939004 0.025722293 0.029455862 0.033132768 0.036746831 0.040291621 0.043761274 0.047148516 0.050447082 0.053651088 0.056753553 0.059748836 0.062632251 0.06539964 0.06804668 0.070569576 0.072966545 0.075234901 0.07737303 0.079380454 0.081256754 0.082999425 0.08460753 0.086081231 0.08741977 0.088623931 0.089694846 0.090634918 0.091447121 0.092135364 0.092705958 0.093167093 0.093528349 0.093792326 0.093963771 0.094043759 0.094034671 0.093938404 0.093756423 0.093491701 0.093145171 0.092717806 0.092211006 0.091627561 0.09097323 0.09025035 0.089461488 0.088609017 0.087695034 0.086724631 0.085701087 0.084627353 0.083505138 0.082338357 0.081130703 0.079884444 0.078601106 0.077283663 0.075935898 0.074560599 0.073160326 0.07173863 0.070298263 0.068842713 0.067374732 0.065896286 0.064410003 0.062918294 0.061423263 0.059927043 0.058430997 0.056936415 0.055444821 0.053957668 0.052475868 0.051000259 0.049531977 0.048072767 0.046624491 0.045188358 0.043765382 0.042356613 0.040962508 0.039583429 0.038220499 0.036874719 0.035546988 0.034238342 0.032949866 0.031682247 0.030435652 0.029209946 0.02800603 0.026824714 0.025667111 0.024533959 0.023424506 0.022338755 0.021277041 0.020239164 0.019225396 0.018235663 0.017270498 0.016330704 0.015415875 0.014524926 0.013658044 0.012815605 0.011997859 0.011204949 0.010436756 0.009692963 0.008973131 0.008277333 0.007605455 0.006956819 0.006331551 0.00572899 0.005148606 0.004589869 0.004051875 0.003534422 0.003037363 0.002560935 0.002105359 0.001670421 0.001255196 0.000859084 0.000482059 0.000124081 -0.00021578 -0.000538996 -0.000846826 -0.001139436 -0.001417199 -0.001680211 -0.001928847 -0.0021634 -0.002384269 -0.002592162 -0.002787537 -0.002970898 -0.00314246 -0.003302559 -0.003451393 -0.003589148 -0.003715955 -0.003831947 -0.003937411 -0.004032949 -0.004118775 -0.004195126 -0.004262325 -0.004320498 -0.004369892 -0.004410921 -0.004443884 -0.004469162 -0.004487304 -0.004498828 -0.004504156 -0.004503633 -0.004497407 -0.004485728 -0.004468944 -0.004447603 -0.004422025 -0.004392505 -0.004359283 -0.004322646 -0.004282724 -0.004239705 -0.004193836 -0.004145167 -0.004093664 -0.004039569 -0.003983085 -0.003924315 -0.003863481 -0.003800793 -0.003736405 -0.003670425 -0.003603025 -0.003534388 -0.00346462 -0.00339381 -0.003322118 -0.003249808 -0.003177089 -0.003104057 -0.003030514 -0.002956581 -0.002882276 -0.002807794 -0.002733187 -0.002658609 -0.002584177 -0.002509998 -0.002436185 -0.002362811 -0.00229001];
        HRFBasisFuncs(:,5)=[3.19E-06 2.87E-05 0.000107287 0.000275361 0.000552341 0.000949557 0.001481071 0.002154624 0.002979557 0.003958665 0.005091261 0.006375626 0.007809018 0.009387893 0.01110715 0.01296653 0.014964118 0.017092884 0.019345872 0.021715386 0.024192642 0.026769702 0.029434332 0.032165104 0.034945693 0.03775974 0.040591583 0.043421391 0.046228221 0.048993263 0.051699261 0.054327996 0.0568605 0.059285265 0.061591525 0.063755966 0.065757487 0.067584713 0.069233196 0.070699787 0.071964258 0.073017602 0.073858355 0.074485525 0.07489858 0.075097455 0.075082547 0.07485472 0.074415308 0.0737661 0.072909357 0.071847778 0.070584538 0.069123244 0.067467937 0.065623107 0.063593636 0.061384858 0.059002488 0.056452626 0.053741768 0.050876792 0.047865382 0.044715271 0.041438832 0.038056237 0.034582946 0.031033337 0.027419394 0.023750316 0.020039268 0.016299278 0.012544034 0.0087903 0.005050108 0.001334303 -0.002345368 -0.005976897 -0.009548396 -0.013047038 -0.016464196 -0.01979081 -0.023018424 -0.026138996 -0.029144282 -0.032025789 -0.034776438 -0.037389119 -0.039857902 -0.042178045 -0.044346349 -0.046359486 -0.04821364 -0.049904622 -0.051430001 -0.05278943 -0.053979569 -0.054995306 -0.055835216 -0.056496861 -0.056981427 -0.057290295 -0.05742466 -0.057385804 -0.057171809 -0.056780485 -0.056210286 -0.055462538 -0.05453486 -0.053428474 -0.052144914 -0.050685465 -0.049054756 -0.047257825 -0.045296783 -0.043176054 -0.040900844 -0.038474013 -0.035899015 -0.033180091 -0.030326009 -0.027345748 -0.024247289 -0.021038874 -0.017727934 -0.01432172 -0.010827982 -0.007255812 -0.003614316 8.58E-05 0.003831941 0.007611257 0.011412238 0.015221811 0.019024489 0.022804837 0.026547428 0.030239049 0.033862762 0.037406007 0.040855163 0.044200366 0.047431611 0.050536898 0.053503221 0.056318345 0.058970983 0.061455268 0.063765358 0.065894508 0.067827247 0.06955087 0.071057402 0.072340565 0.073387056 0.074191548 0.074752565 0.07506938 0.075135648 0.074940093 0.074483457 0.07376045 0.07277345 0.071525006 0.070022578 0.068271384 0.066284094 0.064066953 0.06162764 0.05897641 0.0561219 0.053078233 0.049855829 0.04646749 0.042925269 0.039248301 0.035451208 0.031546142 0.02754807 0.023471502 0.019332734 0.0151507 0.010942944 0.006720855 0.002495563 -0.001717716 -0.005898805 -0.010030811 -0.014097574 -0.018083012 -0.021972741 -0.025750716 -0.029401979 -0.032909426 -0.03626029 -0.039443589 -0.042449944 -0.045267645 -0.047885938 -0.050300317 -0.052503171 -0.054487609 -0.056251233 -0.057790419 -0.059101429 -0.060179978 -0.061023973 -0.061630352 -0.061995566 -0.062118059 -0.061998337 -0.061636964 -0.061039765 -0.060214169 -0.05916625 -0.057900027 -0.056418452 -0.054729383 -0.052845305 -0.050774191 -0.048527605 -0.046124677 -0.043572806 -0.040879545 -0.03805586 -0.035106017 -0.032039588 -0.028868438 -0.025603238 -0.022254761 -0.01883214 -0.01534574 -0.011808442 -0.008235275 -0.004635786 -0.001020038 0.002601251 0.006218972 0.009824988 0.013412425 0.016969796 0.020486997 0.023954158 0.027362896 0.030704056 0.033970928 0.037158506 0.040261327 0.043275909 0.046198612 0.0490271 0.051758644 0.054391235 0.056922549 0.059349953 0.061670801 0.063884759 0.065991785 0.067990571 0.069880387 0.071659883 0.073328667 0.074891835 0.076352981 0.077714473 0.078976096 0.080140829 0.081211334 0.082188354 0.083072643 0.08386667 0.084575876 0.085204715 0.085755074 0.086231105 0.086634779 0.086969239 0.087238308 0.0874434 0.087585895 0.087668483 0.08769327 0.087662311 0.087575947 0.087435894 0.087243835 0.087002033 0.086711002 0.086371051 0.085982627 0.085546797 0.085065835 0.084541579 0.083976323 0.083371952 0.082729067 0.082048021 0.081329497 0.080574705 0.079784927 0.078961525 0.078106009 0.077220556 0.076306877 0.075366769 0.0744018 0.073412957 0.072401354 0.071367609 0.070312313 0.069237077 0.06814427 0.067036126 0.065914358 0.06478024 0.06363503 0.062480106 0.061316554 0.06014532 0.058966996 0.057782044 0.056592153 0.055398703 0.054203548 0.053008403 0.051814628 0.050623849 0.049437659 0.048257083 0.047083065 0.045916112 0.044756995 0.043606823 0.042466543 0.041336809 0.040217931 0.039110593 0.038015578 0.036933925 0.035866389 0.034813398 0.033776058 0.032755207 0.031751059 0.030763083 0.029790678 0.028834278 0.027893936 0.026970427 0.02606393 0.025174971 0.024304445 0.023452803 0.022619447 0.021803234 0.021004229 0.020222712 0.019458629 0.018712288 0.017983876 0.017273277 0.016580315 0.015904526 0.015246075 0.014605202 0.013981839 0.013376396 0.012788826 0.012218501 0.011665195 0.011128017 0.010606245 0.010099193 0.009606971 0.009129414 0.008666509 0.008218239 0.00778452 0.007365275 0.006960161 0.006568394 0.006189774 0.005824133 0.005471355 0.005131082 0.004803041 0.004487133 0.004183456 0.003892095 0.003613126 0.003346634 0.003092388 0.002849745 0.002618265 0.00239802 0.002188984 0.001991184 0.001804316 0.001628067 0.001462372 0.001306848 0.001161175 0.001024863 0.000897218 0.000777776 0.000666339 0.000562565 0.000466302 0.000377229 0.000294941 0.000219021 0.000149123 8.51E-05 2.65E-05];
        
        
        disp('Interpolating design matrix to be at 1000Hz sampling rate (this is thought to improve convolution accuracy)')
        samplingRateInHz=1/TR_INSECONDS;
        time = (1:size(taskdesignmat,1))*1/samplingRateInHz;
        newtime = time(1):1/1000:time(end); % new sampling
        new_designmat = interp1(time,taskdesignmat,newtime,'linear');
        
        disp('Interpolating HRFs to be at 1000Hz sampling rate')
        samplingRateInHz=samplingRateOfHRF;
        time = (1:size(HRFBasisFuncs,1))*1/samplingRateInHz;
        newtime = time(1):1/1000:time(end); % new sampling
        new_HRFBasisFuncs = interp1(time,HRFBasisFuncs,newtime,'linear');
        
        disp('Convolving with hemodynamic response function (HRF) basis set')
        taskdesignmat_hrf1=zeros(size(new_designmat,1),size(new_designmat,2)*numBasisFuncs);
        regressorCount=1;
        for basisNum = 1:numBasisFuncs
            thisBasisHRF=new_HRFBasisFuncs(:,basisNum);
            for regressorNum_ForThisBasis=1:size(new_designmat,2)
                convData=conv(new_designmat(:,regressorNum_ForThisBasis),thisBasisHRF);
                taskdesignmat_hrf1(:,regressorCount)=convData(1:size(new_designmat,1),:);
                regressorCount=regressorCount+1;
            end
        end
        
        disp('Downsampling designmatrix back to data sampling rate')
        taskdesignmat_hrf=downsample(taskdesignmat_hrf1,TR_INSECONDS*1000);
        HRFBasisFuncs_inTRs=downsample(new_HRFBasisFuncs,TR_INSECONDS*1000);
        
        taskTiming.taskdesignmat=taskdesignmat;
        taskTiming.taskdesignmat_hrf=taskdesignmat_hrf;
        %taskTiming.EDATIMPORT=EDATIMPORT;
        %taskTiming.reffunc=reffunc;
        taskTiming.HRFBasisFuncs_inTRs=HRFBasisFuncs_inTRs;
        

        
        %Save task timing variables to file (for more efficient processing in future when EXECUTE=0)
        disp(['Saving results to: ' savedfile])
        save(savedfile, 'taskTiming')
        
        
        
        
        
    else
        
        if ~isempty(TASKRUNS)
            %Load task timing variables from file (for more efficient processing when EXECUTE=0)
            disp(['Loading results from: ' savedfile])
            load(savedfile);
        end
        
    end
    
    if ~isempty(TASKRUNS)
        output.taskTiming{subjIndex}=taskTiming;
    end
    
    
    %% Rest nuisance regression
    %Set to 1 if you want to run this procedure again for subjects that already had it run before (otherwise it will load from previously saved files)
    RERUN_RESTGLM=0;
    
    savedfile=[subjTemporaryAnalysisDir mfilenameRep '_' ANALYSISNAME '_RestNuisanceGLMVars.mat'];
    
    %Only execute this if RERUN_TASKGLM==1 or the savedfile doesn't exist for this subject, but skip if no TASKRUNS (i.e., only rest data are included)
    if and(~isempty(RESTRUNS), or(RERUN_RESTGLM == 1, ~exist(savedfile, 'file')))
        disp('Running rest nuisance regression')
        
        %Specify the number of nuisance regressors
        NUMREGRESSORS_NUISANCE=16;
        %Add 2 regressors for GSR
        if GSR
            NUMREGRESSORS_NUISANCE=NUMREGRESSORS_NUISANCE+2;
        end
        visualizeDesignMatrix=0;
        
        restGLMVars = runGLM(tseriesMatSubj, NUMREGRESSORS_NUISANCE, nuisanceTSVars, [], RUNLENGTHS, RESTRUNS, NPROC, NUMPARCELS, GSR, visualizeDesignMatrix);
        
        
       %Save task GLM variables to file (for more efficient processing in future when EXECUTE=0)
        disp(['Saving results to: ' savedfile])
        save(savedfile, 'restGLMVars')
        
    else
        
        if ~isempty(RESTRUNS)
            %Load task GLM variables from file (for more efficient processing when EXECUTE=0)
            disp(['Loading results from: ' savedfile])
            load(savedfile);
        end
        
    end
    
    if ~isempty(RESTRUNS)
        output.rest_fMRI_preprocTS{subjIndex} = restGLMVars.fMRI_resids;
        output.restGLMVars{subjIndex} = restGLMVars;
    end
    
    
    %% Task nuisance regression and GLM
    
    %Set to 1 if you want to run this procedure again for subjects that already had it run before (otherwise it will load from previously saved files)
    RERUN_TASKGLM=1;
    
    savedfile=[subjTemporaryAnalysisDir mfilename '_' ANALYSISNAME '_TaskGLM.mat'];
    
    %Only execute this if RERUN_TASKGLM==1 or the savedfile doesn't exist for this subject, but skip if no TASKRUNS (i.e., only rest data are included)
    if and(~isempty(TASKRUNS), or(RERUN_TASKGLM == 1, ~exist(savedfile, 'file')))

        disp('Running task nuisance regression and GLM')
        
        %Specify the number of nuisance regressors
        NUMREGRESSORS_NUISANCE=16;
        %Add 2 regressors for GSR
        if GSR
            NUMREGRESSORS_NUISANCE=NUMREGRESSORS_NUISANCE+2;
        end
        visualizeDesignMatrix=0;
        
        taskGLMVars = runGLM(tseriesMatSubj, NUMREGRESSORS_NUISANCE, nuisanceTSVars, taskTiming.taskdesignmat_hrf, RUNLENGTHS, TASKRUNS, NPROC, NUMPARCELS, GSR, visualizeDesignMatrix);
        
        %Save task GLM variables to file (for more efficient processing in future when EXECUTE=0)
        disp(['Saving results to: ' savedfile])
        save(savedfile, 'taskGLMVars')
        
    else
        
        if ~isempty(TASKRUNS)
            %Load task GLM variables from file (for more efficient processing when EXECUTE=0)
            disp(['Loading results from: ' savedfile])
            load(savedfile);
        end
        
    end
    
    if ~isempty(TASKRUNS)
        output.task_fMRI_preprocTS{subjIndex} = taskGLMVars.fMRI_resids;
        output.task_betas{subjIndex} = taskGLMVars.fMRI_betas;
        output.taskGLMVars{subjIndex} = taskGLMVars;
    end
    
    
    
    
    %% Apply temporal filtering
    
    %Set to 1 if you want to run this procedure again for subjects that already had it run before (otherwise it will load from previously saved files)
    RERUN_TEMPORALFILTER=0;
    
    savedfile=[subjTemporaryAnalysisDir mfilenameRep '_' ANALYSISNAME '_TemporalFilter.mat'];
    
    if TEMPORALFILTER==1
        
        %Only execute this if RERUN_TEMPORALFILTER==1 or the savedfile doesn't exist for this subject, but skip if TEMPORALFILTER==0
        if or(RERUN_TEMPORALFILTER == 1, ~exist(savedfile, 'file'))
            
            disp('Applying temporal filter')
            
            %Create temporal filter (based on Jonathan Power's script from Petersen lab)
            lopasscutoff=LOWPASS_HZ/(0.5/TR_INSECONDS); % since TRs vary have to recalc each time
            hipasscutoff=HIGHPASS_HZ/(0.5/TR_INSECONDS); % since TRs vary have to recalc each time
            %Using filter order of 1
            filtorder=1;
            [butta, buttb]=butter(filtorder,[hipasscutoff lopasscutoff]);
            
            %Rest data
            if ~isempty(RESTRUNS)
                
                %Interpolate data if scrubbing
                %Use interpolation to account for gaps in time series due to motion scrubbing
                if IMPLEMENT_MOTIONSCRUBBING==1
                    restfMRIData_scrubbed=restGLMVars.fMRI_resids;
                    fMRIData_rest = interpolateAcrossTSGaps(restfMRIData_scrubbed, restGLMVars.temporalmask, RUNLENGTHS, RESTRUNS, NUMPARCELS);
                else
                    fMRIData_rest=restGLMVars.fMRI_resids;
                end
                
                %Apply temporal filter
                filteredData=filtfilt(butta,buttb,fMRIData_rest');
                filteredData=filteredData';
                
                %Reapply scrubbing
                if IMPLEMENT_MOTIONSCRUBBING==1
                    filteredData=filteredData(:,logical(restGLMVars.temporalmask));
                end
                
                filteredDataOutput.rest_fMRI_preprocTS=filteredData;
                
            end
            
            %Task data
            if ~isempty(TASKRUNS)
                
                %Interpolate data if scrubbing
                %Use interpolation to account for gaps in time series due to motion scrubbing
                if IMPLEMENT_MOTIONSCRUBBING==1
                    taskfMRIData_scrubbed=taskGLMVars.fMRI_resids;
                    fMRIData_task = interpolateAcrossTSGaps(taskfMRIData_scrubbed, taskGLMVars.temporalmask, RUNLENGTHS, TASKRUNS, NUMPARCELS);
                else
                    fMRIData_task=taskGLMVars.fMRI_resids;
                end
                
                %Apply temporal filter
                filteredData=filtfilt(butta,buttb,fMRIData_task');
                filteredData=filteredData';
                
                %Reapply scrubbing
                if IMPLEMENT_MOTIONSCRUBBING==1
                    filteredData=filteredData(:,logical(taskGLMVars.temporalmask));
                end
                
                filteredDataOutput.task_fMRI_preprocTS=filteredData;
                
            end
            
            %Save filtered time series variables to file (for more efficient processing in future when EXECUTE=0)
            disp(['Saving results to: ' savedfile])
            save(savedfile, 'filteredDataOutput')
            
        else
            
            %Load filtered time series from file (for more efficient processing when EXECUTE=0)
            disp(['Loading results from: ' savedfile])
            load(savedfile);
            
        end
        
        if ~isempty(TASKRUNS)
            output.task_fMRI_preprocTS{subjIndex}=filteredDataOutput.task_fMRI_preprocTS;
        end
        if ~isempty(RESTRUNS)
            output.rest_fMRI_preprocTS{subjIndex}=filteredDataOutput.rest_fMRI_preprocTS;
        end
        
    end
end


%% Save final result output

outputfilename=[outputdatadir 'results/' mfilename '_' ANALYSISNAME '_output.mat'];
disp(['Saving final results to: ' outputfilename])
save(outputfilename, 'output', '-v7.3');

disp('==Be sure to delete intermediate temporary files when finished with preprocessing all subjects'' data==')
disp('Shell command to delete intermediate temporary files:')
subjTemporaryAnalysisDir_mod=strrep(subjTemporaryAnalysisDir, num2str(SUBJECTLIST{numSubjs}),'*');
disp(['rm -rfv ' subjTemporaryAnalysisDir_mod])

end


%%%%%%%%%%%%%%%%%%%%%%%%%%

function [output_GLM] = runGLM(tseriesMatSubj,NUMREGRESSORS_NUISANCE,nuisanceTSVars,taskdesignmat_hrf,RUNLENGTHS,runNums,NPROC,NUMPARCELS,GSR, visualizeDesignMatrix)

%Specify the number of task regressors
numregressors_task=size(taskdesignmat_hrf,2);
NUMREGRESSORS=NUMREGRESSORS_NUISANCE+numregressors_task;
%Add 1 regressor for each run (to account for linear trend within each run)
if NUMREGRESSORS_NUISANCE>0
    numregressors_extra=length(runNums);
else
    numregressors_extra=0;
end

%Concatenate runs, Organize nuisance regressors
tseriesMatSubj_fMRIconcat=zeros(NUMPARCELS,sum(RUNLENGTHS(runNums)));
X=zeros(NUMREGRESSORS+numregressors_extra,sum(RUNLENGTHS(runNums)));
tmask=ones(sum(RUNLENGTHS(runNums)),1);
runTimings=zeros(length(runNums),sum(RUNLENGTHS(runNums)));
for taskRunIndex=1:length(runNums)
    if taskRunIndex>1
        priorRunsLength=sum(RUNLENGTHS(runNums(1):runNums(taskRunIndex-1)));
    else
        priorRunsLength=0;
    end
    thisRunLength=RUNLENGTHS(runNums(taskRunIndex));
    runStart=priorRunsLength+1;
    runEnd=priorRunsLength+thisRunLength;
    if NUMREGRESSORS_NUISANCE>0  %Assuming nuisance regression already run if NUMREGRESSORS_NUISANCE == 0, fMRI data already concatenated
        %fMRI data
        tseriesMatSubj_fMRIconcat(:,runStart:runEnd)=tseriesMatSubj(:,1:thisRunLength,runNums(taskRunIndex));
        
        %White matter nuisance regressors
        X(1,runStart:runEnd)=nuisanceTSVars.nuisanceTS_whitematter(1:thisRunLength,runNums(taskRunIndex));
        X(2,runStart:runEnd)=[0; diff(nuisanceTSVars.nuisanceTS_whitematter(1:thisRunLength,runNums(taskRunIndex)))];
        %Ventricle
        X(3,runStart:runEnd)=nuisanceTSVars.nuisanceTS_ventricles(1:thisRunLength,runNums(taskRunIndex));
        X(4,runStart:runEnd)=[0; diff(nuisanceTSVars.nuisanceTS_ventricles(1:thisRunLength,runNums(taskRunIndex)))];
        %Motion (12 regressors)
        X(5:16,runStart:runEnd)=nuisanceTSVars.nuisanceTS_motion(:,1:thisRunLength,runNums(taskRunIndex));
        %Run global signal regression if specified
        if GSR
            X(17,runStart:runEnd)=nuisanceTSVars.nuisanceTS_wholebrain(1:thisRunLength,runNums(taskRunIndex));
            X(18,runStart:runEnd)=[0; diff(nuisanceTSVars.nuisanceTS_wholebrain(1:thisRunLength,runNums(taskRunIndex)))];
        end
        %Linear trend for run regressor
        X(NUMREGRESSORS_NUISANCE+taskRunIndex,runStart:runEnd)=linspace(0,1,thisRunLength);
    end
    %Add task regressors
    if numregressors_task>0
        X((NUMREGRESSORS_NUISANCE+numregressors_extra+1):(NUMREGRESSORS+numregressors_extra),runStart:runEnd)=taskdesignmat_hrf(runStart:runEnd,:)';
    end
    %Run transition regressor
    %X(NUMREGRESSORS+taskRunIndex,runStart:runEnd)=ones(thisRunLength,1);
    runTimings(taskRunIndex,runStart:runEnd)=ones(thisRunLength,1);
    %Temporal mask
    tmask(runStart:runEnd)=nuisanceTSVars.temporalMask(1:thisRunLength,runNums(taskRunIndex));
end

%Zscore the design matrix to make it easier to visualize
if visualizeDesignMatrix
    disp('==Make sure to check over the design matrix visually==')
    Xzscored=zscore(X,0,2);
    Xzscored(logical(eye(size(Xzscored))))=0;
    figure;imagesc(Xzscored);title('Regressors');
    disp('Also showing rank correlation among regressors')
    rankMat=zeros(size(X));
    for ind=1:size(X,1)
        rankMat(ind,:)=tiedrank(X(ind,:));
    end
    rankCorrMat=corrcoef(rankMat');
    rankCorrMat(logical(eye(size(rankCorrMat))))=0;
    figure;imagesc(rankCorrMat);title('Regressor Spearman correlations');
end

%Apply temporal mask
X_orig=X;
X_tmasked=X(:,logical(tmask));
X=X_tmasked;
if NUMREGRESSORS_NUISANCE>0
    tseriesMatSubj_fMRIconcat=tseriesMatSubj_fMRIconcat(:,logical(tmask));
else
    tseriesMatSubj_fMRIconcat=tseriesMatSubj;
end

%Test rank of design matrix
matrixRank=rank(X);
disp(['Number of regressors: ' num2str(size(X,1))])
disp(['Rank of matrix: ' num2str(matrixRank)])
if matrixRank < size(X,1)
    disp('ERROR: Matrix is rank deficient; fix the design matrix before running your regression.')
    disp('Consider using PCA on the nuisance regressors (to orthogonalize them). (Do not orthogonalize the task regressors)')
end

% Instantiate empty arrays
fMRI_resids = zeros(size(tseriesMatSubj_fMRIconcat));
fMRI_betas = zeros(NUMPARCELS, size(X,1));

% Begin for loop
parpool(NPROC);
parfor (regionNum=1:NUMPARCELS, NPROC)
%for (regionNum=1:NUMPARCELS)
    % Get the region's data
    ROITimeseries = tseriesMatSubj_fMRIconcat(regionNum,:);
    
    %Remove mean of time series for each run (do not need a constant term in this case)
    lastRunEnd=0;
    for taskRunIndex=1:length(runNums)
        thisRunEnd=lastRunEnd+sum(tmask(logical(runTimings(taskRunIndex,:))));
        tmaskForRun=(lastRunEnd+1):thisRunEnd;
        ROITimeseries(tmaskForRun)=ROITimeseries(tmaskForRun)-mean(ROITimeseries(tmaskForRun));
        lastRunEnd=thisRunEnd;
    end
    
    % Regress out the nuisance time series, keep the residuals and betas
    [beta,bint,resid] = regress(ROITimeseries', X');
    
    % Collect rest regression results
    fMRI_resids(regionNum, :) = resid;
    fMRI_betas(regionNum,:) = beta';

end
delete(gcp);

output_GLM.fMRI_resids = fMRI_resids;
output_GLM.fMRI_betas = fMRI_betas;
output_GLM.temporalmask = tmask;
output_GLM.X_orig = X_orig;
output_GLM.X_tmasked = X_tmasked;
%Use this for task functional connectivity analyses (e.g., choose time points using a threshold of 0.5)
if numregressors_task>0
    output_GLM.taskdesignmat_hrf_tmasked = taskdesignmat_hrf(logical(tmask),:);
end

end


%%%%%%%%%%%%%%%%%%%

function [output_InterpolatedData] = interpolateAcrossTSGaps(fMRI_tseries, tmask, RUNLENGTHS, runNums, NUMPARCELS)

%Place fMRI data into original-sized matrix, with NaNs in the gaps
fMRIData_scrubbed=fMRI_tseries;
fMRIData_withNans=nan(NUMPARCELS,sum(RUNLENGTHS(runNums)));
fMRIData_withNans(:,logical(tmask))=fMRIData_scrubbed;

fMRIData_interpolated=zeros(NUMPARCELS,sum(RUNLENGTHS(runNums)));

for taskRunIndex = 1:length(runNums)
    
    %Prep run timing and data
    if taskRunIndex>1
        priorRunsLength=sum(RUNLENGTHS(runNums(1):runNums(taskRunIndex-1)));
    else
        priorRunsLength=0;
    end
    thisRunLength=RUNLENGTHS(runNums(taskRunIndex));
    runStart=priorRunsLength+1;
    runEnd=priorRunsLength+thisRunLength;
    %fMRI data
    fMRIdata_thisrun=fMRIData_withNans(:,runStart:runEnd);
    
    %Identify gaps in time series due to motion scrubbing
    bd=isnan(fMRIdata_thisrun);
    nongap_timepoints=find(~bd);
    bd([1:(min(nongap_timepoints)-1) (max(nongap_timepoints)+1):end])=0;
    %Implement linear interpolation across gaps (but not at beginning and end of run)
    fMRIdata_thisrun_interpolated=fMRIdata_thisrun;
    fMRIdata_thisrun_interpolated(bd)=interp1(nongap_timepoints,fMRIdata_thisrun(nongap_timepoints),find(bd));
    
    %Set NaNs to 0
    fMRIdata_thisrun_interpolated(isnan(fMRIdata_thisrun_interpolated))=0;
    
    %Output result
    fMRIData_interpolated(:,runStart:runEnd)=fMRIdata_thisrun_interpolated;
    
end

output_InterpolatedData=fMRIData_interpolated;

end



