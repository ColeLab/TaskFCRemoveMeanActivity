
%% Set parameters
%Subject numbers:
subjNums=[100307 120212 162329 196144 221319 530635 702133 885975 123117 142828 163432 197550 249947 552544 729557 887373 103414 124422 143325 167743 199150 250427 559053 732243 889579 103515 125525 255639 579665 734045 894673 103818 128632 149337 172332 200614 293748 581349 748258 896778 105115 129028 149539 175439 201111 298051 585862 753251 896879 110411 130013 150423 205119 304020 598568 788876 901139 111312 151223 182739 205725 307127 792564 917255 113619 133827 151627 329440 627549 826353 932554 114924 133928 153429 185139 209733 638049 856766 115320 134324 156637 210617 397760 645551 857263 984472 135932 157336 192439 212318 414229 654754 859671 992774 117122 136833 158035 192540 214019 448347 665254 861456 118730 137128 158540 193239 214221 485757 672756 865363 118932 138231 159239 194140 214423 499566 677968 119833 139637 161731 217429 528446];
TRInSec=.72;
rawDataDir='/data/hcp-bluearc/OpenAccess/';
runNames={'tfMRI_EMOTION_RL' 'tfMRI_EMOTION_LR' 'tfMRI_GAMBLING_RL' 'tfMRI_GAMBLING_LR' 'tfMRI_LANGUAGE_RL' 'tfMRI_LANGUAGE_LR' 'tfMRI_MOTOR_RL' 'tfMRI_MOTOR_LR' 'tfMRI_RELATIONAL_RL' 'tfMRI_RELATIONAL_LR' 'tfMRI_SOCIAL_RL' 'tfMRI_SOCIAL_LR' 'tfMRI_WM_RL' 'tfMRI_WM_LR'};
TRsPerRun=[176 176 253 253 316 316 284 284 232 232 274 274 405 405];
%FSL-style timing file names (without .txt)
timingFilenamesByRun={{'fear' 'neut'}...
    {'fear' 'neut'}...
    {'win' 'loss'}...
    {'win' 'loss'}...
    {'story' 'math'}...
    {'story' 'math'}...
    {'cue' 'lf' 'rf' 'lh' 'rh' 't'}...
    {'cue' 'lf' 'rf' 'lh' 'rh' 't'}...
    {'relation' 'match'}...
    {'relation' 'match'}...
    {'mental' 'rnd'}...
    {'mental' 'rnd'}...
    {'0bk_body' '0bk_faces' '0bk_places' '0bk_tools' '2bk_body' '2bk_faces' '2bk_places' '2bk_tools'}...
    {'0bk_body' '0bk_faces' '0bk_places' '0bk_tools' '2bk_body' '2bk_faces' '2bk_places' '2bk_tools'}};
    
%Note from Greg Burgess regarding the Emotion condition:
%The EMOTION scans are in fact only 176TRs long, due to an late-detected bug. The EVs are effectively correct. Even though they code a longer time period, they will yield the correct GLM in FSL, because timepoints occurring after the end of the scan are ignored by feat. What is definitely incorrect is that the user guide should mention that both the EMOTION scan protocols and task paradigm were inadvertently truncated ~15 seconds earlier than initially intended. The good news is that the task design is powerful enough that we only need an additional 0.01% to 0.02% signal change to make up for the final 15 seconds of task. In other words, we're really not taking a hit from the lost data.

%Converting to AFNI format timing files
% HCPBasicBlockFactor:
% [dmat dmatLabels] = ReffuncToAFNI(reffuncsAll,{'100307','120212','162329','196144','221319','530635','702133','885975','123117','142828','163432','197550','249947','552544','729557','887373','103414','124422','143325','167743','199150','250427','559053','732243','889579','103515','125525','255639','579665','734045','894673','103818','128632','149337','172332','200614','293748','581349','748258','896778','105115','129028','149539','175439','201111','298051','585862','753251','896879','110411','130013','150423','205119','304020','598568','788876','901139','111312','151223','182739','205725','307127','792564','917255','113619','133827','151627','329440','627549','826353','932554','114924','133928','153429','185139','209733','638049','856766','115320','134324','156637','210617','397760','645551','857263','984472','135932','157336','192439','212318','414229','654754','859671','992774','117122','136833','158035','192540','214019','448347','665254','861456','118730','137128','158540','193239','214221','485757','672756','865363','118932','138231','159239','194140','214423','499566','677968','119833','139637','161731','217429','528446'},...
%     14, [1:14], [176 176 253 253 316 316 284 284 232 232 274 274 405 405], .72, {'EMOTION_fear','EMOTION_neut','GAMBLINGwin','GAMBLINGloss',...
%     'LANGUAGEstory','LANGUAGEmath','MOTORcue','MOTORlf','MOTORrf','MOTORlh','MOTORrh','MOTORt',...
%     'RELATIONALrelation','RELATIONALmatch','SOCIALmental','SOCIALrnd','WM0bk_body','WM0bk_faces','WM0bk_places','WM0bk_tools','WM2bk_body','WM2bk_faces','WM2bk_places','WM2bk_tools'},...
%     {'tfMRI_EMOTION\w*_fear\w*','tfMRI_EMOTION\w*_neut\w*','tfMRI_GAMBLING\w*win\w*','tfMRI_GAMBLING\w*loss\w*',...
%     'tfMRI_LANGUAGE\w*story\w*','tfMRI_LANGUAGE\w*math\w*','tfMRI_MOTOR\w*cue\w*','tfMRI_MOTOR\w*lf\w*','tfMRI_MOTOR\w*rf\w*','tfMRI_MOTOR\w*lh\w*','tfMRI_MOTOR\w*rh\w*','tfMRI_MOTOR\w*t\w*',...
%     'tfMRI_RELATIONAL\w*relation\w*','tfMRI_RELATIONAL\w*match\w*','tfMRI_SOCIAL\w*mental\w*','tfMRI_SOCIAL\w*rnd\w*','tfMRI_WM\w*0bk_body\w*','tfMRI_WM\w*0bk_faces\w*','tfMRI_WM\w*0bk_places\w*','tfMRI_WM\w*0bk_tools\w*','tfMRI_WM\w*2bk_body\w*','tfMRI_WM\w*2bk_faces\w*','tfMRI_WM\w*2bk_places\w*','tfMRI_WM\w*2bk_tools\w*'},...
%     1, 1, 'HCPBasicBlockFactor');



%% Create reference function by run and subject
reffuncsAll={};
reffuncsBySubjByRun=cell(length(subjNums),1);
for subjNum = 1:length(subjNums)
    for runNum = 1:length(runNames)
        reffuncs_thisrun=repmat({runNames{runNum}},TRsPerRun(runNum),1);
        %Go through each timing file, adding labels to TRs when indicated by file
        for timingFileNum=1:length(timingFilenamesByRun{runNum})
            filename=[timingFilenamesByRun{runNum}{timingFileNum} '.txt'];
            fid=fopen([rawDataDir num2str(subjNums(subjNum)) '/MNINonLinear/Results/' runNames{runNum} '/EVs/' filename]);
            %Extract timing info as 3 x number of events matrix. First of 3 is onset, second is duration, third is modulation (ignored)
            timingInfo=fscanf(fid,'%f %f %f', [3 inf]);
            fclose(fid);
            %Convert timing to reference function labels
            for eventNum=1:size(timingInfo,2)
                onset=timingInfo(1,eventNum);
                onsetTR=round(onset/TRInSec)+1;
                stop=onset+timingInfo(2,eventNum);
                stopTR=round(stop/TRInSec)+1;
                for trNum=onsetTR:stopTR
                    %Append label to this TR in the reffunc
                    if trNum <= TRsPerRun(runNum)
                        reffuncs_thisrun{trNum}=[reffuncs_thisrun{trNum} '_' timingFilenamesByRun{runNum}{timingFileNum}];
                    end
                end
            end
        end
        reffuncsBySubjByRun{subjNum}{runNum}=reffuncs_thisrun;
        reffuncsAll=[reffuncsAll; reffuncs_thisrun];
    end
end

save([mfilename '_' datestr(now) '.mat'], 'reffuncsAll','reffuncsBySubjByRun');