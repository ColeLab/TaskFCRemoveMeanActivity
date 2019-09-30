function [designMatrix designMLabels] = ReffuncToAFNI(RefFunc, SubjectArray, NumRunsTotal, RunsToInclude, RunLengths, TRInSecs, CondNames, CondSearchStrings, RegExpr, OutputToFiles, AnalysisName)
%
%This function takes an array of condition labels (reference function) and creates a series of text files (one per stimulus function) for use in AFNI
%
%Created by Michael W. Cole
%http://www.mwcole.net
%mwcole@mwcole.net
%
%
% Examples:
%
% HCPBasicBlockFactor:
% [dmat dmatLabels] = ReffuncToAFNI(reffuncsAll,{'100307','120212','162329','196144','221319','530635','702133','885975','123117','142828','163432','197550','249947','552544','729557','887373','103414','124422','143325','167743','199150','250427','559053','732243','889579','103515','125525','255639','579665','734045','894673','103818','128632','149337','172332','200614','293748','581349','748258','896778','105115','129028','149539','175439','201111','298051','585862','753251','896879','110411','130013','150423','205119','304020','598568','788876','901139','111312','151223','182739','205725','307127','792564','917255','113619','133827','151627','329440','627549','826353','932554','114924','133928','153429','185139','209733','638049','856766','115320','134324','156637','210617','397760','645551','857263','984472','135932','157336','192439','212318','414229','654754','859671','992774','117122','136833','158035','192540','214019','448347','665254','861456','118730','137128','158540','193239','214221','485757','672756','865363','118932','138231','159239','194140','214423','499566','677968','119833','139637','161731','217429','528446'},...
%     14, [1:14], [176 176 253 253 316 316 284 284 232 232 274 274 405 405], .72, {'EMOTION_fear','EMOTION_neut','GAMBLINGwin','GAMBLINGloss',...
%     'LANGUAGEstory','LANGUAGEmath','MOTORcue','MOTORlf','MOTORrf','MOTORlh','MOTORrh','MOTORt',...
%     'RELATIONALrelation','RELATIONALmatch','SOCIALmental','SOCIALrnd','WM0bk_body','WM0bk_faces','WM0bk_places','WM0bk_tools','WM2bk_body','WM2bk_faces','WM2bk_places','WM2bk_tools'},...
%     {'tfMRI_EMOTION\w*_fear\w*','tfMRI_EMOTION\w*_neut\w*','tfMRI_GAMBLING\w*win\w*','tfMRI_GAMBLING\w*loss\w*',...
%     'tfMRI_LANGUAGE\w*story\w*','tfMRI_LANGUAGE\w*math\w*','tfMRI_MOTOR\w*cue\w*','tfMRI_MOTOR\w*lf\w*','tfMRI_MOTOR\w*rf\w*','tfMRI_MOTOR\w*lh\w*','tfMRI_MOTOR\w*rh\w*','tfMRI_MOTOR\w*t\w*',...
%     'tfMRI_RELATIONAL\w*relation\w*','tfMRI_RELATIONAL\w*match\w*','tfMRI_SOCIAL\w*mental\w*','tfMRI_SOCIAL\w*rnd\w*','tfMRI_WM\w*0bk_body\w*','tfMRI_WM\w*0bk_faces\w*','tfMRI_WM\w*0bk_places\w*','tfMRI_WM\w*0bk_tools\w*','tfMRI_WM\w*2bk_body\w*','tfMRI_WM\w*2bk_faces\w*','tfMRI_WM\w*2bk_places\w*','tfMRI_WM\w*2bk_tools\w*'},...
%     1, 1, 'HCPBasicBlockFactor');

% Versions:
% 1.0b - First functioning version
% 1.1b - Adding stimtiming file creation functionality
% 1.2b - Changing NumRuns to RunsToInclude (to allow for odd-even splits)
% 1.3b - Changing RunLengths to RunLengthss (to allow for runs of different lengths)

NumRuns=length(RunsToInclude);

numTRsTotal = length(SubjectArray) * sum(RunLengths);
designMatrix = zeros(numTRsTotal, length(CondSearchStrings));

%Iterating through conditions
for condNum = 1:length(CondSearchStrings)

    condLabel = CondSearchStrings{condNum};
    trCount = 1;
    dmTRCount = 1;

    %Iterating through subjects
    for subjNum = 1:length(SubjectArray)
        subjName = SubjectArray(subjNum);

        %Iterating through runs
        for runNum = 1:NumRunsTotal
            
            %Only include this run if in RunsToInclude, otherwise increment TR count
            if isempty(find(RunsToInclude==runNum, 1))
                trCount = trCount + RunLengths(runNum);
            else

                %Iterating through TRs
                for runTR = 1:RunLengths(runNum)

                    %If condLabel matches Reffunc at this TR, add a 1 to the output, otherwise 0
                    %If RegExpr is true, then include partial string matches (based on regular expressions), otherwise only include exact matches
                    if RegExpr
                        %partial match
                        %match = ~isempty(strfind(RefFunc{trCount}, condLabel));
                        match=~isempty(regexp(RefFunc{trCount},condLabel,'once'));
                        %match=find(cellfun('isempty',regExpList)==0);
                    else
                        %exact match
                        match = ~isempty(strmatch(RefFunc{trCount}, condLabel,'exact'));
                    end

                    %Save match result to design matrix
                    designMatrix(dmTRCount, condNum) = match;

                    trCount = trCount + 1;
                    dmTRCount = dmTRCount + 1;

                end
            end
        end
    end
end

dmat = designMatrix;
designMLabels = CondNames;

if(OutputToFiles)

    %Writing stimulus file; for use with -stim_file option in 3ddeconvolve
    for colNum = 1:length(CondSearchStrings)
        trCount = 1;

        for subjCount = 1:length(SubjectArray)
            disp(['Saving AFNI unconvolved stim files (1D file) for subject ' SubjectArray{subjCount}]);

            writeName = [SubjectArray{subjCount} '_stimfile_' AnalysisName '_EV' num2str(colNum) '_' CondNames{colNum} '.1D'];
            disp(['Writing ' writeName ' (used with -stim_file option in 3dDeconvolve)'])
            fidW = fopen(writeName, 'wt');

            for runCount = 1:NumRuns
                %thisDesignMatrix = designMatrix{subjCount}{runCount};

                %Write vector
                for row = trCount:(trCount+RunLengths(runCount)-1)
                    rowString = [num2str(designMatrix(row, colNum), '%6f')];
                    fprintf(fidW, '%s\n',rowString);
                    trCount = trCount + 1;
                end

            end
            fclose(fidW);
            %Writing timing file (stim timing file); for use with -stim_times option in 3ddeconvolve
            disp(['Writing stime_' writeName ' (used with -stim_times option in 3dDeconvolve)'])
            if length(unique(RunLengths)) == 1
                eval(['!make_stim_times.py -files ' writeName ' -prefix stime_' writeName ' -tr ' num2str(TRInSecs) ' -nruns ' num2str(NumRuns) ' -nt ' num2str(RunLengths(1))])
            else
                disp('Warning: Run lengths differ, so stime_ files aren''t created')
            end
        end
    end

%     %Writing timing file (t-file); for use with -stim_times option in 3ddeconvolve
%     for colNum = 1:length(CondSearchStrings)
%         trCount = 1;
% 
%         for subjCount = 1:length(SubjectArray)
%             disp(['Saving AFNI stim timing files (1D file; used with -stim_times option in 3ddeconvolve) for subject ' SubjectArray{subjCount}]);
% 
%             writeName = [SubjectArray{subjCount} '_stimfile_' AnalysisName '_EV' num2str(colNum) '_' CondNames{colNum} '_tfile.1D'];
%             disp(['Writing ' writeName])
%             fidW = fopen(writeName, 'wt');
% 
%             for runCount = 1:NumRuns
%                 %thisDesignMatrix = designMatrix{subjCount}{runCount};
% 
%                 %Write vector
%                 for rowNum = trCount:(trCount+RunLengths-1)
%                     rowString = [num2str(designMatrix(row, colNum), '%6.1f')];
%                     if strcmp(rowString,'1')
%                         rowStringOut = num2str(rowNum);
%                         fprintf(fidW, '%s\n',rowStringOut);
%                     end
%                  writeName   trCount = trCount + 1;
%                 end
% 
%             end
%             fclose(fidW);
%         end
%     end
    
    for subjCount = 1:length(SubjectArray)
        writeName = [SubjectArray{subjCount} '_' AnalysisName '_ConcatRef.1D'];
        disp(['Writing ' writeName])
        fidW = fopen(writeName, 'wt');

        %Writing file for -concat option in 3ddevonvolve (removes trends for individual runs)
        lineCount=0;
        for runNum=1:NumRuns
            fprintf(fidW, '%s\n',num2str(lineCount));
            lineCount = lineCount + RunLengths(runNum);
        end

        fclose(fidW);
    end
end



