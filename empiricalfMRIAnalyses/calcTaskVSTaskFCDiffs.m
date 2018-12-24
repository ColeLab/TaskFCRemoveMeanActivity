function [taskFCvTaskFC_diffs, taskFCvTaskFC_diffs_tvals, taskFCvTaskFC_diffs_pvals, percentSigValuesFDR_total, percentSigValuesFDR_pos, percentSigValuesFDR_neg, FDRThreshP05] = calcTaskVSTaskFCDiffs(task1ConnMatrix, task2ConnMatrix)
    taskFCvTaskFC_diffs=zeros(size(task1ConnMatrix,1),size(task1ConnMatrix,2),size(task1ConnMatrix,3));
    taskFCvTaskFC_diffs_tvals=zeros(size(taskFCvTaskFC_diffs));
    taskFCvTaskFC_diffs_pvals=ones(size(taskFCvTaskFC_diffs));
    upperTriangle=(tril(ones(size(task1ConnMatrix,1),size(task1ConnMatrix,2)))==0);
    taskFCvTaskFC_diffs_pvals_vector=ones(sum(sum(upperTriangle))*size(taskFCvTaskFC_diffs_pvals,3),1);
    taskFCvTaskFC_diffs_Diffvals_vector=ones(sum(sum(upperTriangle))*size(taskFCvTaskFC_diffs_pvals,3),1);
    %percentSigValuesFDR_total_bytask=zeros(size(task1ConnMatrix,3),1);
    %percentSigValuesFDR_pos_bytask=zeros(size(task1ConnMatrix,3),1);
    %percentSigValuesFDR_neg_bytask=zeros(size(task1ConnMatrix,3),1);
    %FDRthreshByTask=zeros(size(task1ConnMatrix,3),1);
    pval_count=1;
    
    taskNum=1;
    pvals_vector_thistask=ones(sum(sum(upperTriangle)),1);
    Diffvals_vector_thistask=zeros(sum(sum(upperTriangle)),1);
    pval_count_thistask=1;
    taskFCvTaskFC_diffs(:,:,taskNum)=tanh(mean(atanh(task1ConnMatrix(:,:,taskNum,:)),4)-mean(atanh(task2ConnMatrix(:,:,1,:)),4));
    for rowVal=1:size(task1ConnMatrix,1)
        for colVal=1:size(task1ConnMatrix,2)
            if rowVal>colVal
                [~,p,~,stats]=ttest(atanh(task1ConnMatrix(colVal,rowVal,taskNum,:)),atanh(task2ConnMatrix(colVal,rowVal,1,:)));
                taskFCvTaskFC_diffs_tvals(colVal,rowVal,taskNum)=stats.tstat;
                taskFCvTaskFC_diffs_tvals(rowVal,colVal,taskNum)=stats.tstat;
                taskFCvTaskFC_diffs_pvals(colVal,rowVal,taskNum)=p;
                taskFCvTaskFC_diffs_pvals(rowVal,colVal,taskNum)=p;
                taskFCvTaskFC_diffs_Diffvals_vector(pval_count)=mean(atanh(task1ConnMatrix(colVal,rowVal,taskNum,:)))-mean(atanh(task2ConnMatrix(colVal,rowVal,1,:)));
                taskFCvTaskFC_diffs_pvals_vector(pval_count)=p;
                Diffvals_vector_thistask(pval_count_thistask)=taskFCvTaskFC_diffs_Diffvals_vector(pval_count);
                pvals_vector_thistask(pval_count_thistask)=p;
                pval_count_thistask=pval_count_thistask+1;
                pval_count=pval_count+1;
            end
        end
    end
    %Calculate detection rate for this task only
    %     [p_FDRthresh_thistask]= FDR(pvals_vector_thistask,0.05);
    %     if isempty(p_FDRthresh_thistask)
    %         p_FDRthresh_thistask=0;
    %     end
    %     FDRthreshByTask(taskNum)=p_FDRthresh_thistask;
    %     sigValuesFDR_thistask=pvals_vector_thistask<p_FDRthresh_thistask;
    %     percentSigValuesFDR_total_bytask(taskNum)=100*(sum(sigValuesFDR_thistask))/length(pvals_vector_thistask);
    %     numSigValuesFDR_pos_thistask=sum(sigValuesFDR_thistask.*(Diffvals_vector_thistask>0));
    %     percentSigValuesFDR_pos_bytask(taskNum)=100*numSigValuesFDR_pos_thistask/length(Diffvals_vector_thistask);
    %     numSigValuesFDR_neg_thistask=sum(sigValuesFDR_thistask.*(Diffvals_vector_thistask<0));
    %     percentSigValuesFDR_neg_bytask(taskNum)=100*numSigValuesFDR_neg_thistask/length(Diffvals_vector_thistask);
    %disp(['Percentage of significant task FC changes from rest across all connections, task ' num2str(taskNum) ': ' num2str(percentSigValuesFDR_total_bytask(taskNum))])
    %disp(['Percentage of significant task FC changes from rest across all connections FC INCREASES, task ' num2str(taskNum) ': ' num2str(percentSigValuesFDR_pos_bytask(taskNum))])
    %disp(['Percentage of significant task FC changes from rest across all connections FC DECREASES, task ' num2str(taskNum) ': ' num2str(percentSigValuesFDR_neg_bytask(taskNum))])
    
    %Calculate FDR correction for multiple comparisons
    %[p_FDR] = mafdr(taskFCvTaskFC_diffs_pvals_vector);
    [p_FDRthresh]= FDR(taskFCvTaskFC_diffs_pvals_vector,0.05);
    if isempty(p_FDRthresh)
        disp('No values survived correction for multiple comparisons')
        p_FDRthresh=0;
    end
    FDRThreshP05=p_FDRthresh;
    numSigValuesFDR=sum(taskFCvTaskFC_diffs_pvals_vector<p_FDRthresh);
    percentSigValuesFDR_total=100*numSigValuesFDR/length(taskFCvTaskFC_diffs_pvals_vector);
    disp(['Percentage of significant task FC changes from rest across all connections: ' num2str(percentSigValuesFDR_total)])
    sigValuesFDR=taskFCvTaskFC_diffs_pvals_vector<p_FDRthresh;
    numSigValuesFDR_pos=sum(sigValuesFDR.*(taskFCvTaskFC_diffs_Diffvals_vector>0));
    percentSigValuesFDR_pos=100*numSigValuesFDR_pos/length(taskFCvTaskFC_diffs_pvals_vector);
    disp(['Percentage of significant task FC changes from rest across all connections, FC INCREASES ONLY: ' num2str(percentSigValuesFDR_pos)])
    numSigValuesFDR_neg=sum(sigValuesFDR.*(taskFCvTaskFC_diffs_Diffvals_vector<0));
    percentSigValuesFDR_neg=100*numSigValuesFDR_neg/length(taskFCvTaskFC_diffs_pvals_vector);
    disp(['Percentage of significant task FC changes from rest across all connections, FC DECREASES ONLY: ' num2str(percentSigValuesFDR_neg)])
    
 
end