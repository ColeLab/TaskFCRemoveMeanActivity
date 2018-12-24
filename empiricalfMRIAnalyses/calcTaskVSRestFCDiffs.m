function [taskFCvRestFC_diffs, taskFCvRestFC_diffs_tvals, taskFCvRestFC_diffs_pvals, percentSigValuesFDR_total, percentSigValuesFDR_pos, percentSigValuesFDR_neg, percentSigValuesFDR_total_bytask, percentSigValuesFDR_pos_bytask, percentSigValuesFDR_neg_bytask, FDRthreshByTask] = calcTaskVSRestFCDiffs(taskConnMatrix, restConnMatrix)
    taskFCvRestFC_diffs=zeros(size(taskConnMatrix,1),size(taskConnMatrix,2),size(taskConnMatrix,3));
    taskFCvRestFC_diffs_tvals=zeros(size(taskFCvRestFC_diffs));
    taskFCvRestFC_diffs_pvals=ones(size(taskFCvRestFC_diffs));
    upperTriangle=(tril(ones(size(taskConnMatrix,1),size(taskConnMatrix,2)))==0);
    taskFCvRestFC_diffs_pvals_vector=ones(sum(sum(upperTriangle))*size(taskFCvRestFC_diffs_pvals,3),1);
    taskFCvRestFC_diffs_Diffvals_vector=ones(sum(sum(upperTriangle))*size(taskFCvRestFC_diffs_pvals,3),1);
    percentSigValuesFDR_total_bytask=zeros(size(taskConnMatrix,3),1);
    percentSigValuesFDR_pos_bytask=zeros(size(taskConnMatrix,3),1);
    percentSigValuesFDR_neg_bytask=zeros(size(taskConnMatrix,3),1);
    FDRthreshByTask=zeros(size(taskConnMatrix,3),1);
    pval_count=1;
    for taskNum=1:size(taskConnMatrix,3)
        disp(['Task ' num2str(taskNum)])
        pvals_vector_thistask=ones(sum(sum(upperTriangle)),1);
        Diffvals_vector_thistask=zeros(sum(sum(upperTriangle)),1);
        pval_count_thistask=1;
        taskFCvRestFC_diffs(:,:,taskNum)=tanh(mean(atanh(taskConnMatrix(:,:,taskNum,:)),4))-tanh(mean(atanh(restConnMatrix(:,:,1,:)),4));
        for rowVal=1:size(taskConnMatrix,1)
            for colVal=1:size(taskConnMatrix,2)
                if rowVal>colVal
                    [~,p,~,stats]=ttest(atanh(taskConnMatrix(colVal,rowVal,taskNum,:)),atanh(restConnMatrix(colVal,rowVal,1,:)));
                    taskFCvRestFC_diffs_tvals(colVal,rowVal,taskNum)=stats.tstat;
                    taskFCvRestFC_diffs_tvals(rowVal,colVal,taskNum)=stats.tstat;
                    taskFCvRestFC_diffs_pvals(colVal,rowVal,taskNum)=p;
                    taskFCvRestFC_diffs_pvals(rowVal,colVal,taskNum)=p;
                    taskFCvRestFC_diffs_Diffvals_vector(pval_count)=mean(atanh(taskConnMatrix(colVal,rowVal,taskNum,:)))-mean(atanh(restConnMatrix(colVal,rowVal,1,:)));
                    taskFCvRestFC_diffs_pvals_vector(pval_count)=p;
                    Diffvals_vector_thistask(pval_count_thistask)=taskFCvRestFC_diffs_Diffvals_vector(pval_count);
                    pvals_vector_thistask(pval_count_thistask)=p;
                    pval_count_thistask=pval_count_thistask+1;
                    pval_count=pval_count+1;
                end
            end
        end
        %Calculate detection rate for this task only
        [p_FDRthresh_thistask]= FDR(pvals_vector_thistask,0.05);
        if isempty(p_FDRthresh_thistask)
            p_FDRthresh_thistask=0;
        end
        FDRthreshByTask(taskNum)=p_FDRthresh_thistask;
        sigValuesFDR_thistask=pvals_vector_thistask<p_FDRthresh_thistask;
        percentSigValuesFDR_total_bytask(taskNum)=100*(sum(sigValuesFDR_thistask))/length(pvals_vector_thistask);
        numSigValuesFDR_pos_thistask=sum(sigValuesFDR_thistask.*(Diffvals_vector_thistask>0));
        percentSigValuesFDR_pos_bytask(taskNum)=100*numSigValuesFDR_pos_thistask/length(Diffvals_vector_thistask);
        numSigValuesFDR_neg_thistask=sum(sigValuesFDR_thistask.*(Diffvals_vector_thistask<0));
        percentSigValuesFDR_neg_bytask(taskNum)=100*numSigValuesFDR_neg_thistask/length(Diffvals_vector_thistask);
        disp(['Percentage of significant task FC changes from rest across all connections, task ' num2str(taskNum) ': ' num2str(percentSigValuesFDR_total_bytask(taskNum))])
        disp(['Percentage of significant task FC changes from rest across all connections FC INCREASES, task ' num2str(taskNum) ': ' num2str(percentSigValuesFDR_pos_bytask(taskNum))])
        disp(['Percentage of significant task FC changes from rest across all connections FC DECREASES, task ' num2str(taskNum) ': ' num2str(percentSigValuesFDR_neg_bytask(taskNum))])
    end
    %Calculate FDR correction for multiple comparisons
    %[p_FDR] = mafdr(taskFCvRestFC_diffs_pvals_vector);
    [p_FDRthresh]= FDR(taskFCvRestFC_diffs_pvals_vector,0.05);
    numSigValuesFDR=sum(taskFCvRestFC_diffs_pvals_vector<p_FDRthresh);
    percentSigValuesFDR_total=100*numSigValuesFDR/length(taskFCvRestFC_diffs_pvals_vector);
    disp(['Percentage of significant task FC changes from rest across all connections and tasks: ' num2str(percentSigValuesFDR_total)])
    sigValuesFDR=taskFCvRestFC_diffs_pvals_vector<p_FDRthresh;
    numSigValuesFDR_pos=sum(sigValuesFDR.*(taskFCvRestFC_diffs_Diffvals_vector>0));
    percentSigValuesFDR_pos=100*numSigValuesFDR_pos/length(taskFCvRestFC_diffs_pvals_vector);
    disp(['Percentage of significant task FC changes from rest across all connections and tasks, FC INCREASES ONLY: ' num2str(percentSigValuesFDR_pos)])
    numSigValuesFDR_neg=sum(sigValuesFDR.*(taskFCvRestFC_diffs_Diffvals_vector<0));
    percentSigValuesFDR_neg=100*numSigValuesFDR_neg/length(taskFCvRestFC_diffs_pvals_vector);
    disp(['Percentage of significant task FC changes from rest across all connections and tasks, FC DECREASES ONLY: ' num2str(percentSigValuesFDR_neg)])
    
 
end