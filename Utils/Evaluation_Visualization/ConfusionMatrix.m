function [TPR, FPR, TNR, FNR, F1, ACC] = ConfusionMatrix(y, y_pred, show)

% Confusion Matrix
% Usage: ConfusionMatrix(y, y_pred, Num_class)

m = [y y_pred];

P = nnz(y(:,1)==+1);
N = nnz(y(:,1)==-1);

TP = sum((m(:,1)==+1 & m(:,2)==+1));
FP = sum((m(:,1)==-1 & m(:,2)==+1));
TN = sum((m(:,1)==-1 & m(:,2)==-1));
FN = sum((m(:,1)==+1 & m(:,2)==-1));

TPR = TP/P; TNR = TN/N;
FPR = FP/N; FNR = FN/P;

ACC = (TP+TN)/(TP+TN+FP+FN);
SEN = TP/(FN+TP);
SPE = TN/(TN+FP);

PREC = TP/(TP+FP);
RECA = TP/(TP+FN);

F1 = 2*1/(1/PREC+1/RECA);

CM = [TPR FNR; FPR TNR];
rowNames = {'AP','AN'};
colNames = {'PP','PN'};
Con_Matr = array2table(CM,'RowNames',rowNames,'VariableNames',colNames);

if isequal(show,'on')

    disp(Con_Matr)
    
    disp(['Accuracy ', num2str(ACC)])
    disp(['Specificity ', num2str(SPE)])
    disp(['Sensitivity ', num2str(SEN)])
    disp(['F1 ', num2str(F1)])

end
