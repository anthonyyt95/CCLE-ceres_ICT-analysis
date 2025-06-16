%% Performs DEG between heamtopoietic cancers and all others

data = ccleICT_a.values;
genes = ccleICT_a.genes;
cells = ccleICT_a.cells;
annot = ccleICT_a.annot;
query = 'lymphoid';

tic
% Acquires indices of all cells associated with 'query' term
ind = [];
for i = [1:length(annot(:,1))]
    cell = lower(annot{i,2});
    
    loc = strfind(cell, lower(query));
    if isempty(loc)
        continue
    else
        ind = [ind; i];
    end
end
qdata = data(:,ind);
data(:,ind) = [];

% Performs differential analysis
output = {};
for i = [1:length(genes(:,1))]
    val1 = data(i,:);
    val2 = qdata(i,:);
    
    [~, p] = ttest2(val1, val2);
    logfc = log2(mean(val2)/mean(val1));
    avg_expr = mean([val2 val1]);
    
    output{i,1} = genes{i,1};
    output{i,2} = logfc;
    output{i,3} = p;
    output{i,4} = avg_expr;
end

% Acquires FDR adjusted p values
pvals = cell2mat(output(:,3));
fdr_vals = mafdr(pvals);
output(:,5) = num2cell(fdr_vals);
output = sortrows(output,5);

% Organizes result
heading = {'Genes' 'Log Fold-change' 'P-value' 'Average Expression' 'FDR'};
output = [heading; output];

toc


clear annot cell cells data genes i ind loc qdata query avg_expr logfc p
clear val1 val2 fdr_vals heading pvals
    
