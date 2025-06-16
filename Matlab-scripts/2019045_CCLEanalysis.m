%% Acquires indices for all cell lines that match the query
%
% Input:
%   'annot'             1166-by-7 cell matrix containing cell-line
%                       information
%                           - workspace:
%                           '201804014_workspace_ccleAnnotations.mat'
%
%   'query'             Cell list of terms of interest
%   'col_ind'           Description column in 'annot' to be sifted through
%                           - col_ind = 4 --- Primary Disease
%                           - col_ind = 5 --- Subtype Disease
%
% Output:
%   'output'            X-by-4 cell matrix of cell-types corresponding to
%                       the query + relevant information
%
% Notes:
%   a. query can be a cell of multiple queries


annot = ccleICT_a_hemato.annot;
query = {'Acute Lymphoblastic Leukemia (ALL), T-cell'};
%query = {'Acute Lymphoblastic Leukemia (ALL), B-cell'};
col_ind = 5;

output = {};
step = 1;
for i = [1:length(query)]
    iquery = query{i};
    
    for j = [2:length(annot(:,1))]
        cancertype = annot{j,col_ind};
        
        loc = strfind(lower(string(cancertype)), lower(string(iquery)));
        if isempty(loc)
            continue
        else
            output{step,1} = j;
            output{step,2} = annot{j,4};
            output{step,3} = annot{j,5};
            output{step,4} = annot{j,2};
            step = step + 1;
        end
    end
end
heading = {'Index' 'Primary Disease' 'Subtype Disease' 'CCLE Name'};
output = [heading; output];
indices = cell2mat(output([2:end],1));

clearvars annot cancertype col_ind heading i iquery j loc output query step

%% Performs analysis (expression across cells x expression within cell)

tic

data = ccleICT_a_hemato.values;
genes = ccleICT_a_hemato.genes;
cells = ccleICT_a_hemato.cells;
cellcol = indices;
cell_info = ccle_cellinfo;
gene_info = ccle_geneinfo;
cell_title = 'Leukemia/lymphoma';

% Averages expression values if multiple inputs in 'cellcol' (& also
% divisor values)
if length(cellcol) > 1
    
    avged_data = mean(data(:,cellcol),2);
    
    data(:,cellcol(1)) = avged_data(:,1);
    data(:,cellcol([2:end])) = [];
end

divisions = cell_info(7,cellcol);
div_val = [];
for i = [1:length(divisions(1,:))]
    div_val = [div_val divisions{i}];
end
divisions = mean(div_val,2);

cmeans = mean(cell2mat(cell_info(3,cellcol)));
cstds = mean(cell2mat(cell_info(6,cellcol)));

cellcol = cellcol(1);

ict_data = data;
cell_data = [genes, num2cell(data(:,cellcol))];


%Compares expression of ICTs in cell of interest to that of all other cell-types
output1 = {};
step = 1;
for i = [1:length(ict_data(:,1))]
    gene = genes{i,1};
    [val_list, mu, sigma] = zscore(ict_data(i,:));
    
    cell_zscore = val_list(cellcol);
    withinCell_zscore = (ict_data(i,cellcol) - cmeans)/cstds;
    
        
    output1{step,1} = upper(gene);
    output1{step,2} = cell_zscore;
    output1{step,3} = ict_data(i,cellcol);
    output1{step,4} = withinCell_zscore;
    if cell_zscore < 0 & withinCell_zscore < 0
        output1{step,5} = cell_zscore * withinCell_zscore * -1;
    else
        output1{step,5} = cell_zscore * withinCell_zscore;
    end
    
    step = step + 1;
end
zcombined_output = [output1(:,1) output1(:,2) output1(:,4) output1(:,5)];

output1 = sortrows(output1,2,'desc');
output1(:,[4 5]) = [];
heading = {'Genes' 'Z-scores (across cells)' 'Actual values'};
output1 = [heading; output1];

zcombined_output = sortrows(zcombined_output,4,'desc');
heading = {'Genes' 'Across Cells (z-score)' 'Within Cells (z-score)' 'Multiplied z-score'};
zcombined_output = [heading; zcombined_output];


% Ranks ICT genes within cell of interest according to expression value
output2 = sortrows(cell_data,2,'desc');
for i = [1:length(output2(:,1))]
    val = output2{i,2};
    val_z = (val - cmeans)/cstds;
    
    
    output2{i,3} = val_z;
    

    divisions = cell_info{7,cellcol};
    step = 2;
    bin = 9;
    if val < divisions(2)
        output{i,4} = bin;
        continue
    else
        while val > divisions(step)
            step = step + 1;
            bin = bin - 1;
            
            if step > length(divisions)
                bin = 1;
                break
            end
        end
    end
         
    output2{i,4} = bin;
end
heading = {'Gene' 'Actual Value' 'Z-score (within cell)' 'Tile'};
output2 = [heading; output2];


toc
    
%% Ranks expression of ICTs within a queried cell line class
%
% Note: run 1st block to acquire 'indices'


data = ccleICT_a_hemato.values;
genes = ccleICT_a_hemato.genes;
cells = ccleICT_a_hemato.cells;
cellcol = indices_tall';

% Acquires mean expression values for each gene across corresponding cell
% lines
o_data = data(:,cellcol);
o_data = mean(o_data,2);
output = [genes num2cell(o_data)];
output = sortrows(output,2,'desc');


clear cellcol cells data genes o_data 






