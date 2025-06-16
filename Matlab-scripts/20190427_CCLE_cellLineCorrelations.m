%% Correlation analysis of intra-categorical cell lines

data = ccleICT_a.values;
cells = ccleICT_a.cells;
genes = ccleICT_a.genes;
annot = ccle_annot;
query = 'Diffuse Large B-cell Lymphoma (DLBCL)';
query_ind = 5;

% Finds indices of cell lines corresponding to the query
group = [];
for i = [1:length(annot(:,1))]
    val = annot{i,query_ind};

    present = strfind(lower(string(val)), lower(string(query)));
    if not(isempty(present))
        group = [group i];
        continue
    end
end

cell_val = data(:,group);
cell_line = annot(group,2);
for i = [1:length(cell_line)]
    val = cell_line{i};
    ind = strfind(val,'_');
    val = char(val);
    val = val([1:ind(1)-1]);
    cell_line{i} = val;
end

coeff = corrcoef(cell_val);

output = [cell_line num2cell(coeff)];
heading = [{''} cell_line'];
output = [heading; output];

heatmap(coeff)




            