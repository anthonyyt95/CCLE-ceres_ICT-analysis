%% Queries several cell lines & ranks ICT genes according to CERES score

data = ceres.values;
cells = ceres.cells;
genes = ceres.genes;
annot = ceres.annot;
glist = ICT862;
query = 'Acute Lymphoblastic Leukemia (ALL), B-cell';
%query = 'Acute Lymphoblastic Leukemia (ALL), T-cell';
%query = '_';
query_ind = 5;

tic
% Extracts gene data for all/most ICTs from 'glist' into 'ict_data' &
% 'ict_gene'
ict_data = [];
ict_genes = {};
genes = lower(genes);
missing = {};
for i = [1:length(glist(:,1))]
    icts = strsplit(glist{i,1},' ');
    
    % Searches for all aliases of the current ICT
    match = 0;
    for j = [1:length(icts)]
        ict = lower(string(icts{j}));
        
        loc = find(strcmp(ict, genes));
        if not(isempty(loc))
            ict_genes = [ict_genes; upper(char(ict))];
            ict_data = [ict_data; data(loc,:)];
            match = 1;
            break
        end
    end
    
    % Lists missing ICTs
    if match == 0
        missing = [missing; glist{i,1}];
    end
end

% Acquires z-scores for each gene across all cell lines
zict_data = ict_data;
for i = [1:length(ict_data(:,1))]
    zlist = zscore(ict_data(i,:));
    zict_data(i,:) = zlist;
end

% Acquires associated cancer cell lines
group = [];
output = {};
for i = [1:length(annot(:,1))]
    val = annot{i,query_ind};
    
    present = strfind(lower(string(val)), lower(string(query)));
    if isempty(present)
        continue
    else
        group = [group i];
    end
end

sdata = ict_data(:,group);
sdata = mean(sdata,2);

zdata = zict_data(:,group);
zdata = mean(zdata,2);

output = {};
output = [ict_genes num2cell(sdata) num2cell(zdata)];
sort_by_ceres = sortrows(output,2);
sort_by_zscore = sortrows(output,3);
toc


clear annot cells data genes glist group i ict loc ogene_list present output
clear query query_ind step val sdata zdata zict_data zlist ict_data ict_genes
clear match j icts

%% Acquires overlaps

list1 = a;
list2 = b;

list1 = lower(list1);
list2 = lower(list2);
output = {};
for i = [1:length(list1(:,1))]
    ict1 = list1{i,1};
    
    loc = find(strcmp(ict1, list2));
    if isempty(loc)
        continue
    else
        output = [output; upper(ict1)];
    end
end

clear list1 list2 i ict1 loc 



%% Create an annotation variable for CERES database

cells = ceres.cells;
annot = ccle_annot;

output = {};
for i = [1:length(cells(1,:))]
    dcell = cells{i};
    
    for j = [1:length(annot(:,1))]
        acell = annot{j,1};
        
        if lower(string(dcell)) == lower(string(acell))
            output(i,:) = annot(j,:);
            break
        end
    end
end


%% Acquires CERES information for a cell line & ranks it

glist = ictList_a;
data = ceres.values;
genes = ceres.text(:,1);
cells = ceres.text(1,:);
cell_line = 'ACH-000004';

for i = [1:length(cells)]
    dcells = cells{i};
    if lower(string(cell_line)) == lower(string(dcells))
        ind = i;
        break
    end
end

cdata = data(:,ind);
cdata = [genes num2cell(cdata)];

% Extracts ICT genes
% Extracts gene data for all/most ICTs from 'glist' into 'output'
output = {};
step = 1;
for i = [1:length(glist(:,1))]
    ict = glist{i,1};
    
    for j = [1:length(cdata(:,1))]
        data_gene = cdata{j,1};
        if lower(string(ict)) == lower(string(data_gene))
            output(step,:) = cdata(j,:);
            step = step + 1;
        else
            continue
        end
    end
end
heading = {'Gene' cell_line};
output = sortrows(output,2);
output = [heading; output];


%% Acquires CERES information for a queried gene

data = ceres.values;
genes = ceres.text(:,1);
cells = ceres.text(1,:);

query = 'trpm7';

output = {};
for i = [2:length(genes(:,1))]
    dgene = genes{i,1};
    
    if lower(string(query)) == lower(string(dgene))
        output = num2cell(data(i,:));
        output{1,1} = query;
        break
    end
end


%% Acquires CERES score for items in a 'list' according to 'ref'

list = output;
ref = r;

list = lower(list);
ref_genes = lower(ref(:,1));
output = {};
for i = [1:length(list(:,1))]
    ict = list{i,1};
    
    loc = find(strcmp(ict, ref_genes));
    score = ref{loc,2};
    
    output = [output; upper(ict) num2cell(score)];
end

clear list i ref ref_genes ict loc score





