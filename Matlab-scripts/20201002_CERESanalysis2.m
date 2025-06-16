%% Queries several cell lines & ranks ICT genes according to CERES score

load('20190426_workspace_CRISPRi.mat')
glistName = 'E:\Documents\NYU\NYU Langone\PhD\Feske Lab\Experiments\04.10.19_CancerDatabases&Analysis\Data\2020.10.02_ICTlist\GeneList.ICT.MT.txt';

data = ceres.values;
cells = ceres.cells;
genes = ceres.genes;
annot = ceres.annot;
%query = {'B-cell, hodgkins',
%    'B-cell, Non';
%    'DLBCL';
%    'Burkitts',
%    'multiple myeloma'};
query = {'Acute Lymphoblastic Leukemia (ALL), B-cell'};
%query = '_';

tic
% Acquires gene list
glist = readtable(glistName,'ReadVariableNames',false);
glist = table2cell(glist);

% Extracts gene data for all/most ICTs from 'glist' into 'ict_data' &
% 'ict_gene'
ict_data = [];
ict_genes = {};
ict_class = {};
genes = lower(genes);
missing = {};
for i = [1:length(glist(:,1))]
    icts = strsplit(glist{i,1},',');
    
    % Searches for all aliases of the current ICT
    match = 0;
    for j = [1:length(icts)]
        ict = lower(string(icts{j}));
        
        loc = find(strcmp(ict, genes));
        if not(isempty(loc))
            ict_genes = [ict_genes; icts(2)];
            ict_data = [ict_data; data(loc,:)];
            ict_class = [ict_class; glist{i,2}];
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
for i = [1:length(annot(:,1))]
    annot{i,5} = char(string(annot{i,5}));
end

group = [];
output = {};
cellsOut = {};
for i = [1:length(annot(:,1))]
    val = annot{i,5};
    
    for j = [1:length(query)]
        qTerm = query{j};
        
        present = strfind(lower(string(val)), lower(string(qTerm)));
        if isempty(present)
            continue
        else
            group = [group i];
            cellsOut = [cellsOut; val];
            break
        end
    end
end
cellsOut = unique(cellsOut);

sdata = ict_data(:,group);
sdata = mean(sdata,2);

zdata = zict_data(:,group);
zdata = mean(zdata,2);

output = {};
output = [ict_genes ict_class num2cell(sdata) num2cell(zdata)];
sort_by_ceres = sortrows(output,3);
sort_by_zscore = sortrows(output,4);
toc

clear annot cells data genes glist group i ict loc ogene_list present output
clear query query_ind step val sdata zdata zict_data zlist ict_data ict_genes
clear match j icts glistName ict_class qTerm

%% Acquires CCLE expression data

load('20190426_workspace_ccleAnnotation_Hematopoietic.mat');
%query = {'B-cell';
%    'DLBCL';
%    'Burkitts'};
query = {'Acute Lymphoblastic Leukemia (ALL), B-cell'};
glistName = 'E:\Documents\NYU\NYU Langone\PhD\Feske Lab\Experiments\04.10.19_CancerDatabases&Analysis\Data\2020.10.02_ICTlist\GeneList.ICT.MT.txt';
geneInput = sort_by_ceres(:,1);

% Acquires gene list & prepares expression data
glist = readtable(glistName,'ReadVariableNames',false);
glist = table2cell(glist);

values = ccle_hemato_full.values;
annot = ccle_hemato_full.annot;
genes = ccle_hemato_full.genes;

% Acquires associated cancer cell lines
group = [];
output = {};
for i = [1:length(annot(:,1))]
    val = annot{i,5};
    
    for j = [1:length(query)]
        qTerm = query{j};
        
        present = strfind(lower(string(val)), lower(string(qTerm)));
        if isempty(present)
            continue
        else
            group = [group i];
            break
        end
    end
end
values = values(:,group);

% Associates expression data to CERES data (from above output)
valuesOut = [];
for i = [1:length(geneInput(:,1))]
    gene = geneInput{i,1};
    
    loc = find(strcmp(gene, genes(:,1)));
    if not(isempty(loc))
        valuesOut(i,1) = mean(values(loc(1),:));
    end
end
output = [geneInput, num2cell(valuesOut)];

clear aliases annot gene geneInput genes genesOut glist glistName i j loc
clear match missing query values valuesOut group present qTerm val

%% Acquires top genes whose dependency correlates with that of a query gene

load('20190426_workspace_CRISPRi.mat')
query = 'SLC39A7';

data = ceres.values;
cells = ceres.cells;
genes = ceres.genes;
annot = ceres.annot;

% Acquires query gene info
loc = find(strcmp(genes,query));
qVals = data(loc,:);

% Acquires correlations for all other genes
rhoOut = [];
pOut = [];
for i = [1:length(genes(:,1))]
    val = data(i,:);
    
    [r p] = corr(qVals',val','type','Pearson');
    rhoOut = [rhoOut; r];
    pOut = [pOut; p];
end
output = [genes num2cell(rhoOut) num2cell(pOut)];
output = sortrows(output,2,'desc');

clear query data cells genes annot loc qVals rhoOut pOut i val r p

%% Cleans up ICT list

listName = 'E:\Documents\NYU\NYU Langone\PhD\Feske Lab\Experiments\04.10.19_CancerDatabases&Analysis\Data\2020.10.02_ICTlist\ICT.MT.List.txt';

% Load data
data = readtable(listName,'ReadVariableNames',false);
data = table2cell(data);

dataAliases = data(:,[1:4]);
dataClasses = data(:,5);

uniqueGenes = unique(data(:,1));

% Cleans up list
aliasesOut = {};
classOut = {};
for i = [1:length(uniqueGenes(:,1))]
    gene = uniqueGenes(i,1);
    loc = find(strcmp(data(:,1),gene));
    if not(isempty(loc))
        aliases = dataAliases(loc(1),:);
        aliases = strjoin(aliases,',');
        
        class = unique(dataClasses(loc,1));
        class = strjoin(class,'/');
        
        aliasesOut = [aliasesOut; aliases];
        classOut = [classOut; class];
    end
end

output = [aliasesOut, classOut];

clear listName data dataAliases dataClasses uniqueGenes aliasesOut
clear classOut i gene loc aliases class 




