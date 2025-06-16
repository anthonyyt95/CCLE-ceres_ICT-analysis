%% PCA analysis based on ORAI1, ORAI2, ORAI3, STIM1, & STIM2


data = ccleICT_a_hemato.values;
genes = ccleICT_a_hemato.genes;
cells = ccleICT_a_hemato.cells;
annot = annot_hemato;
gofinterest = {'ORAI1' 'ORAI2' 'ORAI3' 'STIM1' 'STIM2'};
highlight = 'hodgkin';
highlight_ind = 5;      


color_trip = [0 0 0;
    200 200 200];
x_lim = 'N/A';
y_lim = 'N/A';
axis_width = 3;
barline_width = 3;
markersize = 500
fig_size = [10 10 800 700];


tic
% Acquires expression info for genes of relevance
pdata = [];
pgenes = [];
for i = [1:length(gofinterest)]
    ginterest = gofinterest{i};
    
    for j = [1:length(genes(:,1))]
        dgene = genes{j,1};
        
        if lower(string(ginterest)) == lower(string(dgene))
            pdata = [pdata; data(j,:)];
            pgenes = [pgenes; upper(dgene)];
            break
        end
    end
end
data = pdata;
genes = pgenes


% PCA analysis (to acquire 'score')
[coeff, score, latent] = pca(data');

% Plots PCA analyses
figure('Renderer', 'painters', 'Position', fig_size);
hold
xlabel('PC1')
ylabel('PC2') 
zlabel('PC3')
colors = color_trip/255;

if lower(string(highlight)) == 'n/a'
    scatter(score(:,1),score(:,2),markersize,colors(1,:),'.');
else
    
    % Finds indices of cell lines corresponding to the query
    group = [];
    for i = [1:length(annot(:,1))]
        val = annot{i,highlight_ind};
        
        present = strfind(lower(string(val)), lower(string(highlight)));
        if not(isempty(present))
            group = [group i];
            continue
        end
    end
    
    score1x = score(group,1);
    score1y = score(group,2);
    
    score(group,:) = [];
    
    scatter(score(:,1), score(:,2), markersize, colors(2,:), '.');
    scatter(score1x, score1y, markersize, colors(1,:), '.');
end



if isnumeric(x_lim)
    xlim([x_lim]);
end
if isnumeric(y_lim)
    ylim([y_lim]);
end
set(gca,'linewidth',axis_width)
set(gca, 'TickDir', 'out')
title(highlight);

toc
%set(gca, 'xtick', [])


%% Ranks cancers by a cluster of genes

data = ccleICT_a_haemto.values;
genes = ccleICT_a_haemto.genes;
cells = ccleICT_a_haemto.cells;
annot = annot_hemato;

% % Acquires indices of genes of interest
% for i = [1:length(genes(:,1))]
%     val = genes{i,1};
%     
%     if lower
% 
% output = {};
% for i = [1:length(annot(:,1))]
