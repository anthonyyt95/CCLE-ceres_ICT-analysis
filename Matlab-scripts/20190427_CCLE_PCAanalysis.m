%% PCA analysis of cancer data - 2D

data = ceres.values;
cells = ceres.cells;
genes = ceres.genes;
annot = ceres.annot;
highlight = 'lymphoid';
highlight_ind = 2;      

color_trip = [0 0 0;
    200 200 200];
x_lim = 'N/A';
y_lim = 'N/A';
axis_width = 3;
barline_width = 3;
markersize = 500
fig_size = [10 10 800 700];

tic

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
toc


if isnumeric(x_lim)
    xlim([x_lim]);
end
if isnumeric(y_lim)
    ylim([y_lim]);
end
set(gca,'linewidth',axis_width)
set(gca, 'TickDir', 'out')
title(highlight);

%set(gca, 'xtick', [])


%% PCA analysis of cancer data - 3D

data = ccleICT_a_hemato.values;
cells = ccleICT_a_hemato.cells;
genes = ccleICT_a_hemato.genes;
annot = annot_hemato;
highlight = ref;
highlight_ind = 5;      


color_trip = [0 0 0;
    200 200 200];
x_lim = 'N/A';
y_lim = 'N/A';
axis_width = 3;
barline_width = 3;
markersize = 500;
fig_size = [10 10 800 700];


tic

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
    scatter3(score(:,1),score(:,2),score(:,3),markersize,colors(1,:),'.');
    set(gca,'linewidth',axis_width)
    set(gca, 'TickDir', 'out')    
    
    figure
    scatter(score(:,1), score(:,2),markersize,colors(1,:),'.');
    set(gca,'linewidth',axis_width)
    set(gca, 'TickDir', 'out')    
    
    figure
    scatter(score(:,1),score(:,3),markersize,colors(1,:),'.');
    set(gca,'linewidth',axis_width)
    set(gca, 'TickDir', 'out')    
    
    figure
    scatter(score(:,2),score(:,3),markersize,colors(1,:),'.');
    set(gca,'linewidth',axis_width)
    set(gca, 'TickDir', 'out')
else
    
    % Finds indices of cell lines corresponding to the query
    group = [];
    for i = [1:length(highlight)]
        qcell = highlight{i};
        
        for j = [1:length(annot(:,1))]
            val = annot{j,highlight_ind};

            present = strfind(lower(string(val)), lower(string(qcell)));
            if not(isempty(present))
                group = [group j];
                continue
            end
        end
        
    end
    
    score1x = score(group,1);
    score1y = score(group,2);
    score1z = score(group,3);
    score(group,:) = [];
    
    scatter3(score(:,1), score(:,2), score(:,3), markersize, colors(2,:), '.');
    scatter3(score1x, score1y, score1z, markersize, colors(1,:), '.');
    set(gca,'linewidth',axis_width)
    set(gca, 'TickDir', 'out')
    title(highlight)
    
    figure
    hold
    scatter(score(:,1), score(:,2), markersize, colors(2,:), '.');
    scatter(score1x, score1y, markersize, colors(1,:), '.');
    set(gca,'linewidth',axis_width)
    set(gca, 'TickDir', 'out')
    xlabel('PC1')
    ylabel('PC2') 
    
    figure
    hold
    scatter(score(:,1), score(:,3), markersize, colors(2,:), '.');  
    scatter(score1x, score1z, markersize, colors(1,:), '.');
    set(gca,'linewidth',axis_width)
    set(gca, 'TickDir', 'out')
    xlabel('PC1')
    ylabel('PC3') 
    
    figure
    hold
    scatter(score(:,2), score(:,3), markersize, colors(2,:), '.');    
    scatter(score1y, score1z, markersize, colors(1,:), '.');
    set(gca,'linewidth',axis_width)
    set(gca, 'TickDir', 'out')
    xlabel('PC2')
    ylabel('PC3') 
end

if isnumeric(x_lim)
    xlim([x_lim]);
end
if isnumeric(y_lim)
    ylim([y_lim]);
end

%set(gca, 'xtick', [])

toc


%% PCA analysis of cancer data - 2D, multicolored

data = ccleICT_a_hemato.values;
cells = ccleICT_a_hemato.cells;
genes = ccleICT_a_hemato.genes;
annot = annot_hemato;
highlight = ref;
highlight_ind = 5;      


color_trip = [234 92 239;
    101 92 239;
    92 205 239;
    92 239 155;
    239 121 92;
    255 0 0;
    18 137 50;
    204 204 0;
    200 200 200];
x_lim = 'N/A';
y_lim = 'N/A';
axis_width = 3;
barline_width = 3;
markersize = 500
fig_size = [10 10 800 700];

tic

% PCA analysis (to acquire 'score')
[coeff, score, latent] = pca(data');

% Plots PCA analyses
figure('Renderer', 'painters', 'Position', fig_size);
hold
xlabel('PC1')
ylabel('PC2') 
zlabel('PC3')
colors = color_trip/255;

for i = [1:length(highlight(1,:))]
    grouping = highlight(:,i);
    
    % Finds indices of cell lines corresponding to the query
    group = [];
    for j = [1:length(grouping)]
        qcell = grouping{j};
        
        if isempty(qcell)
            break
        end
        
        for k = [1:length(annot(:,1))]
            val = annot{k,highlight_ind};
            
            if lower(string(val)) == lower(string(qcell))
                group = [group k];
            end
%             present = strfind(lower(string(val)), lower(string(qcell)));
%             if not(isempty(present))
%                 group = [group k];
%                 continue
%             end
        end

    end

    score1x = score(group,1);
    score1y = score(group,2);
    score1z = score(group,3);

    scatter3(score1x, score1y, score1z, markersize, colors(i,:), '.');
%     scatter(score1x, score1y, markersize, colors(i,:), '.');
%     xlabel('PC1')
%     ylabel('PC2') 
%     scatter(score1x, score1z, markersize, colors(i,:), '.');
%     xlabel('PC1')
%     ylabel('PC3') 
%     scatter(score1y, score1z, markersize, colors(i,:), '.');
%     xlabel('PC2')
%     ylabel('PC3') 
end
toc


if isnumeric(x_lim)
    xlim([x_lim]);
end
if isnumeric(y_lim)
    ylim([y_lim]);
end
set(gca,'linewidth',axis_width)
set(gca, 'TickDir', 'out')

%set(gca, 'xtick', [])


%% t-SNE analysis of cancer data

data = ccleICT_a_hemato.values;
cells = ccleICT_a_hemato.cells;
genes = ccleICT_a_hemato.genes;
annot = annot_hemato;
highlight = 'n/a';
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

% PCA analysis (to acquire 'score')
score = tsne(data');

% Plots PCA analyses
figure('Renderer', 'painters', 'Position', fig_size);
hold
colors = color_trip/255;

if lower(string(highlight)) == 'n/a'
    scatter(score(:,1),score(:,2),markersize,colors(1,:),'.');
    group1 = [annot(:,1) annot(:,2) annot(:,4) annot(:,5) num2cell(score(:,1)) num2cell(score(:,2))];
else
    
    % Finds indices of cell lines corresponding to the query
    group = [];
    group_cells = {};
    for i = [1:length(annot(:,1))]
        val = annot{i,highlight_ind};
        
        present = strfind(lower(string(val)), lower(string(highlight)));
        if not(isempty(present))
            group = [group i];
            group_cells = [group_cells; annot{i,2}];
            continue
        end
    end
    
    score1x = score(group,1);
    score1y = score(group,2);
    
    score(group,:) = [];
    annot(group,:) = [];
    group1_cells = annot(:,2);
    
    scatter(score(:,1), score(:,2), markersize, colors(2,:), '.');
    scatter(score1x, score1y, markersize, colors(1,:), '.');
       
    group1 = [group_cells num2cell(score1x(:,1)) num2cell(score1y(:,1))];
    group2 = [group1_cells num2cell(score(:,1)) num2cell(score(:,2))];
end
toc


if isnumeric(x_lim)
    xlim([x_lim]);
end
if isnumeric(y_lim)
    ylim([y_lim]);
end
set(gca,'linewidth',axis_width)
set(gca, 'TickDir', 'out')
title(highlight);

%set(gca, 'xtick', [])

clearvars annot axis_width barline_width cells color_trip colors data fig_size
clearvars genes group group1_cells group_cells highlight 
clearvars highlight_ind i markersize present score score1x score1y val x_lim
clearvars y_lim

%% Highlights cell-types of interest (based on a preformed t-SNE)

data = group1;
highlight = {'b-cell'};
highlight_ind = 4;      


color_trip = [0 0 0;
    200 200 200];
x_lim = 'N/A';
y_lim = 'N/A';
axis_width = 3;
barline_width = 3;
markersize = 500
fig_size = [10 10 800 700];


figure('Renderer', 'painters', 'Position', fig_size);
hold
colors = color_trip/255;

score = cell2mat(data(:,[5,6]));
hscore = [];
exclude = [];
for i = [1:length(highlight(:,1))]
    hterm = lower(string(highlight{i,1}));
    
    for j = [1:length(data(:,1))]
        dterm = lower(string(data{j,highlight_ind}));
        
        present = strfind(dterm, hterm);
        if isempty(present)
            continue
        end
        
        hscore = [hscore; data{j,5} data{j,6}];
        exclude = [exclude; j];
    end
end
score(exclude,:) = [];

scatter(score(:,1), score(:,2), markersize, colors(2,:), '.');
scatter(hscore(:,1), hscore(:,2), markersize, colors(1,:), '.');
    


if isnumeric(x_lim)
    xlim([x_lim]);
end
if isnumeric(y_lim)
    ylim([y_lim]);
end
set(gca,'linewidth',axis_width)
set(gca, 'TickDir', 'out')
title(highlight);

clearvars axis_width barline_width color_trip colors data dterm exclude
clearvars fig_size highlight highlight_ind hscore hterm i j markersize present
clearvars score x_lim y_lim


%% Acquires genes sorted according to PC# coefficients

genes = ccleICT_a.genes;
gcoeff = coeff;
pc_num = 2;

data = [genes num2cell(gcoeff(:,pc_num))];
output = sortrows(data,2,'desc');

clear genes gcoeff pc_num data


%% Plot genes that result in the highest variance 

variances = latent;
figure
hold


