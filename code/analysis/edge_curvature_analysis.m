%% Test correlation between behavior and curvature

% Load the curvature matrixes
% load(Load the curvature matrixes);

% Load behavioral data
% behav_data = csvread(Load behavioral data);

column_labels = {'pddbi_change', 'vabs_change', 'eow_change', 'cgi_6'};
vabs_scores = behav_data(:, 2);
eow_scores = behav_data(:, 3);
cgi_scores = behav_data(:, 4);

% Load Atlas
atlas_labels = load_atlas();

% Difference = TP2-TP1
difference_mat = zeros(size(curvature_mat_TP1));
for n = 1:19
    TP1 = curvature_mat_TP1(:, :, n);
    TP2 = curvature_mat_TP2(:, :, n);
    difference = TP2 ./ TP1;
    difference_mat(:, :, n) = difference;
end

count = 0;
region_pair = [];
correlation_mat = zeros(83, 83, 3);
pval_mat = zeros(83, 83, 3);

n = 1;
for r=1:83
    for c=1:83
        differences = squeeze(difference_mat(r, c, :));
        [rho_vabs, pval_vabs] = corr(differences, vabs_scores, 'Type', 'Spearman');
        [rho_eow, pval_eow] = corr(differences, eow_scores, 'Type', 'Spearman');
        [rho_cgi, pval_cgi] = corr(differences, cgi_scores, 'Type','Spearman');
        
        pval_vec = [pval_vabs, pval_eow, pval_cgi];
        pval_mat(r, c, :) = pval_vec;
        
        pval_vec = pval_vec < 0.01;
        rho_vec = [rho_vabs, rho_eow, rho_cgi];
        correlation_mat(r, c, :) = rho_vec;
        if sum(pval_vec) >= 2
            count = count + 1;
            region_pair(n, 1) = r;
            region_pair(n, 2) = c;
            n = n + 1;
        end
    end
end


foo = sort(region_pair, 2);
region_pair = unique(foo, 'rows');


pairs_list = {};
corr_list = {};
pval_list = {};

for n=1:size(region_pair, 1)
    edge = [atlas_labels{region_pair(n, 1)}, ' : ', atlas_labels{region_pair(n, 2)}];
    pairs_list{n} = edge;
    corr_list{n} = squeeze(correlation_mat(region_pair(n, 1), region_pair(n, 2), :));
    pval_list{n} = squeeze(pval_mat(region_pair(n, 1), region_pair(n, 2), :));
    
    differences = squeeze(difference_mat(region_pair(n, 1), region_pair(n, 2), :));
    
    disp(edge); %vabs eow cgi
    disp(strcat('Change in curvature mean: ', num2str(nanmean(differences)), ' std: ', num2str(nanstd(differences))));
    disp(strcat('corr-vabs: ', num2str(corr_list{n}(1)), ' corr-eow: ', ...
        num2str(corr_list{n}(2)), ' corr_cgi: ', num2str(corr_list{n}(3))));
    disp(strcat('pval-vabs: ', num2str(pval_list{n}(1)), ' pval-eow: ', ...
        num2str(pval_list{n}(2)), ' pval_cgi: ', num2str(pval_list{n}(3))));
    disp(' ');
end

corr_list = corr_list';
pval_list = pval_list';
pairs_list = pairs_list';


% Check FDR correction
At = pval_mat(:, :, 1).';
m  = tril(true(size(At)));
v  = At(m).';
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(v,.05,'pdep','yes');

At = pval_mat(:, :, 2).';
m  = tril(true(size(At)));
v  = At(m).';
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(v,.05,'pdep','yes');

At = pval_mat(:, :, 3).';
m  = tril(true(size(At)));
v  = At(m).';
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(v,.05,'pdep','yes');









