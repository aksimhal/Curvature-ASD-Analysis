%% Scalar curvature analysis
% Brain regions that correlate with behavioral tests

% Load the curvature matrixes
% load(Load the curvature matrixes);

% Load Atlas
atlas_labels = load_atlas();

% Compute scalar curvature
scalar_k_mat_TP1 = zeros(83, 19);
for z = 1:19
    for node_n = 1:83
        node_edges_vec = curvature_mat_TP1(node_n, :, z);
        dx = sum(node_edges_vec);
        mu_vec = node_edges_vec / dx;
        k_mu_vec = node_edges_vec.*mu_vec;
        scalar_curvature = sum(k_mu_vec);
        scalar_k_mat_TP1(node_n, z) = scalar_curvature;
    end
end

% Compute scalar curvature
scalar_k_mat_TP2 = zeros(83, 19);
for z = 1:19
    for node_n = 1:83
        node_edges_vec = curvature_mat_TP2(node_n, :, z);
        dx = sum(node_edges_vec); 
        mu_vec = node_edges_vec / dx;
        k_mu_vec = node_edges_vec.*mu_vec;
        scalar_curvature = sum(k_mu_vec);
        scalar_k_mat_TP2(node_n, z) = scalar_curvature;
    end
end

% Difference = TP2/TP1
difference_nodes = scalar_k_mat_TP2 ./ scalar_k_mat_TP1;

% Load behavioral data
% behav_data = csvread(Load behavioral data);

column_labels = {'pddbi_change', 'vabs_change', 'eow_change', 'cgi_6'};
vabs_scores = behav_data(:, 2);
eow_scores = behav_data(:, 3);
cgi_scores = behav_data(:, 4);

count = 0;
region_pair = [];
correlation_mat = zeros(83, 3);
pval_mat = zeros(83, 3);
n = 1;

for node_n=1:83
    differences = squeeze(difference_nodes(node_n, :));
    differences = differences(:);
    [rho_vabs, pval_vabs] = corr(differences, vabs_scores, 'Type', 'Spearman', 'rows','complete');
    [rho_eow, pval_eow] = corr(differences, eow_scores, 'Type', 'Spearman', 'rows','complete');
    [rho_cgi, pval_cgi] = corr(differences, cgi_scores, 'Type','Spearman', 'rows','complete');
    
    pval_vec = [pval_vabs, pval_eow, pval_cgi];
    pval_mat(node_n, :) = pval_vec;
    
    pval_vec = pval_vec < 0.05;
    rho_vec = [rho_vabs, rho_eow, rho_cgi];
    correlation_mat(node_n, :) = rho_vec;
    if sum(pval_vec) >= 2
        count = count + 1;
        region_pair(n, 1) = node_n;
        n = n + 1;
    end
end

disp('Change in VABS Score'); 
disp(mean(vabs_scores)); 
disp(std(vabs_scores)); 

disp('Change in EOW'); 
disp(mean(eow_scores)); 
disp(std(eow_scores)); 

disp('Changes in CGI'); 
disp(mean(cgi_scores));
disp(std(cgi_scores)); 


for n = 1:length(region_pair)
    disp(atlas_labels{region_pair(n)});
    noderatio = difference_nodes(region_pair(n), :);
    disp(nanmean(noderatio))
    disp(nanstd(noderatio));
end


% Test to see if the results survive FDR correction
alpha = 0.05;
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pval_mat(:, 1),alpha,'pdep','yes');
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pval_mat(:, 2),alpha,'pdep','yes');
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pval_mat(:, 3 ),alpha,'pdep','yes');

