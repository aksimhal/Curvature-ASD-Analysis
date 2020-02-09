%% Generate Curvature Matrixes 

% Load the connectivity matrixes
% load(connectivity matrixes);
% load(white matter volume for normalization);

% Time point 1 
curvature_mat_TP1 = zeros(size(total_mat_vox_TP1)); 
for n = 1:size(total_mat_vox_TP1, 3)
    A = total_mat_vox_TP1(:, :, n);
    volume = WM_matrix(n, 1);
    A = A/volume; 
    K = calculate_curvature(A);
    curvature_mat_TP1(:, :, n) = K; 
    disp(n);
end

% Time point 2
curvature_mat_TP2 = zeros(size(total_mat_vox_TP2)); 
for n = 1:size(total_mat_vox_TP2, 3)
    A = total_mat_vox_TP2(:, :, n);
    volume = WM_matrix(n, 2); 
    A = A/volume; 
    K = calculate_curvature(A);
    curvature_mat_TP2(:, :, n) = K; 
    disp(n);
end

save('curvature_ASD_hop.mat', 'curvature_mat_TP1', 'curvature_mat_TP2');

