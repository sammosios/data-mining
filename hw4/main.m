E = csvread('example1.dat');

col1 = E(:,1);
col2 = E(:,2);

col1 = E(:,1);
col2 = E(:,2);
max_ids = max(max(col1,col2));
As= sparse(col1, col2, 1, max_ids, max_ids); 
A = full(As);  %Created Adjacency Matrix (since we only have normal graph edges, this is equivalent to an affinity matrix)

D = diag(sum(A, 2)); %Diagonal Degree Matrix

D_inv_sqrt = diag(1 ./ sqrt(diag(D)));
L = D_inv_sqrt * A * D_inv_sqrt; %Normalized Laplacian

k = 2; % Cluster Number

[X, eigenValues] = eigs(L, k, 'largestabs'); %Extracting k largest eigenvalues

Y = zeros(N, k);

for i = 1:N
    row_norm = sqrt(sum(X(i, :).^2));   % compute the norm of each row, which is the sqare root of the sum of squares of the row values
    Y(i, :) = X(i, :) / row_norm;       % divide each element by row norm
end

[labels, C] = kmeans(Y, k, 'Replicates', 10); %Run k means on the normalized matrix

%Plotting Deatils

% Project to first two dimensions
X_plot = Y(:,1:2);
C_plot = C(:,1:2);

% Colormap
colors = lines(k);

figure;
hold on;

% Plot nodes
h_points = gscatter(X_plot(:,1), X_plot(:,2), labels, colors, 'o', 8);

% Plot centroids
h_centroids = gobjects(k,1);  % preallocate handles
for c = 1:k
    h_centroids(c) = plot(C_plot(c,1), C_plot(c,2), 'x', ...
        'MarkerSize',12, 'LineWidth',2, 'Color', colors(c,:));
end

% Create custom legend
legend_entries = cell(2*k,1);
for c = 1:k
    legend_entries{c} = sprintf('Cluster %d', c);          % points
    legend_entries{k+c} = sprintf('Centroid %d', c);       % centroids
end

% Combine handles for legend: first the scatter group handles, then centroid handles
all_handles = [h_points; h_centroids];

legend(all_handles, legend_entries, 'Location', 'best');

xlabel('First eigenvector dimension');
ylabel('Second eigenvector dimension');
title('Spectral Clustering with Centroids');
grid on;
hold off;

%Cluster Metrics

%Intra-Cluster Average Degree

avg_internal_degree = zeros(k, 1);

for c = 1:k
    nodes_in_c = find(labels == c);     
    A_sub = A(nodes_in_c, nodes_in_c);  
    internal_degrees = sum(A_sub, 2);   
    avg_internal_degree(c) = mean(internal_degrees);  
end

for c = 1:k
    fprintf('Cluster %d average internal degree: %.4f\n', c, avg_internal_degree(c));
end

%Inter-Cluster Average Degree

inter_cluster_degree = zeros(k,1);  

for c = 1:k
    nodes_in_c = find(labels == c);       
    nodes_out_c = find(labels ~= c);      
    A_out = A(nodes_in_c, nodes_out_c);
    out_degrees = sum(A_out, 2);
    inter_cluster_degree(c) = mean(out_degrees);
end

for c = 1:k
    fprintf('Cluster %d average inter-cluster degree: %.4f\n', c, inter_cluster_degree(c));
end

%Overall Averages
fprintf('Average internal cluster degree: %.4f\n', mean(avg_internal_degree));
fprintf('Average inter-cluster degree: %.4f\n', mean(inter_cluster_degree));