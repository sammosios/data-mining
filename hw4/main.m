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

%Number of 0 eigenvalues for graph analysis

GL = D - A; % construct normal Laplacian

e = eig(GL); % get all eigenvalues
num_zero_eigenvalues = sum(abs(e) < 1e-10); % count eigenvalues close to 0
disp(num_zero_eigenvalues)

%Sorted Fiedler Observations

[v, D] = eig(GL);
[vals_sorted, order] = sort(diag(D));
fiedler = v(:, order(2));

[f_sorted, idx_sorted] = sort(fiedler);
figure;
plot(f_sorted, '-o', 'MarkerSize',6);
xlabel('Node (sorted by Fiedler value)');
ylabel('Fiedler value');
title('Sorted Fiedler Vector — community indication');
grid on;

k = 4; % Cluster Number

[X, eigenValues] = eigs(L, k, 'largestabs'); %Extracting k largest eigenvalues

N = size(X, 1);

Y = zeros(N, k);

for i = 1:N
    row_norm = sqrt(sum(X(i, :).^2));   % compute the norm of each row, which is the sqare root of the sum of squares of the row values
    Y(i, :) = X(i, :) / row_norm;       % divide each element by row norm
end

[labels] = kmeans(Y, k, 'Replicates', 10); %Run k means on the normalized matrix

%Plotting Details

% Project to first two dimensions
X_plot = X(:,1:2);

% Colormap
colors = lines(k);

figure('Color', 'w');
ax = gca;      
ax.Color = 'w';  
ax.XColor = 'k'; 
ax.YColor = 'k'; 
hold on;

%Plot edges
for i = 1:length(col1)
    n1 = col1(i);
    n2 = col2(i);
    plot([X_plot(n1,1), X_plot(n2,1)], [X_plot(n1,2), X_plot(n2,2)], 'k-', 'LineWidth', 0.5); 
end

% Plot nodes
h_points = gscatter(X_plot(:,1), X_plot(:,2), labels, colors, 'o', 8);

xlabel('First eigenvector dimension');
ylabel('Second eigenvector dimension');
title('Spectral Clustering');
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