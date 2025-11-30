E = csvread('example2.dat');

col1 = E(:,1);
col2 = E(:,2);

max_ids = max(max(col1,col2));
As= sparse(col1, col2, 1, max_ids, max_ids); 
A = full(As);

D = diag(sum(A,2));
L = D - A;

[v, D] = eig(L);
[vals_sorted, order] = sort(diag(D));
fiedler = v(:, order(2))

[f_sorted, idx_sorted] = sort(fiedler);
figure;
plot(f_sorted, '-o', 'MarkerSize',6);
xlabel('Node (sorted by Fiedler value)');
ylabel('Fiedler value');
title('Sorted Fiedler Vector — community indication');
grid on;