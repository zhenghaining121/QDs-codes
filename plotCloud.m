function [Mat, x, y] = plotCloud(data, maxConduct, minConduct, maxDist, minDist,bins)

n_bins = bins;
nanData = -n_bins;
s_C_bins = (maxConduct - minConduct)/(n_bins);
s_D_bins = (maxDist - minDist)/n_bins;
Conduct_index = linspace(maxConduct, minConduct, n_bins);
Dist_index = [nanData, linspace(minDist, maxDist, n_bins)];
% Conduct_index = roundn(Conduct_index, -3);
% Dist_index = roundn(Dist_index, -3);
Matrix = zeros(n_bins, n_bins);
Matrix = [Conduct_index', Matrix];
Matrix = [Dist_index; Matrix];
% data = roundn(data, -3);
[m, n] = size(data);

for i = 1: m
    
    x_ind = find((data(i, 1) < Matrix(1, :)+s_D_bins/2) & (data(i, 1) >= Matrix(1, :)-s_D_bins/2));
    y_ind = find((data(i, 2) < Matrix(:, 1)+s_C_bins/2) & (data(i, 2) >= Matrix(:, 1)-s_C_bins/2));
    
    Matrix(y_ind, x_ind) = Matrix(y_ind, x_ind) + 1;
    
    
end

Mat = Matrix(2:end, 2:end);
x = Matrix(1, 2:end);
y = Matrix(2:end, 1);

