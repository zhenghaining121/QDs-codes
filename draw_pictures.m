%%
%This code is for drawing pictures of 1D 2D and pleateu steps
%code writer:JF 
%Time:20190724
%notice
%variable parameter:
%num_file\nmu_sampling\conduct_up_limit\conduct_down_limit\bins

%%
%load data
%parameters
clear
filename_str = 'goodtrace-trace-'; 
num_file = 3442;  
fileclass = '.txt'; 
num_sampling = 5000;        
num_count = 1;
logs_summary = zeros(num_file, 2);
logs_summary(:, 1) = 1 : num_file; 
conduct_up_limit = 0.8;
conduct_down_limit = -7;

Matrix_Condu = zeros(20000, num_sampling);
Matrix_Dist = zeros(20000, num_sampling);


%%
tic

for i = 1 : num_file 
    log_count = 0;
    filename = strcat(filename_str, num2str(i), fileclass);
    if (exist(filename) ~= 0)   
        original = load(filename);
        [m, n] = size(original);     
        Origal_maxIndex = m; 
        
        if Origal_maxIndex >num_sampling
            down = min(original(:, 2));
            if down <-7
                conduct_up_limit_index = find(original(:, 2) > conduct_up_limit);
                conduct_down_limit_index = find(original(:, 2) < conduct_down_limit);
                data = original(max(conduct_up_limit_index) : min(conduct_down_limit_index), :);
                maxIndex = size(data, 1);  
                if maxIndex < num_sampling
                    continue;
                end
    
                % cut data store
                ind = floor(linspace(1, maxIndex, num_sampling));
                data_sample = data(ind, :);
                Matrix_Condu(num_count, :) = data_sample(1:num_sampling, 2);
                Matrix_Dist(num_count, :) = data_sample(1:num_sampling, 1); 
        
                num_count = num_count + 1; 
                log_count = log_count + 1;  
    
                logs_summary(i, 2) = log_count;    
            end
        end
    end
end

zero_index = find(all(Matrix_Condu==0, 2));
Matrix_Condu(zero_index, :) = [];
Matrix_Dist(zero_index, :) = [];

TEST_condu = Matrix_Condu;
TEST_dist = Matrix_Dist;

[m1,n1] = size(TEST_condu);
%%
figure(1)
plotcloud_cond_max = 0.8;
plotcloud_cond_min = -6;
plotcloud_dist_max = 2;
plotcloud_dist_min = -0.5;
bins = 300;
conductance_all = reshape(TEST_condu,m1*n1,1);
dist_all = reshape(TEST_dist,m1*n1,1);
overlay_single = [dist_all,conductance_all];
[M,~,~] = plotCloud(overlay_single,plotcloud_cond_max,plotcloud_cond_min,plotcloud_dist_max,plotcloud_dist_min,bins);%You can change the range
imagesc(M,[0,100]);
title('figure 1');
xlabel('distance');
ylabel('conductance');
colorbar;
axis on

figure(2)
Condu_draw = reshape(TEST_condu,m1*n1,1);
histo_bins_condu = 500;
[gap_condu,dist_condu] = Gauss_draw_bins(Condu_draw,histo_bins_condu);
histogram(Condu_draw,gap_condu);
title('figure 2');
xlabel('conductacne(¡÷G)/log(G/G0)');
ylabel('counts');

figure(3)
Origin_Dist_2 = zeros(m1,1);
condu_downlimit = -5.4;
for i = 1:m1
    for j = 1:(n1-1)
        if (TEST_condu(i,j)>condu_downlimit)&&(TEST_condu(i,j+1)<condu_downlimit)
            Origin_Dist_2(i,1) = TEST_dist(i,j);
        end
    end
end
histo_bins_dist = 50;
[gap_dist,dist_dist] = Gauss_draw_fixbins(Origin_Dist_2,histo_bins_dist);%Êä³öµÄºá×Ý×ø±ê
histogram(Origin_Dist_2,gap_dist);
axis([0 2 0 100]);
title('figure 3');
xlabel('distance/nm');
ylabel('counts');






