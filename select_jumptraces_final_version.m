%%
%This code is for selecting jumping conductance traces
%code writer:JF 
%Time:20190309 finalversion:20190405
%notice
%variable parameter: num_file\nmu_sampling\conduct_up_limit\conduct_down_limit


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
conduct_down_limit = -7.0;

Matrix_Condu = zeros(20000, num_sampling);
Matrix_Dist = zeros(20000, num_sampling);

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


file_number = num_count - 1;
length = 400;
forward_length = 150;
Origin_Condu = zeros(20000, num_sampling);
Origin_Dist = zeros(20000, num_sampling);
Select_Condu = zeros(20000, length);
Select_Dist = zeros(20000, length);
Satis_Condu = zeros(20000, length);
Satis_Dist = zeros(20000, length);
four_start = 0.4;
judge_number = 30;
four_c_start = -4.5;
four_c_end = -6;
Record_jumppoint = zeros(file_number, 1);
        

for i = 1:file_number
    for j = 1:(num_sampling - 2-judge_number)
        if ((TEST_condu(i,j+judge_number)-four_start)>TEST_condu(i,j))&&(TEST_condu(i,j+1)>TEST_condu(i,j))&&(TEST_condu(i,j)<four_c_start)&&(TEST_condu(i,j)>four_c_end)
            Origin_Condu(i,:) = TEST_condu(i,:);
            Origin_Dist(i,:) = TEST_dist(i,:);
            Select_Condu(i,:) = TEST_condu(i,j-forward_length:j+length-forward_length-1)-TEST_condu(i,j-forward_length);
            Select_Dist(i,:) = TEST_dist(i,j-forward_length:j+length-forward_length-1)-TEST_dist(i,j-forward_length); 
            Satis_Condu(i,:) = TEST_condu(i,j-forward_length:j+length-forward_length-1)-TEST_condu(i,j);
            Satis_Dist(i,:) = TEST_dist(i,j-forward_length:j+length-forward_length-1)-TEST_dist(i,j); 
            Record_jumppoint(i,1) = j;
        end
        continue;
    end
end           
         
zero_index = find(all(Select_Condu==0, 2));
zero_index2 = find(all(Record_jumppoint==0,2));
Origin_Condu(zero_index, :) = [];
Origin_Dist(zero_index, :) = [];
Select_Condu(zero_index, :) = [];
Select_Dist(zero_index, :) = []; 
Satis_Condu(zero_index, :) = [];
Satis_Dist(zero_index, :) = [];  
Record_jumppoint(zero_index2,:) = [];
[sm,sn] = size(Select_Condu);
[sm1,sn1] = size(Origin_Condu);

jump_conductance = zeros(sm,1);
jump_distance = zeros(sm,1);
for i = 1:sm
    jump_conductance(i,1) = Origin_Condu(i,Record_jumppoint(i,1));
    jump_distance(i,1) = Origin_Dist(i,Record_jumppoint(i,1));
end

figure(1)
plotcloud_cond_max = 1.4;
plotcloud_cond_min = -1.4;
plotcloud_dist_max = 0.16;
plotcloud_dist_min = 0;
bins = 160;
conductance_all = reshape(Select_Condu,sm*length,1);
dist_all = reshape(Select_Dist,sm*length,1);
overlay_single = [dist_all,conductance_all];
[M,~,~] = plotCloud(overlay_single,plotcloud_cond_max,plotcloud_cond_min,plotcloud_dist_max,plotcloud_dist_min,bins);%You can change the range
%GAUSSIAN FIT
Gauss_peak = Gauss_peak_fit(M,bins);
imagesc(M,[0,100]);
title('figure 1');
xlabel('distance');
ylabel('conductance');
hold on
plot(1:bins,Gauss_peak,'Color',[0 0 1],'LineWidth',2.5);
colorbar;
axis off

figure(2)
plotcloud_cond_max2 = 1.4;
plotcloud_cond_min2 = -1.4;
plotcloud_dist_max2 = 0.17;
plotcloud_dist_min2 = -0.1;
bins = 200;
conductance_all2 = reshape(Satis_Condu,sm*length,1);
dist_all2 = reshape(Satis_Dist,sm*length,1);
overlay_single2 = [dist_all2,conductance_all2];
[M2,~,~] = plotCloud(overlay_single2,plotcloud_cond_max2,plotcloud_cond_min2,plotcloud_dist_max2,plotcloud_dist_min2,bins);%You can change the range
%GAUSSIAN FIT
Mx2 = zeros(bins, bins);
My2 = zeros(bins, bins);
Mymax2 = zeros(1, bins);
Gauss_peak2 = zeros(1, bins);
% Gauss_peak2 = Gauss_peak_fit(M2,bins)
imagesc(M2,[0,100]);
title('figure 2');
xlabel('distance');
ylabel('conductance');
% hold on
% plot(1:bins,Gauss_peak2,'Color',[0 0 1],'LineWidth',2.5);
colorbar;
axis off

figure(3)
high_jump_condu = zeros(sm,1);
low_jump_condu = zeros(sm,1);
offset = zeros(sm,1);
condu_histogram = zeros(sm,1);
point_high_jump = zeros(sm,1);
for i = 1:sm
    high_jump_condu(i,:) = max(Select_Condu(i,160:300));%[160/(500/200):300/(500/200)]= 48
    num_high_jump = find(Select_Condu(i,:)== max(Select_Condu(i,120:300)));
    point_high_jump(i,:) = max(num_high_jump);
    low_jump_condu(i,:) = min(Select_Condu(i,1:num_high_jump));
    offset(i,:) = Select_Condu(i,1) - low_jump_condu(i,:);
    condu_histogram(i,:) = high_jump_condu(i,:) + offset(i,:);
end

histo_bins = 70;
[histo_gap,histo_condu] = Gauss_draw_bins(condu_histogram,histo_bins);
  
histogram(condu_histogram,histo_bins);
title('figure 3');
xlabel('relative conductacne(¡÷G)/log(G/G0)');
ylabel('counts');
axis([0 2 0 190]);
a1 =      129.2  ;
b1 =      0.8125  ;
c1 =      0.4468  ;
Gauss_fit = a1.*exp(-((histo_gap-b1)./c1).^2);
hold on
plot(histo_gap,Gauss_fit,'Color',[0.070 0.637 0.922],'LineWidth',2.5);%x:histo_gap  y:histo_condu

figure(4)
max_jump_point = zeros(sm,1);
max_jump_condu = zeros(sm,1);
for i = 1:sm
    max_jump_point(i,1) = Record_jumppoint(i,1) + point_high_jump(i,1) - forward_length;
    max_jump_condu(i,1) = Origin_Condu(i,max_jump_point(i,1));
end
histo_bins2 = 35;
[histo_gap2, histo_condu2] = Gauss_draw_bins(max_jump_condu,histo_bins2);

histogram(max_jump_condu,histo_bins2);
title('figure 4');
xlabel('relative conductacne(¡÷G)/log(G/G0)');
ylabel('counts');
axis([-6 -2 0 250]);
a2 =      175.3  ;
b2 =      -4.318  ;
c2 =      0.6765  ;
Gauss_fit2 = a2.*exp(-((histo_gap2-b2)./c2).^2);
hold on
plot(histo_gap2,Gauss_fit2,'Color',[0.070 0.637 0.922],'LineWidth',2.5);%x:histo_gap2  y:histo_condu2

figure(5)
origin_jump_point = zeros(sm,1);
origin_jump_condu = zeros(sm,1);
for i = 1:sm
    origin_jump_point(i,1) = Record_jumppoint(i,1);
    origin_jump_condu(i,1) = Origin_Condu(i,Record_jumppoint(i,1));
end
histogram(origin_jump_condu);
    
toc
