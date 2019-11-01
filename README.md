# QDs-codes
This code is for selecting jump traces and drawing pictures for QDs

This repository contains 7 MATLAB codes for selecting jump traces for QDs data and drawing 1D conductance histograms and 2D conductance-displacement cloud maps.

The main code is select_jumptraces_final_version.m and draw_pictures.m

The parameters introduction of MATLAB code select_jumptraces_final_version.m
%parameters introduction
%num_file=the total number of files; 
%num_sampling = the number of sampling points for every trace; 
%conduct_up_limit=the uplimit value of conductance for single trace; 
%conduct_down_limit = the downlimit value of conductance for single trace;
% four_start = the threhold value of jump conductance
% judge_number = the number of jumping points
% four_c_start = the uplimit value of conductance for jumping trace; 
% four_c_end = the downlimit value of conductance for jumping trace; 
% length = the display points in picture
% forward_length = the relative zero-displacement position before jumping point 
