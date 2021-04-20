
%% Generate outage file
clear all;
close all;
clc;

%%
save_path = '/home/tony/Desktop/';

start_time = 9685;
end_time = 10558;
save_file = [save_path, 'outages.outage'];
disp(save_file);
outages = GenerateOutages(start_time, end_time, 60, 180, 2, save_file);
