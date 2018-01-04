clear all
close all
clc

%origin time
to = 5105.47;

tt = load('arrival_times.out');

tt_sum = tt + to;
save('-ascii','final_arrival_times.out','tt_sum')
