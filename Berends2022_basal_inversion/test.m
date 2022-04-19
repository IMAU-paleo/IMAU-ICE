clc
clear all
close all

filename = 'exp_II_target_5km/restart_ANT.nc';

Hi_noneq = ncread( filename,'Hi',[1,1, 3],[Inf,Inf,1]);
Hi_eq    = ncread( filename,'Hi',[1,1,21],[Inf,Inf,1]);

imagesc3( Hi_noneq ./ Hi_eq )