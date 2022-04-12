clc
clear all
close all

henk = dir;

for i=1:length(henk)
  if contains(henk(i).name,'exp_II') && contains(henk(i).name,'5_km') && henk(i).isdir
    filename = [henk(i).name '/restart_ANT.nc'];
    
    time = ncread(filename,'time');
    Hi   = ncread(filename,'Hi');
    Hs   = ncread(filename,'Hs');
    for ti = 1: length(time)
      Hi(1,:,ti) = Hi(1,1,ti);
      Hs(1,:,ti) = Hs(1,1,ti);
    end
    ncwrite(filename,'Hi',Hi);
    ncwrite(filename,'Hs',Hs);
  end
end