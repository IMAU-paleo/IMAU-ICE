clc
clear all
close all

henk = dir;

find_replace = {'_IceOcean0_','_IceOcean2ra_'};

for i = 1:length(henk)
  filename_old = henk(i).name;
  filename_new = [];
  if length(filename_old)<length(find_replace{1}); continue; end
  l = length(find_replace{1});
  for ci = 1:length(filename_old)-l+1
    if strcmpi(filename_old(ci:ci+l-1),find_replace{1})
      filename_new = [filename_old(1:ci-1) find_replace{2} filename_old(ci+l:end)];
      break
    end
  end
  if ~isempty(filename_new)
    disp(['Copying ' filename_old ' to ' filename_new])
    copyfile(filename_old,filename_new);
  end
end