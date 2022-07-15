clc
clear all
close all

ensemble_name = 'MISMIPplus';

henk = dir;
foldernames = {};
for i = 1:length(henk)
  if length(henk(i).name) > length(ensemble_name)
    if strcmpi(henk(i).name(1:length(ensemble_name)),ensemble_name) && henk(i).isdir
      foldernames{end+1} = henk(i).name;
    end
  end
end

for fi = 1:length(foldernames)
  generate_MISMIPplus_output( foldernames{fi})
end