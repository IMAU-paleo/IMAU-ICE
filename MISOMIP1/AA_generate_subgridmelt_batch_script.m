clc
clear all
close all

ensemble_name = 'subgridmelt';

%% List all templates and variations

templates  = {};
variations = {};

henk = dir('config-files');
for i = 1:length(henk)
  varname = [ensemble_name '_var'];
  if length(henk(i).name) >= length(varname)
    if strcmpi(henk(i).name(1:length(varname)),varname)
      variations{end+1} = henk(i).name;
    end
  end
  templatename = [ensemble_name '_template'];
  if length(henk(i).name) >= length(templatename)
    if strcmpi(henk(i).name(1:length(templatename)),templatename)
      templates{end+1} = henk(i).name;
    end
  end
end

%% Generate the batch script

filename = ['AA_run_' ensemble_name '_batch.csh'];

if exist(filename,'file')
  delete(filename)
end

fid = fopen(filename,'w');

fprintf(fid,'#! /bin/csh -f\n');
fprintf(fid,'\n');
fprintf(fid,'cd ..\n');
fprintf(fid,'./compile_all.csh\n');
fprintf(fid,'\n');

for ti = 1:length(templates)
  for vi = 1:length(variations)
    fprintf(fid,['mpiexec -n 2 IMAU_ICE_program   ',...'
      'MISOMIP1/config-files/', templates{ti}, '   ', ...
      'MISOMIP1/config-files/', variations{vi},'\n']);
  end
  fprintf(fid,'\n');
end

fclose(fid);