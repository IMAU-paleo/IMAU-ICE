clc
clear all
close all

% Compress and bundle all ISMIP results folders

%% Create folders to hold the compressed files
for mi = 1:8
  foldername = ['PROTECT_projections_IMAUICE' num2str( mi)];
  if exist( foldername,'dir') == 7
    henk = dir( foldername);
    for i = 3: length( henk)
      delete( [foldername '/' henk( i).name])
    end
    rmdir( foldername)
  end
  mkdir( foldername)
end

%% Find, compress, and move all ISMIP results folders
for mi = 1: 8
  foldername_dst = ['PROTECT_projections_IMAUICE' num2str( mi)];
  
  henk = dir;
  for i = 1: length( henk)
    if henk( i).isdir && startsWith( henk( i).name,['IMAUICE' num2str( mi) '_'])
      piet = dir( henk( i).name);
      for j = 1: length( piet)
        if piet( j).isdir && startsWith( piet( j).name,'GIS_IMAU_IMAUICE')
          % This is the ISMIP output folder. Compress it.
          filenames_src = {};
          jan = dir( [henk( i).name '/' piet( j).name]);
          for k = 1: length( jan)
            filename = [henk( i).name '/' piet( j).name '/' jan( k).name];
            if contains( filename,'.nc')
              filenames_src{ end+1} = filename;
            end
          end
          filename_dst = [foldername_dst '/' piet( j).name '.tgz'];
          disp(['Compressing ' henk( i).name '...'])
          tar( filename_dst, filenames_src);
        end
      end
    end
  end
end