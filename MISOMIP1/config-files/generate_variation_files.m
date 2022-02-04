function generate_variation_files( Dims, ensemble_name)

% Determine total number of ensemble members
n = 1;
for di = 1:length(Dims)
  n = n * length(Dims{di}.Opts);
end

% List all possible variations
allcombis = zeros( n,length(Dims));
combi = ones( 1,length(Dims));
for i = 1:n
  allcombis( i,:) = combi;
  combi(end) = combi(end)+1;
  for j = length(Dims):-1:2
    if combi(j) > length(Dims{j}.Opts)
      combi(j) = 1;
      combi(j-1) = combi(j-1)+1;
    end
  end
end

% Create all possible variations
for i = 1:n
  
  combi = allcombis(i,:);
  
  Variation.suffix = '';
  Variation.Vars = {};
  
  for di = 1:length(Dims)
    Opt = Dims{di}.Opts{combi(di)};
    Variation.suffix = [Variation.suffix, Opt.suffix];
    for vi = 1:length(Opt.Vars)
      Variation.Vars{end+1} = Opt.Vars{vi};
    end
  end
  
  generate_variation_file( Variation, ensemble_name)
  
end

  function generate_variation_file( Variation, ensemble_name)
    
    filename = [ensemble_name '_var' Variation.suffix];
    
    emptyline = '';
    for i = 1:40
      emptyline = [emptyline, ' '];
    end
    
    if exist(filename,'file')
      delete(filename)
    end
    
    fid = fopen(filename,'w');
    
    fprintf(fid,'&CONFIG\n');
    fprintf(fid,'\n');
    fprintf(fid,['fixed_output_dir_suffix_config          = ''', Variation.suffix, '''\n']);
    fprintf(fid,'\n');
    
    for vi = 1:length(Variation.Vars)
      strline = emptyline;
      strline(1:length(Variation.Vars{vi}.Name)) = Variation.Vars{vi}.Name;
      
      val = Variation.Vars{vi}.Value;
      valval = str2double(val);
      if isnan(valval)
        val = ['''' val ''''];
      end
      
      strline = [strline, ' = ', val '\n'];
      fprintf(fid,strline);
    end
    
    fprintf(fid,'\n');
    fprintf(fid,'/');
    fprintf(fid,'\n');
    
    fclose(fid);
    
  end

end