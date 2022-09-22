function map = redbluemap(n)
  ncolors=3;
  c=flipud([...
    255    0    0;
    255   255   255;...
    0     0     255]);
  pp=1:(n-1)/(ncolors-1):n;
  r=interp1(pp,c(:,1),1:n);
  g=interp1(pp,c(:,2),1:n);
  b=interp1(pp,c(:,3),1:n);
  map=[r' g' b']/255;
end