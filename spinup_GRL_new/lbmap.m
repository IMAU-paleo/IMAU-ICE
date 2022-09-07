function map = lbmap(n,scheme)
  %valid schemes
  switch lower(scheme)
    case 'blue'
      baseMap = BlueMap;
    case 'bluegray'
      baseMap = BlueGrayMap;
    case 'brownblue'
      baseMap = BrownBlueMap;
    case 'redblue'
      baseMap = RedBlueMap;
    otherwise
      error(['Invalid scheme ' scheme])
  end
  idx1 = linspace(0,1,size(baseMap,1));
  idx2 = linspace(0,1,n);
  map = interp1(idx1,baseMap,idx2);

  function baseMap = BlueMap
  baseMap = [243 246 248;
             224 232 240;
             171 209 236;
             115 180 224;
              35 157 213;
               0 142 205;
               0 122 192]/255;
  end
  function baseMap = BlueGrayMap
  %DivergingBlueGray
  baseMap = [  0 170 227;
              53 196 238;
             133 212 234;
             190 230 242;
             217 224 230;
             146 191 170;
             109 122 129;
              65  79  81]/255;
  end
  function baseMap = BrownBlueMap
  baseMap = [144 100  44;
             187 120  54;
             225 146  65;
             248 184 139;
             244 218 200;
             241 244 245;
             207 226 240;
             160 190 225;
             109 153 206;
              70  99 174;
              24  79 162]/255;
  end
  function baseMap = RedBlueMap
  baseMap = [175  53  71;
             216  82  88;
             239 133 122;
             245 177 139;
             249 216 168;
             242 238 197;
             216 236 241;
             154 217 238;
              68 199 239;
               0 170 226;
               0 116 188]/255;
  end
end