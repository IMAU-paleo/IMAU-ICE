function Plot_MISMIPplus_ISOMIPplus_MISOMIP1_Venn_diagram

close all

wa = 700;
ha = 440;

wf = 25+wa+25;
hf = 25+ha+25;

H.Fig = figure('position',[300,300,wf,hf],'color','w');
H.Ax  = axes('parent',H.Fig,'units','pixels','position',[25,25,wa,ha],'xlim',[0,wa],'ylim',[0,ha],'xtick',[],'ytick',[]);
H.Ax.XAxis.Visible = 'off';
H.Ax.YAxis.Visible = 'off';

r_in  = 160;
r_out = 220;
nc    = 20;
alpha_in = 0;
alpha_out = 0.5;

r_inv  = linspace(r_in,r_out,nc+1); r_inv  = r_inv(  1:end-1);
r_outv = linspace(r_in,r_out,nc+1); r_outv = r_outv( 2:end  );
alphav = linspace(alpha_in,alpha_out,nc);

color_ice  = [0.6,0.6,0.7];
centre_ice = [230,220];

color_sea  = [0.4,0.4,1.0];
centre_sea = [470,220];

for j = 1:nc
  donut_patch( H.Ax,centre_ice(1),centre_ice(2),r_inv(j),r_outv(j),color_ice,alphav(j));
  donut_patch( H.Ax,centre_sea(1),centre_sea(2),r_inv(j),r_outv(j),color_sea,alphav(j));
end

donut_patch( H.Ax,centre_ice(1),centre_ice(2),0,r_out,color_ice,0.1);
donut_patch( H.Ax,centre_sea(1),centre_sea(2),0,r_out,color_sea,0.1);

donut_patch( H.Ax,centre_ice(1),centre_ice(2),r_out*0.98,r_out,color_ice*0.8,1);
donut_patch( H.Ax,centre_sea(1),centre_sea(2),r_out*0.98,r_out,color_sea*0.8,1);

text(H.Ax,centre_ice(1)-40,centre_ice(2)+r_out-50,'Ice'  ,'fontsize',36,'color',color_ice*0.6);
text(H.Ax,centre_sea(1)-30,centre_sea(2)+r_out-50,'Ocean','fontsize',36,'color',color_sea*0.6);

text(H.Ax,centre_ice(1)-170,centre_ice(2),'MISMIP+','fontsize',36,'color',color_ice*0.3);
text(H.Ax,0.5*(centre_ice(1)+centre_sea(1))-85.5,centre_ice(2),'MISOMIP1','fontsize',36,'color',0.5*(color_ice+color_sea)*0.3);
text(H.Ax,centre_sea(1)+25,centre_ice(2),'ISOMIP+','fontsize',36,'color',color_sea*0.3);

  function p = donut_patch( ax,xc,yc,r_in,r_out,color,alpha)
    
    n = 100;
    thetav = linspace(0,2*pi,n);
    
    nV = 2*n;
    V  = zeros(nV,2);
    F  = zeros(n,4);
    
    vi = 0;
    fi = 0;
    for i = 1:n
      % inner vertex
      x_in = xc + r_in * sin(thetav(i));
      y_in = yc + r_in * cos(thetav(i));
      vi = vi + 1;
      V(vi,:) = [x_in,y_in];
      % outer vertex
      x_out = xc + r_out * sin(thetav(i));
      y_out = yc + r_out * cos(thetav(i));
      vi = vi + 1;
      V(vi,:) = [x_out,y_out];
      % face
      vi1 = vi-1;
      vi2 = vi+1;
      vi3 = vi+2;
      vi4 = vi;
      if (vi2>nV); vi2 = vi2-nV; end
      if (vi3>nV); vi3 = vi3-nV; end
      fi = fi+1;
      F(fi,:) = [vi1,vi2,vi3,vi4];
    end
    
    p = patch('parent',ax,'vertices',V,'faces',F,'edgecolor','none','facecolor',color,'facealpha',alpha);
    
  end

end