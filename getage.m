function age = getage(H,cent)

xmin = .2; xmax = .9; ymin = .15; ymax = .4;
           if length(H) > 2000
               im_dim = 2048;
           else
               im_dim = 1024;
           end
           lims =  round([xmin,xmax,ymin,ymax]*im_dim); 
           nd = NucDensity(cent,lims,0);
     
          age = floor(4.8 + log2(nd));