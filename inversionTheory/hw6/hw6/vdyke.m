function gz=vdyke(xobs, zobs, x, b, zt, zb, rho)
%
% function gz=vdyke(xobs, zobs, x, b, zt, zb, rho)
%
% function to calculate the gravity field(s) due to 
% a single vertical dyke.
%  
%                    
%
%         
%-------------------------------------------------------
%
% Inputs:
% xobs: a scalar or an 1D array containing the x-coordinates
%       observation locations
% zobs: a scalar of 1D array of same length and shape as xobs,
%       containing the z-coordinates of the observation locations
%
%       If a single number is give for each variable, 
%       the function returns a sngle value of the gravity 
%       anomaly in mGal
%
%       If a array of numbers (equal lengh) is given to
%       each of these two variable, the function returns 
%       the gravity anomalies at those loctions.
%
%       Note that the z-variable is defined as positive down.
%       So if the surface is z=0 m, than a gravity stations that
%       is 1 m above the surface will have zobs=-1
%
%       Correspondingly, a basement depth of 850 m would have z=850.
%
%
% x:    central x-coordinate of the dyke
% b:    width of the dyke
% zt:   z-coordinate of the top of the dyke
% zb:   z-coordinate of the bottom of the dyke
% rho:  density contrast of the dyke (basin fill)
%
% Output:
% gz:  scalar or 1D array of same length and shape as xobs,
%      containing the calculated gravity value in mGal
% 
% Require a function gploy.m
%
% Coordinate system:
% x-axis: positive easting
% z-axis: positive down vertically
%
%---------------------------------------------------------

% construct the "polygon" representing the vertical dyke
%
xc(1)=x-b/2; zc(1)=zt-0.0001*b;
xc(2)=x+b/2; zc(2)=zt+0.0001*b;
xc(3)=xc(2); zc(3)=zb+0.0001*b;
xc(4)=xc(1); zc(4)=zb-0.0001*b;
nc=4; 

gz=gpoly(xobs, zobs, xc, zc, nc, rho);

% end of function vdyke.m
%

function g=gpoly(xo,zo,xc,zc,nc,rho);
%
% function g=gpoly(xo,zo,xc,zc,nc,rho)
%
% Calculate the gravity field of a 2D prism
% with polygonal cross-section.
%
% input:
%
%    xo, zo: coordinates of the observation point (arrays)
%    xc, zc: coordinates of the vertices of the polygon (arrays)
%            ordered clockwise
%    nc:     number of vertices
%    rho:    density
%
% output:
%
%    g:        calculated gravity data
%
gcons=6.672E-03;

sum=xo.*0.0;

for ic=1:1:nc,

   if ic==nc;
      ic2=1;
   else
      ic2=ic+1;
   end;

   xcb=xc(ic);
   xce=xc(ic2);

   dx=abs(xc(ic2)-xc(ic));
   dz=abs(zc(ic2)-zc(ic));
   if dz > 1.0E-6, 
      zcb=zc(ic);
      zce=zc(ic2);
   else
      zcb=zc(ic)+0.00001*dx;
      zce=zc(ic2);
   end;

   x1=xcb-xo;
   x2=xce-xo;

   z1=zcb-zo;
   z2=zce-zo;

   rt1=x1.*x1 + z1.*z1;
   rt2=x2.*x2 + z2.*z2;

   bot=z2-z1;

   alpha=(x2-x1)./bot;
   beta=(x1.*z2 - x2.*z1)./bot;

   factor=beta./(1.0+alpha.*alpha);

   term1=log(rt2./rt1)./2;
   term2=atan2(z2,x2)-atan2(z1, x1);

   sum=sum + factor.*(term1-alpha.*term2);

end;
     
g=2.0*rho*gcons*sum;

%
% end of function
%