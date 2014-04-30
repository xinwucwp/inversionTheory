%forward model
function gz = grav(zb)
load dat.txt;
%zb=100*ones(120,1);
xobs = dat(:,1);
%xobs=xobs1';
%zobs = dat(:,2);

gra= dat(:,3);

%xobs=-6000:250:6000;
[k,~]=size(xobs);
[m,~]=size(zb);
x = linspace(-6000.0,6000.0,m);
zobs=-2.0*ones(k,1);
b=100;
zt=0;
%zt=0*ones(k,1);
m_st=200*ones(k,1);
gz = zeros(k,1);
rho=-1.0;
for i= 1:m
    gz=gz+vdyke(xobs,zobs,x(i),b,zt,zb(i),rho);
end

%function gz=vdyke(xobs, zobs, x, b, zt, zb, rho)