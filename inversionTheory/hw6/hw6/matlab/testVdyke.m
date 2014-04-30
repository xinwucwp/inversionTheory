load dat.txt
nd = 48;
nm = 120;
xc = -5950:100:6000;
xo = dat(:,1);
zo = -2+zeros(nd,1);
zt = zeros(nm,1);
mk = 100+zeros(nm,1);
mk(120)=1.0;
dk = zeros(nd,1);
wd = 100.0;
dc = 0.1;
for i=1:120
    xci = xc(i);
    zti = zt(i);
    zbi = mk(i);
    dki = vdyke(xo,zo,xci,wd,zti,zbi,dc);
    dk = dk+dki;
end
