dtheta = pi/10;
ntheta = 20;
r1 = zeros(21,1);
r2 = zeros(21,1);
r3 = zeros(21,1);
r4 = zeros(21,1);
theta = zeros(21,1);
for i=0:ntheta
    thetai = i*dtheta;
    theta(i+1) = thetai;
    xi = cos(thetai);
    yi = sin(thetai);
    vi = [xi yi];
    r1(i+1) = norm(vi,1.0);
    r2(i+1) = norm(vi,1.5);
    r3(i+1) = norm(vi,2.0);
    r4(i+1) = norm(vi,3.0);
end
h=figure
subplot(1,1,1)
h1=polar(theta,r1,'r-');
set(h1,'linewidth',3)
hold on
h2=polar(theta,r2,'g-');
set(h2,'linewidth',3)
hold on
h3=polar(theta,r3,'b-');
set(h3,'linewidth',3)
hold on
h4=polar(theta,r4,'c-');
set(h4,'linewidth',3)