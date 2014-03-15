ma = [0.5,1.0,1.5,1.0,0.5,0.0,-0.5,-1.0,-0.5,0.0,10.0];
phi = zeros(4,1);
pv = zeros(4,1);
pv(1) = 1.0;
pv(2) = 1.5;
pv(3) = 2.0;
pv(4) = 3.0;
for i=1:4
    phi(i) = norm(ma,pv(i))^pv(i);
end
%
%fig1 = figure;
%set(fig1,'position',[900,500,600,300]);
%plot(p,phi,'r+','linewidth',3)
%hold on;
%plot(p,phi,'r-','linewidth',2.0)
%ylabel('phi(ma)','fontsize',32)
%xlabel('p=1.0,1.5,2.0,3.0','fontsize',32);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gma = zeros(4,11);
ma = ma + 0.01;
for i=1:4
    p = pv(i);
    for j=1:11
      maj = ma(j);    
      gma(i,j)=p*maj*(abs(maj))^(p-2.0);
      display(p)
    end
end
close all;
%fig2 = figure;
%set(fig2,'position',[900,500,600,300]);
figure('position',[900,500,600,200])
plot(1:11,gma(1,:),'r+','linewidth',3);
hold on;
plot(1:11,gma(1,:),'r-','linewidth',2);
xlabel('index')
ylabel('grad(phi)')
figure('position',[900,500,600,200])
plot(1:11,gma(2,:),'g+','linewidth',3);
hold on;
plot(1:11,gma(2,:),'g-','linewidth',2);
xlabel('index')
ylabel('grad(phi)')
figure('position',[900,500,600,200])
plot(1:11,gma(3,:),'b+','linewidth',3);
hold on;
plot(1:11,gma(3,:),'b-','linewidth',2);
xlabel('index')
ylabel('grad(phi)')
figure('position',[900,500,600,200])
plot(1:11,gma(4,:),'c+','linewidth',3);
hold on;
plot(1:11,gma(4,:),'c-','linewidth',2);
xlabel('index')
ylabel('grad(phi)')