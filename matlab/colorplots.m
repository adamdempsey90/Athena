%n='0001';
%dir='../bin/q2_nu001/';
%dat=load([dir,'PlanetDisk.',n,'.tab']);
%x=dat(:,3); y=dat(:,4); d=dat(:,5); vx=dat(:,6); vy=dat(:,7); vy=vy+1.5.*x; 
x=dat(:,1); y=dat(:,2); d=dat(:,3); vx=dat(:,4); vy=dat(:,5)+1.5.*dat(:,1);
xlist=unique(x); ylist=unique(y);
[xx,yy]=meshgrid(xlist,ylist); dd=gridq(x,y,d); vvx=gridq(x,y,vx); vvy=gridq(x,y,vy);

figure; pcolor(xx,yy,log10(dd)); shading interp; colorbar
hold on; quiver(xx,yy,vvx,vvy)
xlabel('x'); ylabel('y'); title('\Sigma');
hold off
figure; pcolor(xx,yy,vvx); shading interp; colorbar
xlabel('x'); ylabel('y'); title('v_x');

figure; pcolor(xx,yy,vvy); shading interp; colorbar
xlabel('x'); ylabel('y'); title('v_y');


quiver(xx,yy,vvx,vvy)

