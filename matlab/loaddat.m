m='0040';
n='0040';


dat=load(['../bin/PlanetDisk.',n,'.tab']);
datavg=load(['../bin/PlanetDisk.',m,'.flux.tab']);

x=dat(:,3); y=dat(:,4); d=dat(:,5); vx=dat(:,6); vy=dat(:,7); vy = vy + 1.5.*x;
vxavg=datavg(:,6); vyavg=datavg(:,7); densavg=datavg(:,8);
xlist=unique(x); ylist=unique(y);
[xx,yy]=meshgrid(xlist,ylist); dd=gridq(x,y,d); vvx=gridq(x,y,vx); vvy=gridq(x,y,vy);

ddavg=zeros(size(dd));
for i=1:length(ddavg(:,1))
	ddavg(i,:)=densavg(:);
end

nu=0.1;
omega=1;
q=.05;
rs=3.0;
c=1.0;

Ly = 40;
kyi = 1.0;




Lx = max(xlist)-min(xlist);

%Ly = max(ylist)-min(ylist);

Nx = length(xlist);
Ny = length(ylist);


% Take out shear
vy = vy + 1.5.*x;

%[flux,Th]=angflux(x,y,vx,vy,d,nu,omega,rs,q);

davg=datavg(:,8);

%vxavg=datavg(:,6);
%vyavg=datavg(:,7);

%[vyhat,kky,xxy]=transform(x,y,vy,Ly);
%[vxhat,kkx,xxx]=transform(x,y,vx,Ly);
%[dhat,kkd,xxd]=transform(x,y,d-ones(size(d)),Ly);

vyhat=zeros([Ny Nx]);
vxhat=zeros([Ny Nx]);
dhat=zeros([Ny Nx]);

for i=1:Nx
	ind=find(x == xlist(i));
	vyhat(:,i) = fft(vy(ind))/Ny;
	vxhat(:,i) = fft(vx(ind))/Ny;
	dhat(:,i) = fft(d(ind)-1)/Ny;
end

%Tk = 4.*kkd.*besselk(0,abs(xxd.*kkd)).*q .* imag(dhat);

%ind=find(kky(:,1)>=kyi*2*pi/Ly,1);



%normfactor = (1/(Ny));

%normfactor = 1;


x1=x(find(x >= 1,1));
ind2 = find(x == x1);

figure; plot(xlist,-real(vyhat(2,:)),xlist,-imag(vyhat(2,:))); xlabel('x'); ylabel('(v/c)'); legend('Re(v)', 'Im(v)');
xlim([0,max(xlist)]);
%ylim([-1,1]);

figure; plot(xlist,-real(vxhat(2,:)),xlist,-imag(vxhat(2,:))); xlabel('x'); ylabel('(u/c)'); legend('Re(u)', 'Im(u)');
xlim([0,max(xlist)]);
%ylim([-1,1]);

figure; plot(xlist,-real(dhat(2,:)),xlist,-imag(dhat(2,:))); xlabel('x'); ylabel('\delta \Sigma'); 
legend('Re(\delta \Sigma)', 'Im(\delta \Sigma)'); xlim([0,max(xlist)]);


figure; pcolor(xx,yy,dd); shading interp; colorbar
xlabel('x'); ylabel('y'); title('\Sigma');

figure; pcolor(xx,yy,vvx); shading interp; colorbar
xlabel('x'); ylabel('y'); title('v_x');

figure; pcolor(xx,yy,vvy); shading interp; colorbar
xlabel('x'); ylabel('y'); title('v_y');

figure; pcolor(xx,yy,dd-ddavg); shading interp; colorbar
xlabel('x'); ylabel('y'); title('\Sigma - < \Sigma >');

%ylim([-1,1]);

%figure; plot(xxd(ind,:),Tk(ind,:)./q.^2); xlabel('x'); ylabel('( \frac{dT}{dx} )_{k}'); 
%legend('Tk');

%figure; plot(y(ind2) - .75.*x1.^2,(d(ind2)-1)./q,'-o'); xlabel('y - 0.75 x^2'); ylabel(' d1/d ) * Mth/Mp at x=1');
%figure; plot(xxx(ind,:),real(vxhat(ind,:)).*normfactor,'-x'); xlabel('x'); ylabel('(u/c) M_{th}/M_p'); legend('Re(u)');
%ylim([-1,1]);
