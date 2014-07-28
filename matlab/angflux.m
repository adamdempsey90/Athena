function [x,Fp,Twb,Th,xp,deltaFp,intTwb,intTh,dxFp,Pibarout,dxPibarout]=angflux(dat)

[x,vx,vy,d,dxvx,dxvy,dxd]=getmeans2(dat,0);
q=1.5; om=1; c=1; 
nu =.003; k=1; xs = 1;
mp = .1;
Nx=length(x);
dx = min(diff(x));

phik = (-mp./pi).*besselk(0,abs(k.*sqrt(x.^2+xs.^2)));


dbar = d(1,:);
vxbar = vx(1,:);
vybar = vy(1,:);

sig = d(k+1,:);
u = vx(k+1,:);
v = vy(k+1,:);

dxdbar = dxd(1,:);
dxvxbar = dxvx(1,:);
dxvybar = dxvy(1,:);

dxsig = dxd(k+1,:);
dxu = dxvx(k+1,:);
dxv = dxvy(k+1,:);


Pibar=zeros(2,2,Nx);
Pip = zeros(2,2,Nx);
Pipp=zeros(2,2,Nx);
dxPibar=zeros(2,2,Nx);
dxPip=zeros(2,2,Nx);
dxPipp=zeros(2,2,Nx);


Pibar(1,1,:) = dbar.*(-c.^2+(4./3).*nu.*dxvxbar);
Pibar(2,1,:) = dbar.*nu.*(dxvybar-q.*om);
Pibar(1,2,:) = Pibar(2,1,:);
Pibar(2,2,:) = dbar.*(-c.^2-(2./3).*nu.*dxvxbar);

Pip(1,1,:) = -c.^2.*sig + nu.*dbar.*((4./3).*dxu-(2./3).*I.*k.*v)+(4./3).*nu.*sig.*dxvxbar;
Pip(2,1,:) = nu.*dbar.*(dxv+I.*k.*u)+nu.*sig.*(dxvybar-q.*om);
Pip(2,2,:) = -c.^2.*sig+nu.*dbar.*((4./3).*I.*k.*v -(2./3).*dxu)-(2./3).*nu.*sig.*dxvxbar;
Pip(1,2,:) = Pip(2,1,:);

Pipp(1,1,:) = nu.*conj(sig).*((4./3).*dxu-(2./3).*I.*k.*v);
Pipp(2,1,:) = nu.*conj(sig).*(dxv+I.*k.*u);
Pipp(2,2,:) = nu.*conj(sig).*((4./3).*I.*k.*v - (2./3).*dxu);
Pipp(1,2,:) = Pipp(2,1,:);

Pipp = real(Pipp);


for i=1:2
for j=1:2
	dxPibar(i,j,:) = calcderiv(dx,Pibar(i,j,:),2);
	dxPip(i,j,:) = calcderiv(dx,Pip(i,j,:),2);
	dxPipp(i,j,:) = calcderiv(dx,Pipp(i,j,:),2);
	dyPibar(i,j,:) = 0;
	dyPip(i,j,:) = I.*k.*Pip(i,j,:);
	dyPipp(i,j,:) = I.*k.*Pipp(i,j,:);
end
end
Pibarout=reshape(Pibar(1,2,:),[1 Nx]);
dxPibarout=reshape(dxPibar(1,2,:),[1 Nx]);

figure; plot(x,reshape(dxPibar(1,2,:),[1 Nx])); ylabel('bar')
figure; plot(x,reshape(dxPip(1,2,:),[1 Nx])); ylabel('dx prime')	
figure; plot(x,reshape(dyPip(1,2,:),[1 Nx])); ylabel('dy prime')	

Fp = 2.*real(vxbar.*conj(sig).*u +dbar.*conj(u).*v);
Twb = 2.*real(dbar.*conj(u).*dxv - conj(sig).*u.*(dxvybar+(2-q)) ...
	- (conj(sig).*sig)./(dbar.^2).*reshape(dxPibar(1,2,:),[1 Nx]));
%	+ (conj(sig)./dbar).*(reshape(dxPip(1,2,:)+dyPip(2,2,:),[1 Nx])));

Th = 2.*real(conj(sig).*I.*k.*phik);

dxFp = calcderiv(dx,Fp,2);
ind = find(x>=0);
xp = x(ind);
intTh = cumtrapz(x(ind),Th(ind));
intTwb = cumtrapz(x(ind),Twb(ind));
deltaFp = Fp(ind)- ones([1,length(ind)]).*Fp(ind(1));

figure; plot(xp,deltaFp,xp,intTh+intTwb); legend('\Delta Fp','\int Th+Twb'); xlabel('x');

figure; plot(x,dxFp,x,Twb+Th); legend('dxFp','Twb+Th'); xlabel('x');


