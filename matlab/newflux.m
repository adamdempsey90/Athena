% dat file already loaded %

c=1; om=1.0002; q=1.5; xs=.6*.05;
nu=.003;
mp=4;
k=1;

[x,vx1,vy1,d1,dxvx1,dxvy1,dxd1]=getmeans2(dat,0);
x=transpose(x);
vx = transpose(vx1(1,:));
vy = transpose(vy1(1,:));
d = transpose(d1(1,:));
u = transpose(vx1(k+1,:));
v = transpose(vy1(k+1,:));
sig = transpose(d1(k+1,:));

dxvx = transpose(dxvx1(1,:));
dxvy = transpose(dxvy1(1,:));
dxd = transpose(dxd1(1,:));
dxu = transpose(dxvx1(k+1,:));
dxv = transpose(dxvy1(k+1,:));
dxsig = transpose(dxd1(k+1,:));


% Smooth vx and dxvx

dlmwrite('vx.dat',vx,'delimiter','\t');
dlmwrite('dxvx.dat',dxvx,'delimiter','\t');

system('module load matlab');
system('matlab -nodisplay -nosplash -nodesktop -r "try, smoothing;end, exit"');

vx=load('smoothedvx.dat');
dxvx=load('smootheddxvx.dat');


dx = min(diff(x));
Nx = length(x);


Pibar=zeros(2,2,Nx);
Pip=zeros(2,2,Nx);
Pipp=zeros(2,2,Nx);
Pippb=zeros(2,2,Nx);

dxPibar=zeros(2,2,Nx);
dxPip=zeros(2,2,Nx);
dxPipp=zeros(2,2,Nx);
dxPippb=zeros(2,2,Nx);


Pibar(1,1,:) = -c*c.*d+(4./3).*nu.*d.*dxvx;
Pibar(1,2,:) = nu.*d.*(dxvy-q*om);
Pibar(2,2,:) = -c*c*d - (2./3).*nu.*d.*dxvx;
Pibar(2,1,:) = Pibar(1,2,:);

Pip(1,1,:) = -c*c*sig + nu*d.*(2./3).*(2.*dxu-I.*k.*v)+(4./3).*nu.*sig.*dxvx;
Pip(1,2,:) = nu.*d.*(dxv+I.*k.*u)+nu.*sig.*(dxvy-q.*om);
Pip(2,2,:) = -c*c*sig + nu.*d.*(2./3).*(2.*I.*k.*v-dxu)-(2./3).*nu.*sig.*dxvx;
Pip(2,1,:) = Pip(1,2,:);


Pipp(1,1,:) = (4./3).*dxu - (2./3).*I.*k.*v;
Pipp(1,2,:) = dxv + I.*k*u;
Pipp(2,2,:) = (4./3).*I.*k*v - (2./3).*dxu;
Pipp(2,1,:) = Pip(1,2,:);

for i=1:2
for j=1:2
	
	Pipp(i,j,:)  = nu.*sig.*reshape(Pipp(i,j,:),[Nx 1]);
	Pippb(i,j,:) = nu*conj(sig).*reshape(Pipp(i,j,:),[Nx 1]);
	dxPibar(i,j,:) = centderiv(dx,reshape(Pibar(i,j,:),[Nx 1]),10);
	dxPip(i,j,:) = centderiv(dx,reshape(Pip(i,j,:),[Nx 1]),10);
	dxPipp(i,j,:) = centderiv(dx,reshape(Pipp(i,j,:),[ Nx 1]),10);
	dxPippb(i,j,:) = centderiv(dx,reshape(Pippb(i,j,:),[Nx 1]),10);
end
end

phik=@(x) (-2*mp./pi).*besselk(0,abs(k.*sqrt(x.^2+xs.^2)));

divP = dxPip(1,2,:)+I.*k.*Pip(2,2,:);


twb1 = d.*conj(u).*dxv;
twb2 = -conj(sig).*u.*(dxvy+(2-q).*om);
twb3 = -(conj(sig).*sig./(d.*d)).*reshape(dxPibar(1,2,:),[Nx 1]);
twb4 = (conj(sig)./d).*reshape(divP,[Nx 1]);

Twb = 2*real(twb1+twb2+twb3+twb4);

Fp = 2.*real(vx.*conj(sig).*v+d.*conj(u).*v);
Th = 2.*real(conj(sig).*I.*k.*phik(x));

ind = find(x >=0);
xp = x(ind);
iTh = cumtrapz(xp,Th(ind));
iTwb = cumtrapz(xp,Twb(ind));

dxFp = centderiv(dx,Fp,10);

Fpn = Fp(ind) - Fp(ind(1));

figure; plot(x,dxFp+Th,x,Twb); xlim([-10 10]); ylim([-.0001 .0001])
xlabel('x'); legend('dxFp-Th','Twb');
figure; plot(xp,Fpn+iTh,xp,iTwb); xlabel('x'); 
legend('\Delta F - \int Th','\int Twb');
