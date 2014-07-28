function [Twb,Th,Fp,Fpn,dxFp,iTh,iTwb,x,y,vxbar,vybar,dbar,u,v,sig]=calctwb(dat,plotflag)


q=1.5; nu=.003; mp=4; om=sqrt(1+4./(20.^3)); c=1; k=20; xs=.6*.05;

x=dat(:,1); y=dat(:,2); d=dat(:,3); vx=dat(:,4); vy=dat(:,5)+1.5.*x;
xlist=unique(x); ylist=unique(y);
Nx=length(xlist); Ny=length(ylist);
dx = min(diff(xlist)); dy=min(diff(ylist));
xnew = zeros([Nx Ny]);
ynew = zeros([Nx Ny]);
dnew = zeros([Nx Ny]);
vxnew = zeros([Nx Ny]);
vynew = zeros([Nx Ny]);

for i=1:Ny
	ind=find(y == ylist(i));
	xnew(:,i) = x(ind);
	ynew(:,i) = y(ind);
	vxnew(:,i) = vx(ind);
	vynew(:,i) = vy(ind);
	dnew(:,i) = d(ind);

end

x=xnew; y=ynew; vx=vxnew; vy=vynew; d=dnew;

vxbar = mean(vx,2);
vybar = mean(vy,2);
dbar = mean(d,2);
vybarf = vybar.*ones([1 Ny]);
dbarf = dbar.*ones([1 Ny]);

dxdbar = D(dx,dbar,10,1,'out');
dxvxbar = D(dx,vxbar,10,1,'out');
dxvybar = D(dx,vybar,10,1,'out');
dxvxbarf = dxvxbar.*ones([1 Ny]);
dxvybarf = dxvybar.*ones([1 Ny]);
dxdbarf = dxdbar.*ones([1 Ny]);


% Smooth vx and dxvx
%dlmwrite('vx.dat',vxbar,'delimiter','\t');
%dlmwrite('dxvx.dat',dxvxbar,'delimiter','\t');

%system('module load matlab');
%system('matlab -nodisplay -nosplash -nodesktop -r "try, smoothing;end, exit"');

%vxbar=load('smoothedvx.dat');
%dxvxbar=load('smootheddxvx.dat');
vxbarf = vxbar.*ones([1 Ny]);
dxvxbarf = dxvxbar.*ones([1 Ny]);

u = vx - vxbarf;
v = vy - vybarf;
sig = d - dbarf;

disp(mean(mean(u,2)))
disp(mean(mean(v,2)))
disp(mean(mean(sig,2)))
dxu = D(dx,u,10,1,'out');
dyu = D(dy,u,2,2,'periodic');
dxv = D(dx,v,10,1,'out');
dyv = D(dy,v,2,2,'periodic');
dxsig = D(dx,sig,10,1,'out');
dysig = D(dy,sig,2,2,'periodic');

Pipxy = nu.*dbarf.*(dxv+dyu)+nu.*sig.*(dxvybarf-1.5);
Pipyy = -c*c*sig+nu.*dbarf.*((4./3).*dyv-(2./3).*dxu)-(2./3).*nu.*sig.*dxvxbarf;

Pibarxy = dbarf.*nu.*(dxvybarf-1.5);

dxPibarxy = D(dx,Pibarxy,10,1,'out');
divP = D(dx,Pipxy,10,1,'out')+D(dy,Pipyy,2,2,'periodic');



phi = (-2*mp./pi).*besselk(0,abs(k.*sqrt(x.^2+xs.^2))).*cos(k.*y);
phip = phi - mean(phi,2).*ones([1 Ny]);

dyphi = (2*k.*mp./pi).*besselk(0,abs(k.*sqrt(x.^2+xs.^2))).*sin(k.*y);
dyphip = D(dy,phip,2,2,'periodic');

Th1 = sig.*dyphi;
Th = mean(Th1,2);
Thp=mean(sig.*dyphip,2);

twb1=dbarf.*u.*dxv - sig.*u.*(dxvybarf+.5)-(sig./dbarf).^2.*dxPibarxy+(sig./dbarf).*divP;

Fp1 = vxbarf.*sig.*v+dbarf.*u.*v;
Fp=mean(Fp1,2); 
dxFp=D(dx,Fp,10,1,'out');
Twb=mean(twb1,2);
plot(x,Fp,x,dxFp); legend('Fp','dxFp'); xlim([-10 10]);
ind = find(xlist>=0);
xp=xlist(ind);
iTwb = cumtrapz(xp,Twb(ind));
iTh  = cumtrapz(xp,Th(ind));
Fpn = Fp(ind)-Fp(ind(1));



if (plotflag)
	figure; plot(x,dxFp+Th,x,Twb); xlabel('x'); legend('dxFp-Th','Twb');
	xlim([-10 10]);

	figure; plot(xp,Fpn+iTh,xp,iTwb); xlabel('x'); legend('\Delta Fp - \int Th','\int Twb');
	figure; plot(x,dxFp+Thp,x,Twb); xlabel('x'); legend('dxFp-Thp','Twb');
	xlim([-10 10]);
end


% Test D with sin(x)
%
%q=linspace(0,2*pi,100); q=q';
%a=D(min(diff(q)),sin(q),10,1,'out');
%figure; plot(q,a,'-x',q,cos(q))

end

function out=D(dx,in,order,dim,bc)

if mod(order,2)~=0
	disp('Must have an even order');
	out=0;
	return;
end
ng = order/2;
out=zeros(size(in));

s=size(in);
ind=1:s(dim);

if order==2
	coefs=[-.5 .5];
	pos = [-1 1];
elseif order==4
	coefs=[1./12 -2./3 2./3 -1./12];
	pos = [-2 -1 1 2];
elseif order==6
	coefs=[-1./60 3./20 -1./4 1./4 -3./20 1./60];
	pos = [-3 -2 -1 1 2 3];
elseif order==8
	coefs=[1./280 -4./105 1./5 -4./5 4./5 -1./5 4./105 -1./280];
	pos = [-4 -3 -2 -1 1 2 3 4];
elseif order==10
	coefs=[-2. 25. -150. 600. -2100. 2100. -600. 150. -25. 2.];
	coefs=coefs./2520;
	pos=[-5 -4 -3 -2 -1 1 2 3 4 5];
end



if dim==1
if strcmp(bc,'out')
	bcl = in(ng:-1:1,:);
	bcr = in(end:-1:end-ng+1,:);
	temp=[bcl;in;bcr];
elseif strcmp(bc,'periodic')
	bcl = in(end-ng+1:end,:);
	bcr = in(1:ng,:);
	temp=[bcl;in;bcr];
else
	disp('BC not supported');
	return;
end
elseif dim==2
if strcmp(bc,'out')
	bcl = in(:,ng:-1:1);
	bcr = in(:,end:-1:end-ng+1);
	temp=[bcl in bcr];
elseif strcmp(bc,'periodic')
	bcl = in(:,end-ng+1:end);
	bcr = in(:,1:ng);
	temp=[bcl in bcr];
else
	disp('BC not supported');
	return;
end
end
temp2 = zeros(size(temp));
for i=1:s(2-mod(dim+1,2))
	for j=1:2*ng
		if dim==1
			temp2(ind+ng,i)+=temp(ind+ng+pos(j),i).*coefs(j);
		elseif dim==2
			temp2(i,ind+ng)+=temp(i,ind+ng+pos(j)).*coefs(j);
		end
	end
	if dim==1
		out = temp2(ng+1:end-ng,:)./dx;
	elseif dim==2
		out=temp2(:,ng+1:end-ng)./dx;
	end

end

end
