function [out1,out2]=submean(dat) 

x=dat(:,1); y=dat(:,2); d=dat(:,3); vx=dat(:,4); vy=dat(:,5);
dat(:,5)=dat(:,5)+1.5.*dat(:,1);
xlist=unique(x); ylist=unique(y);
Nx=length(xlist); Ny=length(ylist); Ly=max(ylist)-min(ylist);

datbar=zeros(size(xlist)+[0 3]);
datprime = zeros(size(dat(:,3:5)));
disp(size(datbar));
for i=1:Nx
	ind=find(x == xlist(i));
	for j=3:5
		datbar(i,j-2) = trapz(ylist,dat(ind,j))./Ly;
		datprime(ind,j-2) = dat(ind,j) - datbar(i,j-2); 	
	end
end


dx = min(diff(xlist));

dxvx = centderiv(dx,vx,10);
dxvy = centderiv(dx,vy,10);
dxd = centderiv(dx,vy,10);

ux = datprime(:,1); uy=datprime(:,2); sig=datprime(:,3);

dxux = centderiv(dx,ux,10);
dxvx = centderiv(dx,vx,10);
dxsig = centderiv(dx,sig,10);

input = [x d
[dxPibarxy, divP] = calcPi(x,input);

dxTwb = -dbar.*ux.*dxvx+sig.*ux.*(dxvy+.5)+(sig./d).^2 .*dxPibarxy-(sig./d).*divP; 


