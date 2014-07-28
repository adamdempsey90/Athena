function [xlist, vxhat,vyhat,dhat]=getmeans(fname)

dat=load(fname);
x=dat(:,3); y=dat(:,4); d=dat(:,5); vx=dat(:,6); vy=dat(:,7)+1.5.*x;
xlist=unique(x); ylist=unique(y);
Nx=length(xlist); Ny=length(ylist);
for i=1:Nx
  ind = find(x == xlist(i));
  vxhat(:,i) = fft(vx(ind))./Ny;
  vyhat(:,i) = fft(vy(ind))./Ny;
  dhat(:,i) = fft(d(ind))./Ny;
end

