function [xlist, vxhat,vyhat,dhat,dxvx,dxvy,dxd]=getmeans2(dat,plotflag)

x=dat(:,1); y=dat(:,2); d=dat(:,3); vx=dat(:,4); 
vy=dat(:,5)+1.5.*x;
%vy =dat(:,5);
xlist=unique(x'); ylist=unique(y');
Nx=length(xlist); Ny=length(ylist);
Lx = max(xlist)-min(xlist);
for i=1:Nx
  ind = find(x == xlist(i));
  vxhat(:,i) = fft(vx(ind))./Ny;
  vyhat(:,i) = fft(vy(ind))./Ny;
  dhat(:,i) = fft(d(ind))./Ny;

end

dx = mean(diff(xlist));
for i=1:Ny
	dxvx(i,:) = centderiv(dx,vxhat(i,:),10);
	dxvy(i,:) = centderiv(dx,vyhat(i,:),10);
	dxd(i,:) = centderiv(dx,dhat(i,:),10);
end



if plotflag==1
h1=figure;
subplot(2,2,1); plot(xlist,vxhat(1,:)); title('y-averged vx');
subplot(2,2,2); plot(xlist,vyhat(1,:)); title('y-averaged vy');
subplot(2,2,3);plot(xlist,dhat(1,:)); title('y averaged d');

h2=figure;
subplot(2,2,1); plot(xlist,-real(vxhat(2,:)),'-b',xlist,-imag(vxhat(2,:)),'-g');
title('k=1 component for vx');xlim([-Lx/8,Lx/8]);
subplot(2,2,2); plot(xlist,-real(vyhat(2,:)),'-b',xlist,-imag(vyhat(2,:)),'-g'); 
title('k=1 component for vy');xlim([-Lx/8,Lx/8]);
subplot(2,2,3);plot(xlist,-real(dhat(2,:)),'-b',xlist,-imag(dhat(2,:)),'-g');
title('k=1 component for d');xlim([-Lx/8,Lx/8]);

end
