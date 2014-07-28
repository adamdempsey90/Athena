function [vhat,KK,XX]=transform(x,y,v,Ly);
% 2D data of velocity
% Want vhat(ky,x). Do FFT in y direction for a given x
% End with an array of ky's and x and FFT(v);
% First index of vhat specifies ky
% Second index of vhat specifies x
% vhat(1,:) -> vhat(ky(1)) as function of x
% vhat(:,1) -> vhat(x(1)) as function of ky

Nx=length(find( y == min(y)));;
Ny=length(find( x == min(x))); 

xlist = unique(x);
xvals = xlist(find(xlist >= 0));


ky = (0:(Ny/2 - 1)).* (2*pi/Ly);
%xvals = x(Nx/2:Nx); % Take positive half of x, should be symmetric


for j=1:length(xvals)

	ind = find(x == xvals(j));
	vfft = fft(v(ind));
	vhat(:,j) = vfft(1:Ny/2);
end		
	
[XX,KK]=meshgrid(xvals,ky);



