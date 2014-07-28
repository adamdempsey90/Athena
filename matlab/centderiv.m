function out=centderiv(dx,in,order)

% centered diff with zero gradient b.c %
out = zeros(size(in));
ind = 1:length(in);

s = size(in);
ng = order/2;

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



bcl = in(ng:-1:1);
bcr = in(end:-1:end-ng+1);

if s(1)>s(2)
	temp=[bcl;in;bcr];
else
	temp=[bcl in bcr];
end
temp2=zeros(size(temp));

for j=1:2*ng
	temp2(ind+ng) += temp(ind+ng+pos(j)).*coefs(j);
end

out = temp2(ng+1:end-ng)./dx;


