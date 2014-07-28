function out=centderiv(dx,in,order)

% centered diff with zero gradient b.c %
out = zeros(size(in));
ind = 1:length(in);

s = size(in);
ng = order/2;
pos = -ng:ng;
if order==2
	coefs=[-.5 0 .5];
elseif order==4
	coefs=[1./12 -2./3 0 2./3 -1./12];
elseif order==6
	coefs=[-1./60 3./20 -1./4 0 1./4 -3./20 1./60];

elseif order==8
	coefs=[1./280 -4./105 1./5 -4./5 0 4./5 -1./5 4./105 -1./280];

elseif order==10
	coefs=[-2. 25. -150. 600. -2100. 0 2100. -600. 150. -25. 2.];
	coefs=coefs./2520;
end

D = spdiags(coefs.*ones([max(s) 1])./dx,pos,max(s),max(s));
D = [zeros([ng max(s)]); D; zeros([ng max(s)])];



bcl = in(ng:-1:1,:);
bcr = in(end:-1:end-ng+1,:);

in = [bcl; in; bcr];

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


