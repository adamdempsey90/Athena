function out=gridq(x,y,in)

xlist=unique(x);
ylist=unique(y);

out = zeros(length(ylist),length(xlist));

for i=1:length(xlist)
	ind=find(x == xlist(i));
	out(:,i) = in(ind);
end

