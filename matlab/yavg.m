function out=yavg(x,y,q)

% Average along strips at a fixed x
	Nx = length(find(y==min(y)));
	xvals=unique(x);
	dx=xvals(2)-xvals(1);
	out=zeros(Nx,1);
	for j=1:Nx
		ind=find(x==xvals(j));
		dy=max(y(ind))-min(y(ind));
		out(j)= trapz(y(ind),q(ind))./dy;
	end



end	


