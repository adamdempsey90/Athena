function out=reorder(in)
% Reorder fields if they were done in parallel
% Sort in x given y
	x=in(:,1); y=in(:,2);
	Ny=length(find(y==y(1)));
	Nx=length(find(x==x(1)));
	out = zeros(size(in));

	ylist=unique(y);	% Sorted list of unique y values
	for j=1:Ny
		ind=find(y==ylist(j));
		unsortedx=in(ind,:);
		[d1,d2]=sort(unsortedx(:,1));		% Sort by x
		out((j-1)*Nx+1:j*Nx,:)=unsortedx(d2,:);	
	end
end
