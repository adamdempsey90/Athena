function [x,d,vx,vy]=gathersnaps(dir,np,dim,tlim) 

nt=length(tlim)

for i=1:nt

	dat=loadmpi(dir,np,tlim(i),dim);
	if i==1
		x=dat(:,1);
	end

	d(:,i) = dat(:,2);
	vx(:,i) = dat(:,3);
	vy(:,i) = dat(:,4)+1.5*x;

end	
	

