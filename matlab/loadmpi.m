function out=loadmpi(dir,np,time,dim) 

if dim==1
	indstart=2;
elseif dim==2
	indstart=3;
else
	indstart=4;
end


if time<10
	tstr=['000',num2str(time)];
elseif time<100
	tstr=['00',num2str(time)];
elseif time<1000
	tstr=['0',num2str(time)];
else
	tstr=num2str(time);
end

mat=cell(np,1);
for id=0:np-1
	if id~=0
		fname=[dir,'id',num2str(id),'/PlanetDisk-id',num2str(id),'.',tstr,'.tab'];
	else
		fname=[dir,'id0','/PlanetDisk.',tstr,'.tab'];
	end
	dum=load(fname);

	mat{id+1} = dum(:,indstart:end);
	if id==0
		out=mat{id+1};
	else
		out=[out;mat{id+1}];
	end
end


	

