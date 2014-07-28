
t='0067';
massvals={'2.0' '1.5' '1.0' '0.5' '.1' '.05'};
nuvals={'0','.0001' '.001' '.01' '.1'};

for i=1:length(massvals)
	mass(i) = str2num(massvals{i});
end
for j=1:length(nuvals)
	nu(j) = str2num(nuvals{j});
end


dir = cell(length(massvals),length(nuvals));
FWHM = zeros(length(massvals),length(nuvals));
gapdepth = zeros(length(massvals),length(nuvals));

for i=1:length(massvals)
	for j=1:length(nuvals)
		disp('Working on');
		disp(['mass=',massvals{i},' nu=',nuvals{j}])
		dir{i,j} = ['m',massvals{i},'nu',nuvals{j}];
		
		temp=importdata(['../bin/',dir{i,j},'/PlanetDisk.',t,'.tab'],' ', 6);
		dat=temp.data;
%		temp=importdata(['../bin/',dir{i,j},'/PlanetDisk.',t,'.flux.tab'],' ', 1);
%		datavg=temp.data;
		
 		x=dat(:,3); y=dat(:,4); d=dat(:,5); vx=dat(:,6); vy=dat(:,7); vy = vy + 1.5.*x;
%		vxavg=datavg(:,6); vyavg=datavg(:,7); densavg=datavg(:,8);
		xlist=unique(x); ylist=unique(y);
		[xx,yy]=meshgrid(xlist,ylist); dd=gridq(x,y,d); vvx=gridq(x,y,vx); vvy=gridq(x,y,vy);

		for n=1:length(xlist)
			ind = find( x == xlist(n));
			densavg(:,n) = trapz(y(ind),d(ind))./40.0;
		end
		posx = find(xlist >= 0);
		negx = find( xlist < 0);
		
		HM = (max(densavg) + min(densavg))/2;
		indl = find( densavg(negx) <= HM,1);
		indr = posx(1)+find( densavg(posx) >= HM,1);

		FWHM(i,j) = xlist(indr) - xlist(indl);
		gapdepth(i,j) = HM;





	end
end

