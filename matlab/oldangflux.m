function [out,Th]=angflux(x,y,vx,vy,d,nu,omega,rs,m)

	Nx = length(find(y==min(y)));
	Ny = length(find(x==min(x)));
	xlist = unique(x);
	ylist = unique(y);
	dx = xlist(2)-xlist(1);
	dy = ylist(2) - ylist(1);
	Ly = max(ylist)-min(ylist);
    Fh = d.*vx.*vy;	
	
	Fx = .5.*d.*omega.*x.*vx;
	
	
% Get partial derivs
	disp('Computing derivatives with forward 4th order method')
%	dyvx=zeros(size(vx));
	dyphi=zeros(size(vx));
	dxvy=zeros(size(vy));	
	for j = 1:Nx
%		disp(['Working on x=',num2str(xlist(j))])
		ind = find(x == xlist(j));
%		dyvx(ind) = deriv1D(dy,vx(ind),'forward',4,'periodic');
%		dyphi(ind) = deriv1D(dy,phi(ind),'forward',4,'periodic');
		dyphi(ind) = potential(xlist(j),y(ind),rs,m,'y');

	end
	for j=1:Ny
%		disp(['Working on y=',num2str(ylist(j))])
		ind = find(y == ylist(j));
		dxvy(ind) = deriv1D(dx,vy(ind),'forward',4,'zero');
	end


	Fnu = -nu.*d.*dxvy;
%	Fnuy = -nu.*d.*dyvx;
	dxdTphi = -d.*dyphi;
	
	disp('Computing integrals')
	Fhavg=yavg(x,y,Fh).*Ly;
	Fxavg=yavg(x,y,Fx).*Ly;
%	Fnuomegaavg=yavg(x,y,Fnuomega).*Ly;
	Fnuavg = yavg(x,y,Fnu).*Ly;
%	Fnuyavg = yavg(x,y,Fnuy).*Ly;
	dxdTphiavg = yavg(x,y,dxdTphi).*Ly;

	inds = find(xlist>=0);
	Th = cumtrapz(xlist(inds),dxdTphiavg(inds));

	out = [Fhavg Fxavg Fnuavg dxdTphiavg, ...
		  Fhavg+Fnuavg, Fhavg+Fnuavg+Fxavg];

%	LHSTorque = deriv1D(dx,out(:,6),'forward',4,'zero');

% Plot
	labels={'< \Sigma u_x u_y > q^{-2}', ...
			'< 0.5 \Sigma \Omega x u_x > q^{-2}', ...
			'< - \nu \Sigma \partial_x u_y > q^{-2}', ...
			'< - \Sigma \partial_y \phi >', ...
			'Ang. Momentum Flux + Viscous flux', ...
			'Tot Ang Momentum Flux' };
			
	for i=1:length(labels)
		figure(i)
		plot(xlist,out(:,i)./m.^2,'-.'); 
		xlim([0,max(xlist)]);
		xlabel('x'); ylabel(labels{i}); title(labels{i});
	end
		figure(length(labels)+1)
		plot(xlist(inds),Th./m.^2,'-.'); 
		xlabel('x');
		ylabel('< \int_0^x dT/dx dx > q^{-2}'); title('< \int_0^x dT/dx dx > q^{-2}');
	
	
	figure(length(labels)+2)
	plot(xlist(inds),Th./m.^2,xlist(inds),out(inds,1)./m.^2);
	legend('-Th','Fh'); xlabel('x');
	
%	figure(length(labels)+3)
%	plot(xlist(inds),LHSTorque(inds)./m.^2); xlabel('x'); title('\partial_x < Tot LHS Flux >')
	
end

