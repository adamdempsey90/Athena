function [uvel,vvel,sig,x,time]=pdesolver(N,Lx,k,dt,endt)
% 
%
% (I + A dt ) X^(n+1) = X^n - dt B
%
%
	global om c nu q k Mp xsoft Lx N
	 
	om=1; c=1; nu=.1; q-1.5; xsoft=1.0;

	h = Lx/(N-1);
	
	x = -Lx/2 + (0:N-1)'.*h;
	
% Make A matrix
%	Make Derivative Matrices using 4th order finite differences

	Dc = [1/12 -2/3 0 2/3 -1/12];
	D2c = [-1/12 4/3 -5/2 4/3 -1/12];
	
	ng = (length(Dc)-1);
	
	
	A = zeros(3*(N+4),3*(N+4));
	D = zeros(N+ng,N+ng);		% 2 ghost zones on either boundary. 
	D2 = zeros(N+ng,N+ng);
	ident = eye(N+ng);
	zmat = zeros(N+ng,N+ng);
	
	for j=1:length(Dc)
% Room to change the stencil for higher/lower order methods
% Note # of ghost zones has to change if order increased.
		D += diag(Dc(j).*ones(N+ng,1),j-3);
		D2 += diag(D2c(j).*ones(N+ng,1),j-3);
	end
	

	ident=clearghost(ident,ng); 
	D=clearghost(D,ng);
	D2=clearghost(D2,ng);
	pmat=zmat; dpmat=zmat;
	pmat(1+ng/2:N+ng/2,1+ng/2:N+ng/2) = diag(phi(x,0));
	dpmat(1+ng/2:N+ng/2,1+ng/2:N+ng/2) = diag(phi(x,1));
	pmat=clearghost(pmat,ng);
	dpmat=clearghost(dpmat,ng);

	A = [ -nu.*(D2 - k.^2.*ident) , -2.*om.*ident, D ; ...
		  2.*(1-.5*q)*om.*ident, -nu.*(D2 - k.^2.*ident), 1.0i.*k.*ident; ...
		  -c.^2.*D, -1.0i .*k.*c.^2.*ident, zmat];
	
		  
	B = [ pdmat, zmat, zmat; ...
		  zmat, 1.0i.*k.*pmat, zmat; ...
		  zmat, zmat, zmat ];
		  
		  
% Set up time stepping and solution matrix

	nt = endt/dt + 1;
	time = (0:nt-1)'.*dt;
	sol = zeros(3*(N+ng),nt);
	ident2 = clearghost(eye(3*(N+ng)),ng); 
	for t=2:nt
		disp(['Step ',num2str(t-1),'of ',num2str(nt),', time=',num2str(time(t)),', dt=',num2str(dt)];	
		sol(:,t) = (ident2 + A.*dt) \ (sol(:,t-1) - dt.*B);
	end
	
	uvel=sol(ng/2+1:N+ng/2,:);
	vvel=sol(N+3*ng/2+1:2*N+3*ng/2);
	dens=sol(2*N+5*ng/2+1:3*N+5*ng/2);

		  
		  
end

function out=phi(x1,flag)
	global om c nu q k Mp xsoft Lx N

	x = sqrt(xsoft+x1.^2);
	
	if flag==0 % Phi
		out= Mp.*besselk(0,abs(k.*x))./pi;
	end
	
	if flag==1 % Phi'
		out= sign(x1).*Mp.*k.*besselk(1,abs(k.*x))./pi;
	end
	
end	
	
function mat=clearghost(mat,ng)

	for i=1:ng/2
		mat(i,:)=0; 
		mat(end-i+1,:)=0;
	end
	
end
