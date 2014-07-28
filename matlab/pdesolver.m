function [uvel,vvel,dens,x,time]=pdesolver(N,Lx,k,dt,endt,restartflag,ics)
% 
%
% (I + A dt ) X^(n+1) = X^n - dt B
%
%
	global om c nu q Mp xsoft
	 
	om=1; c=1; nu=1e-1; q=1.5; xsoft=1e-2; Mp=1;

	h = Lx/(N-1);
	x = -Lx/2 + (0:N-1)'.*h;
	
% Make A matrix
%	Make Derivative Matrices using 4th order finite differences

%	Dc = [1/12 -2/3 0 2/3 -1/12];
%	D2c = [-1/12 4/3 -5/2 4/3 -1/12];

	Dc = [-1/2 0 1/2];
	D2c = [1 -2 1];
	
	ng = (length(Dc)-1);
	
	ghostind = [ 1:ng/2, 1+N+ng/2:N+3*ng/2, 1+2*N+3*ng/2:2*N+5*ng/2, 1+3*N+5*ng/2:3*(N+ng)]; 
	
	nghostind = setdiff(1:(3*(N+ng)),ghostind);
	
	xvec=zeros(3*(N+ng),1);
	
	xvec(ng/2+1:N+ng/2)=x(:);
	xvec(1+N+3*ng/2:2*N+3*ng/2)=x(:);
	xvec(2*N+5*ng/2+1:3*N+5*ng/2)=x(:);
	xvec=xvec.*-q*om*1.0i*k;
	xmat=diag(xvec);
	
	
	A = zeros(3*(N+ng),3*(N+ng));
	B = zeros(3*(N+ng),1);
	D = zeros(N+ng,N+ng);		% 2 ghost zones on either boundary. 
	D2 = zeros(N+ng,N+ng);
	ident = eye(N+ng);
	zmat = zeros(N+ng,N+ng);
	
	for j=1:length(Dc)
% Room to change the stencil for higher/lower order methods
% Note # of ghost zones has to change if order increased.
		D += diag(Dc(j).*ones(N+ng-abs(j-ng/2-1),1),j-ng/2-1);
		D2 += diag(D2c(j).*ones(N+ng-abs(j-ng/2-1),1),j-ng/2-1);
	end
	D /= h; D2 /= h^2;

	ident=clearghost(ident,ng); 
	D=clearghost(D,ng);
	D2=clearghost(D2,ng);
	pmat=zmat; dpmat=zmat;
	B(1+ng/2:N+ng/2) = -phi(x,k,1);
	B(N+1+3*ng/2:2*N+3*ng/2) = -1.0i.*k.*phi(x,k,0);
	B(ghostind)=0;

	A = [ nu.*(D2 - k.^2.*ident) , 2.*om.*ident, -D ; ...
		  -2.*(1-.5.*q).*om.*ident, nu.*(D2 - k.^2.*ident) , -1.0i.*k.*ident; ...
		  c.^2.*D, 1.0i.*k.*c^2.*ident, zmat];
	A -= xmat;
		  
% Set up time stepping and solution matrix

	nt = endt/dt + 1;
	time = (0:nt-1)'.*dt;
	sol = zeros(3*(N+ng),nt);
	ident2 = clearghost(eye(3*(N+ng)),ng); 
	
	bcvals = zeros(length(ghostind),1);

	bcvals(2*ng+1:3*ng)=1;
	sol(ghostind,1)=bcvals(:);
%	disp(sol(ghostind,1));
%	disp(ident2(ghostind,1))

if restartflag==1
% ics in [ u; v; d ] order 
	sol(nghostind,1)=ics(:);
end



	for t=2:nt
		disp(['Step ',num2str(t-1),' of ',num2str(nt),', time=',num2str(time(t))]);	
		sol(:,t) = (ident2 - A.*dt/2) \ ((ident2+A.*dt/2)*sol(:,t-1) + 2.*B);
		sol(ghostind,t)=sol(ghostind,1);
%		disp(sol(ghostind,t));
	end
	
	uvel=sol(ng/2+1:N+ng/2,:);
	vvel=sol(N+3*ng/2+1:2*N+3*ng/2,:);
	dens=sol(2*N+5*ng/2+1:3*N+5*ng/2,:);

	figure; plot(x,real(uvel(:,end)),'-b',x,imag(uvel(:,end)),'-r'); 
	legend('Re(u)','Im(u)');
	figure; plot(x,real(vvel(:,end)),'-b',x,imag(vvel(:,end)),'-r'); 
	legend('Re(u)','Im(u)');
	figure; plot(x,real(dens(:,end)),'-b',x,imag(dens(:,end)),'-r'); 
	legend('Re(\delta\Sigma)','Im(\delta\Sigma)');	  
		  
end

function out=phi(x1,k,flag)
	global om c nu q Mp xsoft
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
