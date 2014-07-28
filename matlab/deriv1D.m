function out=deriv1D(dx,y,method,order,bc)
% method = center, forward, backward
% order = 2,4,6,8 for center
%		= 1,2,3,4,5,6 for forward/backward
% bc = periodic or zero
% x has uniform spacing

	out = zeros(size(y));
if strcmp(method,'center')
%	disp(['Method is ',method])
%	disp(['Order is ',num2str(order)])
	coeffs = [ 0 0 0 -1/2 ; ...
			   0 0 1/12 2/3; ...
			   0 1/60 3/20 -3/4; ...
			   1/280 -4/105 1/5 -4/5];
	
	coeffs = [ coeffs [0;0;0;0] fliplr(-coeffs)];

	out = zeros(size(y));
% Pad y for boundary terms based b.c's

	if strcmp(bc,'periodic')
%		disp('Boundary condition is periodic')
		newy = [ y(end-3); y(end-2); y(end-1); y(end); y(:); y(1); y(2); y(3); y(4)];
	else
		if strcmp(bc,'zero')
%		disp('Boundary condition not specified, defaulting to zero values');
		newy = [ 0; 0; 0; 0; y(:); 0; 0; 0; 0];
		else
		newy = [1;1;1;1;y(:);1;1;1;1];
		end
	end
	

	for j=5:length(y)+4
		out(j-4)=sum(coeffs(order/2,:).*newy(j-4:j+4)')./dx;
	end

elseif strcmp(method,'forward')
%	disp(['Method is ',method])
%	disp(['Order is ',num2str(order)])
	coeffs = [ -1 1 0 0 0 0 0; ...
			   -3/2 2 -1/2 0 0 0 0; ...
			   -11/6 3 -3/2 1/3 0 0 0; ...
			   -25/12 4 -3 4/3 -1/4 0 0; ...
			   -137/60 5 -5 10/3 -5/4 1/5 0; ...
			   -49/20 6 -15/2 20/3 -15/4 6/5 -1/6];
		
	
	out = zeros(size(y));
% Pad y for boundary terms based b.c's

	if strcmp(bc,'periodic')
%		disp('Boundary condition is periodic')
		newy = [ y(:); y(1:6)];
	else
%		disp('Boundary condition not specified, defaulting to zero values');
		newy = [ y(:); 0; 0; 0; 0; 0; 0];
	end
	
	

	for j=1:length(y)
		out(j)=sum(coeffs(order,:).*newy(j:j+6)')./dx;
	end

elseif strcmp(method,'backward')
%	disp(['Method is ',method])
%	disp(['Order is ',num2str(order)])
	coeffs = [ -1 1 0 0 0 0 0; ...
			   -3/2 2 -1/2 0 0 0 0; ...
			   -11/6 3 -3/2 1/3 0 0 0; ...
			   -25/12 4 -3 4/3 -1/4 0 0; ...
			   -137/60 5 -5 10/3 -5/4 1/5 0; ...
			   -49/20 6 -15/2 20/3 -15/4 6/5 -1/6];
		
	coeffs = - coeffs;
	
	out = zeros(size(y));
% Pad y for boundary terms based b.c's

	if strcmp(bc,'periodic');
%		disp('Boundary condition is periodic')
		newy = [y(end); y(end-1); y(end-2); y(end-3); y(end-4); y(end-5); y(:)];
	else
%		disp('Boundary condition not specified, defaulting to zero values');
		newy = [ 0; 0; 0; 0; 0; 0; y(:)];
	end
	

	for j=7:length(y)+6
		out(j-6)=sum(coeffs(order,:).*newy(j-6:j)')./dx;
	end
elseif strcmp(method,'combo')
%	disp(['Method is ',method])
%	disp(['Order is ',num2str(order)]);
	
	coeffsf = [ -1 1 0 0 0 0 0; ...
			   -3/2 2 -1/2 0 0 0 0; ...
			   -11/6 3 -3/2 1/3 0 0 0; ...
			   -25/12 4 -3 4/3 -1/4 0 0; ...
			   -137/60 5 -5 10/3 -5/4 1/5 0; ...
			   -49/20 6 -15/2 20/3 -15/4 6/5 -1/6];
	coeffsc = [ 0 0 0 -1/2 ; ...
			   0 0 1/12 2/3; ...
			   0 1/60 3/20 -3/4; ...
			   1/280 -4/105 1/5 -4/5];
	
	coeffsc = [ coeffsc [0;0;0;0] fliplr(-coeffsc)];

	coeffsb = - coeffsf;	

	if mod(order,2)==0
		orderc=order;
	else
		orderc=order-1;
	end

	Ny = length(y);
	disp(Ny)
	for j=1:Ny
	disp(j)
		if (j<Ny/3+1)
			out(j) = sum(coeffsf(order,:).*y(j:j+6)')./dx;
		else
			if (j<2*Ny/3+1)
			out(j) = sum(coeffsc(orderc/2,:).*y(j-4:j+4)')./dx;
			else
				out(j) = sum(coeffsb(order,:).*y(j-6:j)')./dx;
			end
		end
	end 

else
	disp('Not a valid method')
	out=0;
end
end
