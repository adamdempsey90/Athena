function out=calcderiv(dx,y,order)
% method = center, forward, backward
% order = 2,4,6,8 for center
%		= 1,2,3,4,5,6 for forward/backward
% bc = periodic or zero
% x has uniform spacing

	coeffs = [ 0 0 0 -1/2 ; ...
			   0 0 1/12 2/3; ...
			   0 1/60 3/20 -3/4; ...
			   1/280 -4/105 1/5 -4/5];
	
	coeffs = [ coeffs [0;0;0;0] fliplr(-coeffs)];

	out = zeros(size(y));
% Pad y for boundary terms based b.c's

	newy = [y(5); y(4); y(3); y(2); y(:); y(end-2); y(end-3); y(end-4); y(end-5)];



	for j=5:length(y)+4
		out(j-4)=sum(coeffs(order/2,:).*newy(j-4:j+4)')./dx;
	end

end
