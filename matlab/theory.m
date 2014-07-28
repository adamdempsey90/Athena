function [x,u,v,d]=theory(k)

	global om c nu q Mp xsoft

	v=




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