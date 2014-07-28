function out=potential(x,y,rs,m,d)
	if strcmp(d,'y')
		out = m .* y ./ (x.^2 + y.^2+rs).^(1.5); 
	elseif strcmp(d,'x')
		out = m.*x./(x.^2+y.^2+rs).^(1.5);
	else
		out = -m ./ sqrt(x.^2+y.^2+rs);
end
