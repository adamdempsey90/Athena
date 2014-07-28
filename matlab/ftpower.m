function [kh,out]=ftpower(x,ft)


Nx=length(x);
nk = length(ft(:,1));

NFFT = 2.^nextpow2(nk);
k=[0:NFFT/2+1 -1*(NFFT/2:-1:1)];

%k = [ nk/2:-1:1 0 1:floor((nk-1)/2)];
%k = [0:nk/2 -1.*(floor((nk-1)/2):-1:1)];
kh = fftshift(k);

power=zeros([nk 1]);

for i=1:nk
	power(i) = trapz(x',conj(ft(i,:)).*ft(i,:))./(max(x)-min(x));
end
out=fftshift(power);
figure; semilogy(out./power(1,:),'-x');
ylabel('Integrated Power/Zero Mode');
