densarg=1;
fluxarg=0;
loadarg=0;

%direct = {'../bin/256_q3_nu0/','../bin/256_q3_nu001/','../bin/256_q3_nu01/', ...
%		  '../bin/256_q3_nu1/'};
%direct = {'../bin/512_q3_nu0/','../bin/512_q3_nu001/','../bin/512_q3_nu01/', ...
%		  '../bin/512_q3_nu1/'};

direct={'../bin/'};
nd=length(direct);
%labels = {'\nu = 0'; '\nu = 0.001'; '\nu = 0.01'; '\nu = 0.1'};
%labels = {'\nu = 0.0';'\nu = 0.001'; '\nu = 0.01'; '\nu = 0.1'};
labels = {'\nu = 0.1'};
%linetype={'-k', '-b', '-m','-r'};
%linetype={'-k','-b', '-m','-r'};
linetype = {'-k'};
if (loadarg == 0) 
n=length(dir([direct{1},'*.flux.tab']));
files=cell(n,nd); dat=cell(n,nd);
for k=1:nd

base=[direct{k},'PlanetDisk.'];
suff='.flux.tab';
for i=0:n-1
	if i<10
		files{i+1,k} = [base,'000',num2str(i),suff];
	elseif i<100
		files{i+1,k} = [base,'00',num2str(i),suff];
	elseif i<1000
		files{i+1,k} = [base,'0',num2str(i),suff];
	else
		files{i+1,k} = [base,num2str(i),suff];
	endif
	
	dat{i+1,k} = load(files{i+1,k});
	s = size(dat{i+1,k});
	dat{i+1,k}(:,end+1)=zeros(s(1),1);
	ind = find(dat{i+1,k}(:,1) >= 0); ind=ind(2:end);
	dat{i+1,k}(ind,end)=cumtrapz(dat{i+1,k}(ind,1),dat{i+1,k}(ind,5));
	ind = flipud(find(dat{i+1,k}(:,1) < 0)); ind=ind(2:end);
	dat{i+1,k}(ind,end)=cumtrapz(dat{i+1,k}(ind,1),dat{i+1,k}(ind,5));
end
end
end

if (densarg==0)
f=0;
while f<3          
	for i=1:n
	for k=1:nd
		plot(dat{i,k}(:,1),dat{i,k}(:,8),linetype{k}); 
%		plot(dat{i,k}(:,1),dat{i,k}(:,4),'-b',dat{i,k}(:,1),dat{i,k}(:,2),'-r', ... 
%		dat{i,k}(:,1),dat{i,k}(:,11),'-k',dat{i,k}(:,1),dat{i,k}(:,2)+dat{i,k}(:,4),'-m'); 
	if k==1
		hold on;
	end
	end
		ylim([.90,1.05]); %xlim([0,20]);
		legend(labels);
		title(['t = ',num2str(i)]); xlabel('x'); ylabel('\Sigma - \Sigma_0');
		hold off;
		pause(.1)

	
	end
	f+=1;
end
end

if (fluxarg==0)
f=0; 
k=1;
while f<3          
	for i=1:n	
		plot(dat{i,k}(:,1),dat{i,k}(:,4),'-b',dat{i,k}(:,1),dat{i,k}(:,2),'-r', ... 
		dat{i,k}(:,1),dat{i,k}(:,11),'-k',dat{i,k}(:,1),dat{i,k}(:,2)+dat{i,k}(:,4),'-m'); 
	
	
	%	ylim([-.0002,.0002]); 
		xlim([0,25]);
		legend('F_{\nu}','F_H','T_h','F_{\nu}+F_H');
		title(['t = ',num2str(i)]); xlabel('x'); ylabel('\Sigma - \Sigma_0');
		pause(.1)
	end
		f+=1;
end
end