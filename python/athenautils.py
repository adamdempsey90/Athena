class Data():
	def __init__(self,dir,time,writeflag=False):

		self.time = time
		if time < 10:
			tstr='000'+str(time)
		elif time < 100:
			tstr = '00' + str(time)
		elif time < 1000:
			tstr = '0' + str(time)
		else:
			tstr = str(time)
		
		if dir[-1] != '/':
			dir += '/'

		print "Loading all files at output # "+tstr+"from" + dir +" ..."
		procstr = "ls "+dir+" | grep id | wc -l" 

		try:
			self.np = int(subprocess.check_output(procstr,shell=True)) 
		except UnboundLocalError:
			import subprocess
			self.np = int(subprocess.check_output(procstr,shell=True)) 

		if self.np == 0:
			print "Working in serial..."
			dat = loadtxt(dir+"PlanetDisk."+tstr+".tab")
		else:
			print "Working with ",self.np," cores..."
			for p in range(self.np):
				if p==0:
					fname=dir+'id0/PlanetDisk.'+tstr+'.tab'
					temp=subprocess.check_output('grep Time '+fname,shell=True)
					self.time=float([temp.split()[i+2] for i in range(len(temp.split())) if temp.split()[i]=='Time'][0][:-1])
					print 'This Output Time is ',self.time
					print '\tLoading '+fname
					dat=loadtxt(fname)
				else:		
					fname=dir+'id'+str(p)+'/PlanetDisk-id'+str(p)+'.'+tstr+'.tab'
					print '\tLoading '+fname
					dat=vstack((dat,loadtxt(fname)))
		self.raw = dat
		print 'Cleaning the Data...'
		self.x,self.y,self.dens,self.vx,self.vy,self.vort,self.nx,self.ny = cleandata(dat)
		self.vxbar,self.vybar,self.dbar,self.u,self.v,self.sig = calcMeans(self,smoothed=True)
		self.Lx = max(self.x[0,:])-min(self.x[0,:])
		self.Ly = max(self.y[:,0])-min(self.y[:,0])
		self.dx = self.Lx/self.nx
		self.dy = self.Ly/self.ny
# Sloppy way of doing this...
		print 'Loading the parameters...'


		temp=dir.split('m')[-1].split('nu')
		self.mp = float(temp[0])
		self.nu = float(temp[-1].split('/')[0].split('_')[0])
		

#		self.mp =  float(subprocess.check_output("grep Mplanet "+dir+"athinput.new",shell=True).split()[-1])
#		self.nu =  float(subprocess.check_output("grep nu_iso "+dir+"athinput.new",shell=True).split()[-1])
		self.om = float(subprocess.check_output("grep omega "+dir+"athinput.new",shell=True).split()[-1])
		self.xs = float(subprocess.check_output("grep Rsoft "+dir+"athinput.new",shell=True).split()[-1])
		self.c = float(subprocess.check_output("grep iso_csound "+dir+"athinput.new",shell=True).split()[-1])
		self.q = float(subprocess.check_output("grep qshear "+dir+"athinput.new",shell=True).split()[-1])
#		self.xx,self.yy = meshgrid(self.x,self.y)
		self.phi = -self.mp/sqrt(self.x**2 + (self.xs*ones(self.x.shape))**2 + self.y**2)
		self.vortens = (self.vort + 2*self.om*ones(self.dens.shape))/self.dens
		self.vortbar = self.vortens.mean(axis=0)
		self.xi = self.vortens -  tile(self.vortbar,[dat.ny,1])
		if (writeflag):
			print 'Dumping data'
			writedata(self,tstr,dir)


		
class FTData():
	def __init__(self,dat):
		self.xx = dat.x
		self.x = dat.x[0,:]
		self.nu = dat.nu
		self.q = dat.q
		self.mp = dat.mp
		self.om = dat.om
		self.xs = dat.xs
		self.c = dat.c
	
		if ~mod(dat.ny,2):
			norm = dat.ny/2+1
		else:
			norm = (dat.ny+1)/2
#		temp=fft.rfft(dat.vx.transpose()).shape
#		norm = dat.ny
		self.vxhat = fft.rfft(dat.vx.transpose())/norm
		self.vyhat = fft.rfft(dat.vy.transpose())/norm
		self.dhat = fft.rfft(dat.dens.transpose())/norm
		self.k = fft.fftfreq(dat.vx.shape[0])[:norm]

		self.ik = self.k / self.k[1]
		self.k = self.ik * 2*pi/dat.Ly
		self.kk = tile(self.k,[dat.nx,1])
		self.ikk = self.kk / self.k[1]
		
				
		self.sk = fft.fftshift(self.k)
		self.sik = fft.fftshift(self.ik)
		self.skk = fft.fftshift(self.kk)
		self.sikk = fft.fftshift(self.ikk)
		self.svxhat = fft.fftshift(self.vxhat)
		self.svyhat = fft.fftshift(self.vyhat)
		self.sdhat = fft.fftshift(self.dhat)

		Phi = -dat.mp/sqrt(dat.x*dat.x + dat.xs*dat.xs+dat.y*dat.y)
		self.phi = fft.rfft(Phi.transpose())/norm
		self.sphi = fft.fftshift(self.phi)

def writedata(dat,time,dir):	
	with open(dir+'m'+str(dat.mp)+'_nu'+str(dat.nu)+'_dump_'+str(time)+'.dat',"w") as f:
		for line in dat.raw:
			line = [str(i)+'\t' for i in list(line)]
			f.write(''.join(line))
			f.write('\n')


		
def cleandata(dat):
	xlist = unique(dat[:,2])
	ylist = unique(dat[:,3])
	Nx = len(xlist)
	Ny = len(ylist)
	
	fx = dat[:,2]
	fy = dat[:,3]
	fd = dat[:,4]
	fvx = dat[:,5]
	fvy = dat[:,6]

	x,y=meshgrid(xlist,ylist)
# x[0,:] gives x 		y[:,0] gives y
	dens=zeros(x.shape)
	vx = zeros(x.shape)
	vy = zeros(x.shape)
	for i,j in enumerate(xlist):
		dens[:,i] = fd[fx==j]
		vx[:,i] = fvx[fx==j]
		vy[:,i] = fvy[fx==j]
	
#	vy += 1.5*x
	_,dyvx = gradient(vx,mean(diff(ylist)),mean(diff(xlist)))
	dxvy,_ = gradient(vy,mean(diff(ylist)),mean(diff(xlist)))
	vort = dxvy - dyvx	
	return x,y,dens,vx,vy,vort,Nx,Ny 

def calcMeans(dat,smoothed=True):
	
	if(smoothed):
		vxbar = tile(smooth(dat.vx.mean(axis=0)),[dat.ny,1])
		vybar = tile(smooth(dat.vy.mean(axis=0)+1.5*dat.x[0,:])-1.5*dat.x[0,:],[dat.ny,1])
 		dbar = tile(smooth(dat.dens.mean(axis=0)),[dat.ny,1])
	else:
		vxbar = tile((dat.vx.mean(axis=0)),[dat.ny,1])
		vybar = tile((dat.vy.mean(axis=0)+1.5*dat.x[0,:])-1.5*dat.x[0,:],[dat.ny,1])
 		dbar = tile((dat.dens.mean(axis=0)),[dat.ny,1])


	u = dat.vx - vxbar
	v = dat.vy - vybar
	sig = dat.dens - dbar

	return vxbar,vybar,dbar,u,v,sig
def pspec(dat,var='density'):
	k=fft.fftshift(dat.ik)
	
	if var=='density':
		QQ = dat.dhat*conj(dat.dhat)
		ystr = '< \sigma|^2 >_x'
	elif var=='vx':
		QQ = dat.vxhat*conj(dat.vxhat)
		ystr = '< |u|^2 >_x'
	elif var=='vy':
		QQ = dat.vyhat*conj(dat.vyhat)
		ystr = '< |v|^2 >_x'
	elif var=='vxvy':
		QQ = dat.vxhat*conj(dat.vyhat)
		ystr = '< u v^* >'
	elif var=='dvxvy':
		QQ = dat.dhat*conj(dat.vxhat)*dat.vyhat
		ystr = '< sig u^* v >'

	Q = fft.fftshift(QQ).mean(axis=0)
	figure()
	semilogy(k[(k>=0)&(Q>1e-6)],Q[(k>=0)&(Q>1e-6)],'-x')
	xlabel('n (k=n*2*pi/Ly)')
	ylabel(ystr)

		
def calcFluxes(dat):

	xfiglim = 6*(dat.c/dat.om)
	tstr = 'm='+str(dat.mp)+', nu='+str(dat.nu)+', res='+str(dat.nx)+'x'+str(dat.ny)
#	dyd,dxd = gradient(dat.dens,dat.dy,dat.dx)
	dyvx,dxvx = gradient(dat.vx,dat.dy,dat.dx)
	dyvy,dxvy = gradient(dat.vy,dat.dy,dat.dx)
	dyphi,dxphi = gradient(dat.phi,dat.dy,dat.dx)
	
	_,dxFtot = gradient(dat.dens*dat.vx*(dat.vy+2*dat.om*dat.x)-dat.nu*dat.dens*(dyvx+dxvy),dat.dy,dat.dx)

	dxFtot = dxFtot.mean(axis=0)
	Ptot = (dat.dens*dyphi).mean(axis=0)

	x = dat.x[0,:]
	ind = (x>=-xfiglim) & (x<=xfiglim)
	vxbar = dat.vxbar 

	vybar = dat.vybar
	dbar = dat.dbar

	u = dat.u 
	v = dat.v
	sig = dat.sig

	_,dxFnlin = gradient(sig*u*v,dat.dy,dat.dx)
	dxFnlin = dxFnlin.mean(axis=0)

	_,dxFlin = gradient(sig*vxbar*vybar+dbar*vxbar*v+u*vybar*dbar,dat.dy,dat.dx)
	dxFlin = dxFlin.mean(axis=0)

	dysig,dxsig = gradient(sig,dat.dy,dat.dx)
	dyu,dxu = gradient(u,dat.dy,dat.dx)
	dyv,dxv = gradient(v,dat.dy,dat.dx)

	_,dxdbar = gradient(dbar,dat.dy,dat.dx)
	_,dxvxbar = gradient(vxbar,dat.dy,dat.dx)
	_,dxvybar = gradient(vybar,dat.dy,dat.dx)


	dy2sig,dxdysig= gradient(dysig,dat.dy,dat.dx)
	_,dx2sig = gradient(dxsig,dat.dy,dat.dx)

	dy2u,dxdyu= gradient(dyu,dat.dy,dat.dx)
	_,dx2u = gradient(dxu,dat.dy,dat.dx)

	dy2v,dxdyv= gradient(dyv,dat.dy,dat.dx)
	_,dx2v = gradient(dxv,dat.dy,dat.dx)

	_,dx2dbar = gradient(dxdbar,dat.dy,dat.dx)
	_,dx2vxbar = gradient(dxvxbar,dat.dy,dat.dx)
	_,dx2vybar = gradient(dxvybar,dat.dy,dat.dx)

	Pibarxy = dat.nu*dbar*dxvybar

	Pipxy = dat.nu*(dbar*(dxv+dyu)+sig*dxvybar)
	Pipyy = -dat.c*dat.c*sig + (2.0/3.0)*dat.nu*(dbar*(2*dyv-dxu) - sig*dxvxbar)
	Pippxy = dat.nu*sig*(dxv+dyu)



	divPip = dxdbar*(dxv+dyu)+dbar*(dx2v+dxdyu)+dxsig*dxvybar+sig*dx2vybar
	divPip += dbar*(2*dy2v-(2.0/3.0)*(dxdyu+dy2v))+dysig*(-2.0/3.0)*dxvxbar
	divPip *= dat.nu
	divPip -= dat.c*dat.c*dysig

	divPibar = dat.nu*(dxdbar*dxvybar+dbar*dx2vybar)

	divPipp = dat.nu*(dxsig*(dxv+dyu) + sig*(dx2v + dxdyu))

#	_,divPibar = gradient(Pibarxy,dat.dy,dat.dx)
	_,divPipp = gradient(Pippxy,dat.dy,dat.dx)
#	_,dxPipxy = gradient(Pipxy,dat.dy,dat.dx)
#	dyPipyy,_ = gradient(Pipyy,dat.dy,dat.dx)	
#	divPip = dxPipxy + dyPipyy

	TTwd = -dbar*u*dxv + sig*u*(dxvybar+2*dat.om) + (sig/dbar)**2 * divPibar - (sig/dbar)*divPip

	Twd = TTwd.mean(axis=0)

#	figure();
#	plot(x,(-dbar*u*dxv).mean(axis=0))
#	plot(x,(sig*u*(dxvybar+2*dat.om)).mean(axis=0))
#	plot(x,(-dbar*u*dxv).mean(axis=0)+(sig*u*(dxvybar+2*dat.om)).mean(axis=0))
#	plot(x,((sig/dbar)**2 *divPibar).mean(axis=0))
#	plot(x,(-(sig/dbar)*divPip).mean(axis=0))
#	plot(x,Twd)
#	title('Twd breakdown')
#	legend(('1','2','1+2','3','4','sum'))

	FFp = vxbar*sig*v + dbar*u*v
	Fp=FFp.mean(axis=0)
#	dxFp = gradient(Fp,dat.dx)
	
	dxFp = dxvxbar*sig*v+vxbar*dxsig*v+vxbar*sig*dxv
	dxFp += dxdbar*u*v + dbar*dxu*v + dbar*u*dxv
	dxFp = dxFp.mean(axis=0)

	FFb = (vybar+2*dat.om*dat.x)*(dbar*vxbar+sig*u) - Pibarxy - Pippxy
	Fb=FFb.mean(axis=0)
#	dxFb = gradient(Fb,dat.dx)	
	dxFb = (vybar+2*dat.om*x)*(dxdbar*vxbar + dbar*dxvxbar + dxsig*u+sig*dxu)
	dxFb += (dxvybar + 2*dat.om)*(dbar*vxbar + sig*u)
	dxFb += -divPibar-divPipp
	dxFb = dxFb.mean(axis=0)

#	phi = -dat.mp*pow(dat.x**2 + dat.y**2 + dat.xs**2,-0.5)
#	phip = dat.phi-	tile(dat.phi.mean(axis=0),[dat.ny,1])

	dyPhi,dxPhi = gradient(dat.phi,dat.dy,dat.dx)

#	dyPhi = dat.mp*dat.y*pow(dat.x*dat.x + dat.xs*dat.xs+dat.y*dat.y,-1.5)

#	dyPhi -= tile(dyPhi.mean(axis=0),[dat.ny,1])
	
	TTh = sig*dyPhi
	Th=TTh.mean(axis=0)	

#	figure();
#	plot(x,gradient((vxbar*sig*v).mean(axis=0),dat.dx))
#	plot(x,gradient((dbar*u*v).mean(axis=0),dat.dx))
#	plot(x,dxFp)
#	plot(x,Th)
#	plot(x,dxFp+Th)
#	title('dxFp+Th breakdown')
#	legend(('F1','F2','Fsum','Th','total'))


#	figure(); plot(x,Th,x,-Twd);
#	xlabel('x'); legend(('Th','-Twd'))

	figure(); plot(x,-(dxFp+Th),x,-(dxFp+Th+dxFnlin),x,Twd); xlim([-xfiglim,xfiglim]);
	#ylim([2*min(hstack((Twd[10:-10],dxFp[10:-10]+Th[10:-10]))),2*max(hstack((Twd[10:-10],dxFp[10:-10]+Th[10:-10])))])
#	ylim([-.0005,.0005])
	xlabel('x');

	y0=min(-(dxFp+Th)[ind])
	y0=min([y0,min(-(dxFp+Th+dxFnlin)[ind])])
	y0=min([y0,min(Twd[ind])])
	y1=max(-(dxFp+Th)[ind])
	y1=max([y1,max(-(dxFp+Th+dxFnlin)[ind])])
	y1=max([y1,max(Twd[ind])])

	ylim([y0*1.5,y1*1.5])
	title('WSS in Rotating Frame: '+tstr)
	legend(('dxFp+Th','dxFp + Th + non-linear','Twd'),loc='lower right')

	figure(); plot(x,dxFb,x,Twd); xlim([-xfiglim,xfiglim]);
	y0=min(dxFb[ind])
	y0=min([y0,min(Twd[ind])])
	y1=max(dxFb[ind])
	y1=max([y1,max(Twd[ind])])
	ylim([y0*1.5,y1*1.5])

	xlabel('x'); title('VSS in Rotating Frame: '+tstr)
	legend(('dxFb','Twd'))

#	figure(); plot(x,dxFb+dxFp,x,-Th); xlim([-xfiglim,xfiglim]);
#	xlabel('x'); title('Total Flux Balance: '+tstr)
#	ylim([-.005,.005])
#	legend(('dxFb+dxFp','-Th'))
	
	figure(); plot(x,-(dxFtot-dxFb+Th),x,-(dxFp+dxFnlin+Th),x,Twd);xlim([-xfiglim,xfiglim]);
	y0=min(-(dxFtot-dxFb+Th)[ind])
	y0=min([y0,min(-(dxFp+Th+dxFnlin)[ind])])
	y0=min([y0,min(Twd[ind])])
	y1=max(-(dxFtot-dxFb+Th)[ind])
	y1=max([y1,max(-(dxFp+Th+dxFnlin)[ind])])
	y1=max([y1,max(Twd[ind])])

	ylim([y0*1.5,y1*1.5])
	xlabel('x'); title('WSS from Total Flux: '+tstr)
	legend(('dxFtot-dxFb+Th','dxFp+Th','Twd'),loc='lower right')

	figure(); plot(x,dxFtot-dxFb,'.',x,dxFp,x,dxFnlin,x,dxFp+dxFnlin)
	y0=min((dxFtot-dxFb)[ind])
	y0=min([y0,min(dxFp[ind])])
	y0=min([y0,min(dxFnlin[ind])])
	y0=min([y0,min((dxFp+dxFnlin)[ind])])
	y1=max((dxFtot-dxFb)[ind])
	y1=max([y1,max(dxFp[ind])])
	y1=max([y1,max(dxFnlin[ind])])
	y1=max([y1,max((dxFp+dxFnlin)[ind])])
	ylim([y0*1.5,y1*1.5])
	xlim([-xfiglim,xfiglim])
	xlabel('x')
	title('Non-linear term: '+tstr)
	legend(('Ftot-Fb','Fp','NL','Fp+NL'))

	figure(); plot(x,dxFlin)
	return Twd,dxFp,Th,dxFb,Fp,Fb,dxFtot




def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=numpy.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    y=numpy.convolve(w/w.sum(),s,mode='valid')
    window_len+=1
    return y[(window_len/2-1):-(window_len/2-1)]

def calcFluxesft(ft):
	
	dx = mean(diff(ft.x))
	dk = ft.k[1] - ft.k[0]
	kk = ft.kk[:,1:]
	vx=ft.vxhat
	vy = ft.vyhat
	d = ft.dhat

	print '1'
	vxbar = vx[:,0]
	u = vx[:,1:]
	vybar = vy[:,0]
	v = vy[:,1:]
	dbar = d[:,0]
	sig = d[:,1:]
	print '2'
#	Phi = ft.phi		
#	phibar = Phi[:,0]
#	phi = Phi[:,1:]
	print '3'
	Fp = vxbar*Ksum(sig,v) + dbar*Ksum(u,v)
	dxFp = gradient(Fp,dx) 
	print '4'
	Th = Ksum(sig,ft.phi[:,1:]*1j*kk)
	print '5'

	dxvybar = gradient(vybar,dx)
	dxvxbar = gradient(vxbar,dx)
	dxdbar = gradient(dbar,dx)
	dxv,_ = gradient(v,dx,dk)
	dxu,_ = gradient(u,dx,dk)
	dxsig,_ = gradient(sig,dx,dk)
	print '6'

	xyterm = dxv + 1j*kk*u
	dxxyterm,_ = gradient(xyterm,dx,dk)
	print '7'
	divPip = dxdbar*Ksum(sig,xyterm) + dbar*Ksum(sig,dxxyterm) + Ksum(sig,sig)*dxvybar
	print '7.1'
	divPip += (dbar*Ksum(sig,-2*kk*kk*v-(2.0/3.0)*1j*kk*(dxu+1j*kk*v))-(2.0/3.0)*Ksum(sig,sig*1j*kk)*dxvxbar)
	print '7.2'
	divPip *= ft.nu/dbar
	print '7.3'
	divPip -= ft.c*ft.c*Ksum(sig,sig*1j*kk)/dbar	
	print '8'
	divPibar = ft.nu*gradient(dbar*dxvybar,dx)
	divPibar *= Ksum(sig,sig)/(dbar**2)
	print '9'
	Twd = -dbar*Ksum(u,dxv) + (dxvybar+2*ft.om)*Ksum(sig,u) + divPibar + divPip

	print '10'
	
	Fb = (vybar+2*ft.om*ft.xx)*(dbar*vxbar+Ksum(sig,u))
	Fb -= ft.nu*dbar*dxvybar
	Fb -= ft.nu*Ksum(sig,dxv+1j*kk*u)
	dxFb = gradient(Fb,dx)
	print '11'

	x = ft.x
	
	figure(); plot(x,dxFp+Th,x,-Twd); xlim([x[10],x[-10]]);
#	ylim([2*min(hstack((Twd[10:-10],dxFp[10:-10]+Th[10:-10]))),2*max(hstack((Twd[10:-10],dxFp[10:-10]+Th[10:-10])))])
#	ylim([-.0005,.0005])
	xlabel('x');  title('WSS in Rotating Frame')
	legend(('dxFp+Th','-Twd'))

#	figure(); plot(x,dxFb,x,Twd); xlim([x[10],x[-10]]);
#	ylim([2*min(hstack((-Twd[10:-10],dxFb[10:-10]))),2*max(hstack((-Twd[10:-10],dxFb[10:-10])))])
#	ylim([-.0005,.0005])
#	xlabel('x'); title('VSS in Rotating Frame')
#	legend(('dxFb','Twd'))

#	figure(); plot(x,dxFb+dxFp,x,-Th); xlim([x[10],x[-10]]);
#	xlabel('x'); title('Total Flux Balance')
#	ylim([-.005,.005])
#	legend(('dxFb+dxFp','-Th'))

	return Twd, dxFp, Th, Fp, dxFb, Fb




def Ksum(q1,q2,a=1,b=1):
	
	T = (a*q1.conj()*b*q2).real

	return T.sum(axis=1)


def centderiv(Q,dy,dx,order=2,bc=('periodic','out'),dbc=(0,0)):
	
	ng = order/2
 	
	s=Q.shape
	if Q.ndim==1:
		work=zeros(s[0]+order)
	elif Q.ndim==2:
		work = zeros((s[0]+order,s[1]+order))


	if order==2:
		coefs=array([-.5,.5])
		pos = array([-1, 1])
	elif order==4:
		coefs=array([1./12,-2./3,2./3,-1./12])
		pos = array([-2,-1,1,2])
	elif order==6:
		coefs=array([-1./60, 3./20, -1./4, 1./4, -3./20, 1./60])
		pos = array([-3, -2, -1, 1, 2,3])
	elif order==8:
		coefs=array([1./280, -4./105, 1./5, -4./5, 4./5, -1./5, 4./105, -1./280])
		pos = array([-4, -3, -2, -1, 1, 2, 3, 4])
	elif order==10:
		coefs=array([-2., 25., -150., 600., -2100., 2100., -600., 150., -25., 2.])
		coefs /= 2520
		pos=array([-5, -4, -3, -2, -1, 1, 2, 3, 4, 5])

	if Q.ndim==2:	
		work[ng:-ng,ng:-ng] = Q
	
		if bc[0]=='periodic':
			work[ng:-ng,:ng] = Q[:,-ng:] 
			work[ng:-ng,-ng:] = Q[:,:ng]	
		elif bc[0]=='out':
			work[ng:-ng,:ng] = Q[:,ng-1::-1]
			work[ng:-ng,-ng:] = Q[:,:-ng-1:-1]
		elif bc[0]=='Dirichlet':
			work[ng:-ng,:ng] = dbc[0]
			work[ng:-ng,-ng:] = dbc[0]

	
		if bc[1]=='periodic':
			work[:ng,ng:-ng] = Q[-ng:,:] 
			work[-ng:,ng:-ng] = Q[:ng,:]	
		elif bc[1]=='out':
			work[:ng,ng:-ng] = Q[ng-1::-1,:]
			work[-ng:,ng:-ng] = Q[:-ng-1:-1,:]
		elif bc[1]=='Dirichlet':
			work[:ng,ng:-ng] = dbc[1]
			work[-ng:,ng:-ng] = dbc[1]

		dyQ = zeros(s)
		dxQ = zeros(s)

		for j,c in enumerate(coefs):
			i = pos[j]	
			if -ng+i == 0:
				dyQ+=work[i+ng:,ng:-ng]*c
				dxQ+=work[ng:-ng,i+ng:]*c
			else:
				dyQ+=work[i+ng:-ng+i,ng:-ng]*c
				dxQ+=work[ng:-ng,i+ng:-ng+i]*c


		dyQ /= dy
		dxQ /= dx

	elif Q.ndim==1:
		work[ng:-ng] = Q
	
		if bc=='periodic':
			work[:ng] = Q[-ng:] 
			work[-ng:] = Q[:ng]	
		elif bc=='out':
			work[:ng] = Q[ng-1::-1]
			work[-ng:] = Q[:-ng-1:-1]
		elif bc=='Dirichlet':
			work[:ng] = dbc
			work[-ng:] = dbc

	
		dyQ = zeros(s)
		dxQ = zeros(s)

		for j,c in enumerate(coefs):
			i = pos[j]	
			if -ng+i == 0:
				dyQ+=work[i+ng:]*c
			else:
				dyQ+=work[i+ng:-ng+i]*c

		dyQ /= dy



	return dyQ,dxQ

def write_mat(fname,mat):
	with open(fname,'w') as f:
		for line in mat:
			for x in line:
				f.write(str(x)+'\t')
			f.write('\n')

def write_means(fname,dat,interp_flag=False,imult=2):
	dbar= dat.dbar[0,:]
	vxbar = dat.vxbar[0,:]
	vybar = dat.vybar[0,:]+1.5*dat.x[0,:]
	x = dat.x[0,:]
	
	if interp_flag:
		nx = imult*len(x)
		xn = linspace(x[0],x[-1],nx)
		dbar = interp(xn,x,dbar)
		vxbar = interp(xn,x,vxbar)
		vybar = interp(xn,x,vybar)
		x = xn

	out = vstack((x,dbar,vxbar,vybar))
	out = out.transpose()
	write_mat(fname,out)


def write_k(fname,ft,ind):
	
	out = vstack((ft.x,real(ft.vxhat[:,ind]),imag(ft.vxhat[:,ind]), \
			real(ft.vyhat[:,ind]),imag(ft.vyhat[:,ind]), \
			real(ft.dhat[:,ind]), imag(ft.dhat[:,ind])))
	out = out.transpose()
	write_mat(fname,out)



