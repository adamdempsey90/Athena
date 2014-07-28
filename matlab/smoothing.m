
vx=load('vx.dat');
dxvx=load('dxvx.dat');
outvx=sgolayfilt(vx,2,101);
outdxvx=sgolayfilt(dxvx,2,101);

dlmwrite('smoothedvx.dat',outvx,'delimiter','\t');
dlmwrite('smootheddxvx.dat',outdxvx,'delimiter','\t');


exit
