#include "athena.h"

void calc_deriv1D(GridS *, double ***, double ***, int, int, int); 

Real **Fh, **Fv, **Fx, **Th;

int main(void) {
	
	GridS *pG;
	double ***in, ***out;
	
	calc_deriv1D(pG, in, out, 1,0,2);
	
	return 0;
}
void init_yavg_output(void) {
	
	Fh = calloc_2D(sizeof(Real));
	Fv = calloc_2D(sizeof(Real));
	
	...
	
	
	return;
}
void calc_yavg_output(GridS *pG, Prims *pV) {
	int i,j,k;
	
	for (k=pG->ks; k<=pG->ke; k++) {
		for (i=pG->is; i <= pG->ie; i++) {
			cc_pos( je and js, x1)
			Fh[k][i]=.5*((pV->d[k][js][i])*(pV->V1[k][js][i])*(pV->V2[k][js][i]) 
						+ (pV->d[k][je][i])*(pV->V1[k][je][i])*(pV->V2[k][je][i]));
			Fx[k][i]=.5*( .... )
			Fv[k][i]=.5*(pG->d[k][js][i]*(pG->V2[k][js][i-1]-pG->V2[k][js][i+1])
						+pG->d[k][je][i]*(pG->V2[k][je][i-1]-pG->V2[k][je][i+1]));
			Th[k][i]=.5*(pG->d[k][js][i]*dyPotplanet(k,js,i) 
						+ pG->d[k][je][i]*dyPotplanet[k,je,i]);
	
			for (j=pG->js+1; j< pG->je; j++) {
				cc_pos( j,x1 )
				Fh[k][i] += (pV->d[k][j][i])*(pV->V1[k][j][i])*(pV->V2[k][j][i]);
				Fx[k][i] += ( ... );
				Fv[k][i] += pG->d[k][j][i]*(pG->V2[k][j][i-1]-pG->V2[k][j][i+1]);
				Th[k][i] += pG->d[k][j][i]*dyPotplanet(k,j,i);
			}
			
			Fh *= pG->dx2;
			Fx *= .5*omega*(pG->dx2);
			Fv *= nu_iso * (pG->dx2) * 0.5 / (pG->dx1);
			Th *= pG->dx2;
			
				
		}
	}
	
// Output
	
	for(k=pG->ks; k<=pG->ke; k++) {
		for(i=pG->is;i<=pG->ie;i++) {
			cc_pos( j k x1 x2 x3)
			fprintf(outfile,"%g %g %g %g %g %g	\n", &x1,Fh[k][i],Fx[k][i],Fv[k][i],Th[k][i]);
		}
	}
	
	return;
}