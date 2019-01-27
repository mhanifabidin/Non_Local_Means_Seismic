/* -------------------------------------------------------
/* Matlab C/mex routine for 2D Non Local Means filtering 
/*  
/* t:	radius of the search window                       
/* f:	radius of the averaging window 
/* h:   weight parameter	
/* m,n    array size in x,y directions 							  
/* -------------------------------------------------------*/
#include "mex.h"
#include <math.h> 
#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))
	/* computational subroutine */     
 void mexnlm2D( double *output, double *input, int t, int f, double h, 
		  int m, int n){ 
		  int i,j;
		  int m2,n2;
		  int i1,j1; 
		  int r,s;
		  int r1,s1;
		  int rmin, rmax, smin, smax;
		  int wsiz, wind;	 
		  double dd, w;
		  double sweight, average, wmax; 
		  double *W1, *W2, *W3; 
		  double *zkernel;
		  double *input2; 
		  void make_kernel(double *zkernel, int f); 
		  void pad_input(double *input, double *input2, int m, int n, int f); 
		  
		
		  m2 = m+(2*f); 
		  n2 = n+(2*f); 
		  wsiz = (2*f)+1;
		
		  input2 = mxCalloc(m2*n2,sizeof(double));
		  zkernel= mxCalloc(wsiz*wsiz,sizeof(double)); 
		  W1     = mxCalloc(wsiz*wsiz,sizeof(double)); 
		  W2     = mxCalloc(wsiz*wsiz,sizeof(double)); 
		  W3     = mxCalloc(wsiz*wsiz,sizeof(double)); 		  
         	      
          pad_input(input,input2,m,n,f); 
		 
		  make_kernel(zkernel,f); 
				 
		  /* Main loops */ 
        
          for(i=0;i<m;i++){ 			 
		      for(j=0;j<n;j++){ 
				  
					  i1 = i + f; 
					  j1 = j + f; 
					 
			   
         
					  /* fill window */
		
					  for(r1=0;r1<wsiz;r1++){ 
						  for(s1=0;s1<wsiz;s1++){
						
					  wind = (r1*wsiz)+s1; 
					  W1[wind] = 
						 input2[(i1+r1-f)*n2+j1+s1-f]; 	
					 
				
							  }
						  }
					
					   wmax    = 0.0; 
		               average = 0.0; 
		               sweight = 0.0; 

					   
                       rmin = max(i1-t,f);
					   rmax = min(i1+t,m+f-1);
					   smin = max(j1-t,f); 
					   smax = min(j1+t,n+f-1); 

					   for(r=rmin;r<rmax+1;r++){
						   for(s=smin;s<smax+1;s++){ 
							  
		              dd = 0.0; 
                      w  = 0.0; 
								  
						if(r==i1 && s==j1) continue; 				
								   					
				          for(r1=0;r1<wsiz;r1++){ 
						     for(s1=0;s1<wsiz;s1++){
							  
								   wind = r1*wsiz+s1;
								  
					  W2[wind] =
						  input2[(r+r1-f)*n2+s+s1-f];

						  }
					  }
						 

						  
						 for(r1=0;r1<wsiz;r1++){ 
						     for(s1=0;s1<wsiz;s1++){
						 wind =  r1*wsiz+s1; 
						 W3[wind] = zkernel[wind]*((W1[wind]-W2[wind])*(W1[wind]-W2[wind])); 
                         dd = dd + W3[wind];     
							
						  }
					  }

					  w = exp(-dd/h); 

						
						 if(w>wmax) wmax = w; 

						 sweight = sweight+w; 
						 average = average+(w*input2[r*n2+s]); 
							   }
						   }
					 
							average = average + wmax*input2[i1*n2+j1]; 
							sweight = sweight + wmax; 
				
                      
					   if(sweight>0.00000001) output[i*n+j] = average/sweight; 
					   else output[i*n+j] = input[i*n+j]; 
					  
 
				  }
			  }
		
         
		  return;
	 }



 void make_kernel(double *zkernel,int f)
{
	double aa,cc, dd, ff;
	double value; 
	int i,j,k,d; 
	int m;


		m = 2*f + 1; 
		ff = (double) f;
		for(d=0;d<m*m;d++){
			zkernel[d] = 0.; 
		}

 for(d=1;d<f+1;d++){
           dd = (double) d; 
	        value = 1./((2*dd+1)*(2*dd+1));
		for(i=-d;i<d+1;i++){
			for(j=-d;j <d+1;j++){ 	
					zkernel[(f-i)*m + (f-j)] = 
				    zkernel[(f-i)*m + (f-j)] + (value/ff); 
			}	
			
			}
	}
     aa = 0;
    for(k=0;k<m*m;k++){
		aa = aa+zkernel[k]; 
	} 

    for(d=0;d<m*m;d++){
		zkernel[d] = zkernel[d]/aa; 
			cc = zkernel[d]; 
			
	}   
     

return;  
 }


 void pad_input(double *input, double *input2, int m, int n, int f)
{
     int m2,n2;
	 int ico; 
	 int i,j; 
	 int iw1,iw2,wind;

	 m2 = m + 2*f;
	 n2 = n + 2*f; 
  
	  for(ico=0;ico<m2*n2;ico++){
			  input2[ico] = 0.; 
		  }    		  

		 /* fill padded array */ 
		  	  
		  for(i=0;i<m;i++){ 
			  for(j=0;j<n;j++){ 
				  
					  wind = (i+f)*n2+j+f;
				
					  input2[wind] = input[i*n+j]; 
				  }	  
				  for(j=0;j<f;j++){ 
					   iw1 = ((i+f)*n2)+j;
					   iw2 = ((i+f)*n2)+j+n+f;
					   input2[iw1] = input[(i*n)+f-j-1]; 
					   input2[iw2] = input[(i*n)+n-j-1];
				  }
			  }
			     
		 
		   for(i=0;i<f;i++){ 
		   for(j=0;j<n2;j++){ 
					  iw1 = i*n2+j; 
					  iw2 = (i+m+f)*n2+j; 
					  input2[iw1] = input2[(f+f-i-1)*n2+j]; 
					  input2[iw2] = input2[(f+m-i-1)*n2+j]; 
					 
				  } 
			  } 
		  

	 return; 
	 
 }


void mexFunction(
	int nlhs,
	mxArray *plhs[],
	int nrhs,
	const mxArray *prhs[] )

{ 

	double *input, *outp;
	double h; 
	int t,f; 
	int m, n;
	int ndm[2];

	if(nrhs ==0) 
	mexErrMsgTxt("Usage: mexnlm2D(input,t,f,h,m,n");
    if(nrhs !=6 && nrhs !=0)
	mexErrMsgTxt("Six inputs Required"); 
	if( nlhs > 1) 
	mexErrMsgTxt("Too many output arguments."); 
	

	/*get the length of the input vectors */ 

    input  = mxGetPr(prhs[0]); 
	t  = (int) mxGetScalar(prhs[1]); 
	f  = (int) mxGetScalar(prhs[2]);
	h  =   mxGetScalar(prhs[3]); 
	m  = (int) mxGetScalar(prhs[4]); 
	n  = (int) mxGetScalar(prhs[5]);
	 
    outp = mxCalloc(m*n,sizeof(double));

         ndm[0] =  m;
	ndm[1] =  n; 

	/* create new array set output pointer to */ 
   
	plhs[0] = mxCreateNumericArray(2,ndm,mxDOUBLE_CLASS,mxREAL);
	outp  = mxGetPr(plhs[0]); 
	/* call the C subroutine */
	/* watch out: permuted dimensions ! */ 
	mexnlm2D(outp,input,t,f,h,n,m);
	return; 
}
