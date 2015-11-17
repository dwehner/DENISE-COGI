/*------------------------------------------------------------------------
 * Pure forward modelling code
 *
 * D. Koehn
 * Kiel, 3.2.2012
 *  ----------------------------------------------------------------------*/

#include "fd.h"


float step_length_est1(FILE *fprec, float ** waveconv, float ** waveconv_rho, float ** waveconv_u, float ** prho, float ** prhonp1, float ** ppi, float ** ppinp1, int iter, int nfstart,
	int nsrc, float ** puipjp, float ** prip, float ** prjp, float L2, int partest, float ** srcpos_loc, float ** srcpos, float ** srcpos1, float ** signals, int ns,
	int nd, float ** pvx, float ** pvy, float ** psxx, float ** psyy, float ** psxy, float ** ux, float ** uy, float ** pvxp1, float ** pvyp1, float ** psi_sxx_x, float ** psi_sxy_x,
	float ** psi_vxx, float ** psi_vyx, float ** psi_syy_y, float ** psi_sxy_y, float ** psi_vyy, float ** psi_vxy, float ** psi_vxxs, float ** pvxm1, float ** pvym1, float ** uttx, 
	float ** utty, float ** absorb_coeff, float *hc, float * K_x, float * a_x, float * b_x, float * K_x_half, float * a_x_half, float * b_x_half, float * K_y, float * a_y, float * b_y,  
	float * K_y_half, float * a_y_half, float * b_y_half, float ** uxy, float ** uyx, int ntr, int **recpos_loc, float **sectionvx, float **sectionvy, float **sectionp, float **sectioncurl, 
	float **sectiondiv, float **sectionread, int ntr_glob, float ** sectionvxdata, float ** sectionvxdiff, float ** sectionvxdiffold, float ** sectionvydata, float ** sectionvydiff, 
	float ** sectionvydiffold, float * epst1, float * L2t, float L2sum, float energy_sum, float ** bufferlef_to_rig, float ** bufferrig_to_lef, 
	float ** buffertop_to_bot, float ** bufferbot_to_top, float **pu, float **punp1, int itest, int nsrc_glob, int nsrc_loc, MPI_Request * req_send, MPI_Request * req_rec){

	extern int MYID;
	
	float opteps_vp, eps_scale;

	int h, i, j, n, nshots, ishot, nt, lsnap, lsamp, nsnap, infoout;

	/* Variables for step length calculation */
        int step1, step2, itests, iteste, stepmax, countstep;
        float scalefac, al, ar, ac, fl, fr, fc, a0, an, fn, a2, f2;
	float acr, bcr, arl, brl, alc, blc, num, den;

	/* Parameters for step length calculations */
        stepmax = 20; /* number of maximum misfit calculations/steplength 2/3*/
	itest=2;

/* ============================================================================================================================*/
/* =============================================== test loop L2 ===============================================================*/
/* ============================================================================================================================*/

/* start with first guess for step length alpha */
eps_scale=0.01; /* maximum model change = 1% of the maximum model value */
countstep=0;     /* count number of forward calculations */
a0=eps_scale;

scalefac = 2.0; /* Interval bound expansion rate */

/* Left edge */
al = 0;  

/* calculate L2 norm */
eps_scale=al; 
forward_mod(fprec,waveconv,waveconv_rho,waveconv_u,prho,prhonp1,ppi,ppinp1,iter,eps_scale,nfstart,nsrc,puipjp,prip,prjp,L2,partest,srcpos_loc,srcpos,srcpos1,signals,ns,
	    nd,pvx,pvy,psxx,psyy,psxy,ux,uy,pvxp1,pvyp1,psi_sxx_x,psi_sxy_x,psi_vxx,psi_vyx,psi_syy_y,psi_sxy_y,psi_vyy,psi_vxy,psi_vxxs,pvxm1,pvym1,uttx,utty,absorb_coeff,hc,K_x, 
	    a_x,b_x,K_x_half,a_x_half,b_x_half,K_y,a_y,b_y,K_y_half,a_y_half,b_y_half,uxy,uyx,ntr,recpos_loc,sectionvx,sectionvy,sectionp,sectioncurl,sectiondiv,sectionread,ntr_glob,
	    sectionvxdata,sectionvxdiff,sectionvxdiffold,sectionvydata,sectionvydiff,sectionvydiffold,epst1,L2t,L2sum,energy_sum,bufferlef_to_rig,bufferrig_to_lef, 
	    buffertop_to_bot,bufferbot_to_top,pu,punp1,itest,nsrc_glob,nsrc_loc,req_send,req_rec);
fl = L2t[2];

/* Right edge */
ar = fabs(a0);      

/* calculate L2 norm */
eps_scale=ar; 
forward_mod(fprec,waveconv,waveconv_rho,waveconv_u,prho,prhonp1,ppi,ppinp1,iter,eps_scale,nfstart,nsrc,puipjp,prip,prjp,L2,partest,srcpos_loc,srcpos,srcpos1,signals,ns,
	    nd,pvx,pvy,psxx,psyy,psxy,ux,uy,pvxp1,pvyp1,psi_sxx_x,psi_sxy_x,psi_vxx,psi_vyx,psi_syy_y,psi_sxy_y,psi_vyy,psi_vxy,psi_vxxs,pvxm1,pvym1,uttx,utty,absorb_coeff,hc,K_x, 
	    a_x,b_x,K_x_half,a_x_half,b_x_half,K_y,a_y,b_y,K_y_half,a_y_half,b_y_half,uxy,uyx,ntr,recpos_loc,sectionvx,sectionvy,sectionp,sectioncurl,sectiondiv,sectionread,ntr_glob,
	    sectionvxdata,sectionvxdiff,sectionvxdiffold,sectionvydata,sectionvydiff,sectionvydiffold,epst1,L2t,L2sum,energy_sum,bufferlef_to_rig,bufferrig_to_lef, 
	    buffertop_to_bot,bufferbot_to_top,pu,punp1,itest,nsrc_glob,nsrc_loc,req_send,req_rec);
fr = L2t[2];


/* Center */
ac = (ar + al)/2;  

/* calculate L2 norm */
eps_scale=ac; 
forward_mod(fprec,waveconv,waveconv_rho,waveconv_u,prho,prhonp1,ppi,ppinp1,iter,eps_scale,nfstart,nsrc,puipjp,prip,prjp,L2,partest,srcpos_loc,srcpos,srcpos1,signals,ns,
	    nd,pvx,pvy,psxx,psyy,psxy,ux,uy,pvxp1,pvyp1,psi_sxx_x,psi_sxy_x,psi_vxx,psi_vyx,psi_syy_y,psi_sxy_y,psi_vyy,psi_vxy,psi_vxxs,pvxm1,pvym1,uttx,utty,absorb_coeff,hc,K_x, 
	    a_x,b_x,K_x_half,a_x_half,b_x_half,K_y,a_y,b_y,K_y_half,a_y_half,b_y_half,uxy,uyx,ntr,recpos_loc,sectionvx,sectionvy,sectionp,sectioncurl,sectiondiv,sectionread,ntr_glob,
	    sectionvxdata,sectionvxdiff,sectionvxdiffold,sectionvydata,sectionvydiff,sectionvydiffold,epst1,L2t,L2sum,energy_sum,bufferlef_to_rig,bufferrig_to_lef, 
	    buffertop_to_bot,bufferbot_to_top,pu,punp1,itest,nsrc_glob,nsrc_loc,req_send,req_rec);
fc = L2t[2];


countstep = 2;

while (countstep<stepmax){
    if (al>=0){
        if (fc>=fr){
            al = ac;   
			fl = fc;
            ac = ar;   
			fc = fr;
            ar = ar*scalefac; 
			/*fr = L2func(x+dx*ar,y+dy*ar);*/

	    /* calculate L2 norm */
        eps_scale=ar; 
        forward_mod(fprec,waveconv,waveconv_rho,waveconv_u,prho,prhonp1,ppi,ppinp1,iter,eps_scale,nfstart,nsrc,puipjp,prip,prjp,L2,partest,srcpos_loc,srcpos,srcpos1,signals,ns,
	    nd,pvx,pvy,psxx,psyy,psxy,ux,uy,pvxp1,pvyp1,psi_sxx_x,psi_sxy_x,psi_vxx,psi_vyx,psi_syy_y,psi_sxy_y,psi_vyy,psi_vxy,psi_vxxs,pvxm1,pvym1,uttx,utty,absorb_coeff,hc,K_x, 
	    a_x,b_x,K_x_half,a_x_half,b_x_half,K_y,a_y,b_y,K_y_half,a_y_half,b_y_half,uxy,uyx,ntr,recpos_loc,sectionvx,sectionvy,sectionp,sectioncurl,sectiondiv,sectionread,ntr_glob,
	    sectionvxdata,sectionvxdiff,sectionvxdiffold,sectionvydata,sectionvydiff,sectionvydiffold,epst1,L2t,L2sum,energy_sum,bufferlef_to_rig,bufferrig_to_lef, 
	    buffertop_to_bot,bufferbot_to_top,pu,punp1,itest,nsrc_glob,nsrc_loc,req_send,req_rec);
        fr = L2t[2];
        
		countstep    = countstep + 1;
		
		}else if (fc>=fl){
            ar = ac;    
			fr = fc;
            ac = al;    
			fc = fl;
            al = al/10; 
			/*fl = L2func(x+dx*ar,y+dy*ar);*/

		/* calculate L2 norm */
        eps_scale=al; 
        forward_mod(fprec,waveconv,waveconv_rho,waveconv_u,prho,prhonp1,ppi,ppinp1,iter,eps_scale,nfstart,nsrc,puipjp,prip,prjp,L2,partest,srcpos_loc,srcpos,srcpos1,signals,ns,
	    nd,pvx,pvy,psxx,psyy,psxy,ux,uy,pvxp1,pvyp1,psi_sxx_x,psi_sxy_x,psi_vxx,psi_vyx,psi_syy_y,psi_sxy_y,psi_vyy,psi_vxy,psi_vxxs,pvxm1,pvym1,uttx,utty,absorb_coeff,hc,K_x, 
	    a_x,b_x,K_x_half,a_x_half,b_x_half,K_y,a_y,b_y,K_y_half,a_y_half,b_y_half,uxy,uyx,ntr,recpos_loc,sectionvx,sectionvy,sectionp,sectioncurl,sectiondiv,sectionread,ntr_glob,
	    sectionvxdata,sectionvxdiff,sectionvxdiffold,sectionvydata,sectionvydiff,sectionvydiffold,epst1,L2t,L2sum,energy_sum,bufferlef_to_rig,bufferrig_to_lef, 
	    buffertop_to_bot,bufferbot_to_top,pu,punp1,itest,nsrc_glob,nsrc_loc,req_send,req_rec);
        fl = L2t[2];

            countstep    = countstep + 1;        
		}else{  /* fc<fr & fc<fl,*/
			break;}
	}
        
    else{
        if (fc>=fr){
            ac = ar;   
			fc = fr;
            ar = ar*scalefac; 
			/*fr = L2func(x+dx*ar,y+dy*ar);*/

		/* calculate L2 norm */
        eps_scale=ar; 
        forward_mod(fprec,waveconv,waveconv_rho,waveconv_u,prho,prhonp1,ppi,ppinp1,iter,eps_scale,nfstart,nsrc,puipjp,prip,prjp,L2,partest,srcpos_loc,srcpos,srcpos1,signals,ns,
	    nd,pvx,pvy,psxx,psyy,psxy,ux,uy,pvxp1,pvyp1,psi_sxx_x,psi_sxy_x,psi_vxx,psi_vyx,psi_syy_y,psi_sxy_y,psi_vyy,psi_vxy,psi_vxxs,pvxm1,pvym1,uttx,utty,absorb_coeff,hc,K_x, 
	    a_x,b_x,K_x_half,a_x_half,b_x_half,K_y,a_y,b_y,K_y_half,a_y_half,b_y_half,uxy,uyx,ntr,recpos_loc,sectionvx,sectionvy,sectionp,sectioncurl,sectiondiv,sectionread,ntr_glob,
	    sectionvxdata,sectionvxdiff,sectionvxdiffold,sectionvydata,sectionvydiff,sectionvydiffold,epst1,L2t,L2sum,energy_sum,bufferlef_to_rig,bufferrig_to_lef, 
	    buffertop_to_bot,bufferbot_to_top,pu,punp1,itest,nsrc_glob,nsrc_loc,req_send,req_rec);
        fr = L2t[2];

            countstep    = countstep + 1; 
		}else if (fc>=fl){
            ac = al;   
		    fc = fr;
            al = al*scalefac; 
			/*fl = L2func(x+dx*al,y+dy*al);*/

		/* calculate L2 norm */
        eps_scale=al; 
        forward_mod(fprec,waveconv,waveconv_rho,waveconv_u,prho,prhonp1,ppi,ppinp1,iter,eps_scale,nfstart,nsrc,puipjp,prip,prjp,L2,partest,srcpos_loc,srcpos,srcpos1,signals,ns,
	    nd,pvx,pvy,psxx,psyy,psxy,ux,uy,pvxp1,pvyp1,psi_sxx_x,psi_sxy_x,psi_vxx,psi_vyx,psi_syy_y,psi_sxy_y,psi_vyy,psi_vxy,psi_vxxs,pvxm1,pvym1,uttx,utty,absorb_coeff,hc,K_x, 
	    a_x,b_x,K_x_half,a_x_half,b_x_half,K_y,a_y,b_y,K_y_half,a_y_half,b_y_half,uxy,uyx,ntr,recpos_loc,sectionvx,sectionvy,sectionp,sectioncurl,sectiondiv,sectionread,ntr_glob,
	    sectionvxdata,sectionvxdiff,sectionvxdiffold,sectionvydata,sectionvydiff,sectionvydiffold,epst1,L2t,L2sum,energy_sum,bufferlef_to_rig,bufferrig_to_lef, 
	    buffertop_to_bot,bufferbot_to_top,pu,punp1,itest,nsrc_glob,nsrc_loc,req_send,req_rec);
        fl = L2t[2];


            countstep    = countstep + 1; 
		}else{
            break;
        }
    }    
}
       
while (countstep<stepmax){
    if ((fc<fl)&&(fc<fr)){ /* Criteria okay - do quadratic search iteration*/
        acr = ac - ar; 
	    bcr = ac*ac - (ar*ar);
        arl = ar - al; 
		brl = ar*ar - (al*al);
        alc = al - ac; 
		blc = al*al - (ac*ac);
        
        num = (bcr*fl + brl*fc + blc*fr);
        den = (acr*fl + arl*fc + alc*fr);

        if(den==0){ /* Optimal solution found*/
            break;
		}

        an = 0.5*num/den;
        /*fn = L2func(x+dx*an,y+dy*an);*/

		/* calculate L2 norm */
        eps_scale=an; 
        forward_mod(fprec,waveconv,waveconv_rho,waveconv_u,prho,prhonp1,ppi,ppinp1,iter,eps_scale,nfstart,nsrc,puipjp,prip,prjp,L2,partest,srcpos_loc,srcpos,srcpos1,signals,ns,
	    nd,pvx,pvy,psxx,psyy,psxy,ux,uy,pvxp1,pvyp1,psi_sxx_x,psi_sxy_x,psi_vxx,psi_vyx,psi_syy_y,psi_sxy_y,psi_vyy,psi_vxy,psi_vxxs,pvxm1,pvym1,uttx,utty,absorb_coeff,hc,K_x, 
	    a_x,b_x,K_x_half,a_x_half,b_x_half,K_y,a_y,b_y,K_y_half,a_y_half,b_y_half,uxy,uyx,ntr,recpos_loc,sectionvx,sectionvy,sectionp,sectioncurl,sectiondiv,sectionread,ntr_glob,
	    sectionvxdata,sectionvxdiff,sectionvxdiffold,sectionvydata,sectionvydiff,sectionvydiffold,epst1,L2t,L2sum,energy_sum,bufferlef_to_rig,bufferrig_to_lef, 
	    buffertop_to_bot,bufferbot_to_top,pu,punp1,itest,nsrc_glob,nsrc_loc,req_send,req_rec);
        fn = L2t[2];

		countstep = countstep + 1;
        
	    if (an>ac){
            if(fn>=fc){
                ar = an; 
				fr = fn;
			}else{
                al = ac; 
				fl = fc;
                ac = an; 
				fc = fn;
			}
		}else{
            if(fn>=fc){
                al = an; 
				fl = fn;
			}else{
                ar = ac; 
				fr = fc;
                ac = an; 
				fc = fn;
			}   
		}

	}else { /* Safety switchover to bisection - safeguard technique */
        a2  = ac + 0.001*(ar-al);
        /*f2  = L2func(x+dx*a2,y+dy*a2);*/

		/* calculate L2 norm */
        eps_scale=a2; 
        forward_mod(fprec,waveconv,waveconv_rho,waveconv_u,prho,prhonp1,ppi,ppinp1,iter,eps_scale,nfstart,nsrc,puipjp,prip,prjp,L2,partest,srcpos_loc,srcpos,srcpos1,signals,ns,
	    nd,pvx,pvy,psxx,psyy,psxy,ux,uy,pvxp1,pvyp1,psi_sxx_x,psi_sxy_x,psi_vxx,psi_vyx,psi_syy_y,psi_sxy_y,psi_vyy,psi_vxy,psi_vxxs,pvxm1,pvym1,uttx,utty,absorb_coeff,hc,K_x, 
	    a_x,b_x,K_x_half,a_x_half,b_x_half,K_y,a_y,b_y,K_y_half,a_y_half,b_y_half,uxy,uyx,ntr,recpos_loc,sectionvx,sectionvy,sectionp,sectioncurl,sectiondiv,sectionread,ntr_glob,
	    sectionvxdata,sectionvxdiff,sectionvxdiffold,sectionvydata,sectionvydiff,sectionvydiffold,epst1,L2t,L2sum,energy_sum,bufferlef_to_rig,bufferrig_to_lef, 
	    buffertop_to_bot,bufferbot_to_top,pu,punp1,itest,nsrc_glob,nsrc_loc,req_send,req_rec);
        f2 = L2t[2];

        countstep    = countstep + 1;

        if(fc<f2){
            ar = a2; 
			fr = f2;
		}else{ /* fc>fr*/
            al = ac; 
			fl = fc;
		}

        ac     = (al+ar)/2; 
        /*fc     = L2func(x+dx*ac,y+dy*ac);*/

		/* calculate L2 norm */
        eps_scale=ac; 
        forward_mod(fprec,waveconv,waveconv_rho,waveconv_u,prho,prhonp1,ppi,ppinp1,iter,eps_scale,nfstart,nsrc,puipjp,prip,prjp,L2,partest,srcpos_loc,srcpos,srcpos1,signals,ns,
	    nd,pvx,pvy,psxx,psyy,psxy,ux,uy,pvxp1,pvyp1,psi_sxx_x,psi_sxy_x,psi_vxx,psi_vyx,psi_syy_y,psi_sxy_y,psi_vyy,psi_vxy,psi_vxxs,pvxm1,pvym1,uttx,utty,absorb_coeff,hc,K_x, 
	    a_x,b_x,K_x_half,a_x_half,b_x_half,K_y,a_y,b_y,K_y_half,a_y_half,b_y_half,uxy,uyx,ntr,recpos_loc,sectionvx,sectionvy,sectionp,sectioncurl,sectiondiv,sectionread,ntr_glob,
	    sectionvxdata,sectionvxdiff,sectionvxdiffold,sectionvydata,sectionvydiff,sectionvydiffold,epst1,L2t,L2sum,energy_sum,bufferlef_to_rig,bufferrig_to_lef, 
	    buffertop_to_bot,bufferbot_to_top,pu,punp1,itest,nsrc_glob,nsrc_loc,req_send,req_rec);
        fc = L2t[2];

        countstep    = countstep + 1;        
	}

    if((ar-al==0)||((fl==fr)&&(fl==fc))){
        break;
	}
}

/*if(MYID==0){
printf("MYID = %d \t L2t[1] = %e \t L2t[2] = %e \t L2t[3] = %e \n",MYID,L2t[1],L2t[2],L2t[3]);
printf("MYID = %d \t epst1[1] = %e \t epst1[2] = %e \t epst1[3] = %e \n",MYID,epst1[1],epst1[2],epst1[3]);}*/

L2t[1] = fl;
L2t[2] = fc;
L2t[3] = fr;

epst1[1] = al;
epst1[2] = ac;
epst1[3] = ar;

if(countstep>=stepmax){
  err(" Steplength estimation failed!");
}

if(ac==0.0){
  err(" Steplength estimation failed!");
}
  
return ac;
}

