options FULLSTIMER THREADS /*CPUCOUNT=ACTUAL*/;

* HERE WE RUN THE PROPOSED MODEL ON THE FIRST SIMULATED DATA SET FROM DGP 1; 

* THE DATA WAS GENERATED IN R - IMPORT INTO SAS USING PROC IML (MAKE SURE YOUR SAS VERSION IS SET-UP TO DO THIS) ;
options set=R_HOME="C:\Users\lovblom\Documents\R\R-3.6.3";
proc iml;
	submit / R;

	setwd('M:/lovblom/PHD/DATASETS/Erik Data Processing/GITHUB/CHAPTER4')
	sasdata <- readRDS('sasdata_DGP1_n1.rds')

	endsubmit;
	call ImportDataSetFromR("sasdata", "sasdata");
quit;

* NEXT, SOME DATA-WRANGLING STEPS;
proc sort data=sasdata;
	by sim_number id outcome_number year;
run;

data sasdata2;
	set sasdata;
	by sim_number id outcome_number;

	if outcome_number = 1 then do;
		lag_year=lag(year);
		lag_state=lag(state);

		if last.outcome_number then last=1;
		else last=0;
	end;
run;

data sasdata3;
	set sasdata2;
		if outcome_number = 1 then do;
		if last=0 then do;
			deltaI=1*(state=1 and lag_state=1)+2*(state=2 and lag_state=1)+3*(state=2 and lag_state=2);
		end;
		else if last=1 then do;
			deltaC = 1*(state=1 and lag_state=1) + 2*(state=2 and lag_state=1) + 3*(state=2 and lag_state=2) + 
				4*(state=999 and lag_state=1) + 5*(state=999 and lag_state=2) + 
				6*(state=3 and lag_state=1) + 7*(state=3 and lag_state=2);
		end;
	end;
run;

data sasdata4;
	set sasdata3;
	if outcome_number=1 and year=0 then delete;*THE FIRST MULTISTATE VISIT MUST BE DELETED;
run;

*LASTLY, RUN THE MODEL USING PROC NLMIXED;
ods output ParameterEstimates=oute
			ConvergenceStatus=mix
			FitStatistics=lrt;
proc nlmixed data=sasdata4 itdetails method=gauss tech=newrap gconv=0 EBOPT
	NTHREADS=6/*NTHREADS=4*/ /*ADJUST THIS OPTION AS NEEDED, OR REMOVE IT*/
;
	*by sim_number;

	title 'PROPOSED JOINT MODEL, DGP 1 SIMULATION #1';

	parms 
	lambda121=-3.5 lambda122=-5 lambda123=-3.9 
	lambda231=-5.7 lambda232=-5.9 lambda233=-6.5

	gamma121 = 0.08 gamma122 = -0.32
	gamma231 = 0.13 /*gamma232 = 0*/

	alpha12=0.12
	alpha23=0.06

	beta0=0.39 beta1=0.19 beta2=0.21
	d11=1.81659 d21=-0.05 d22=0.2236
	sigma=1.0977
	;

	if outcome_number=1 then do;

		pwc_knot1 = 4.8;
		pwc_knot2 = 17.4182;

		array xx[15] x1-x15;
		array wt[15] wt1-wt15;
		x1=0.991455371120813	;
		x2=0.949107912342759	;
		x3=0.864864423359769	;
		x4=0.741531185599394	;
		x5=0.586087235467691	;
		x6=0.405845151377397	;
		x7=0.207784955007898	;
		x8=0;
		x9=-1*x1; x10=-1*x2; x11=-1*x3; x12=-1*x4; x13=-1*x5; x14=-1*x6; x15=-1*x7;
		wt1=0.022935322010529	;
		wt2=0.063092092629979	;
		wt3=0.104790010322250	;
		wt4=0.140653259715525	;
		wt5=0.169004726639267	;
		wt6=0.190350578064785	;
		wt7=0.204432940075298	;
		wt8=0.209482141084728	;
		wt9=wt1; wt10=wt2; wt11=wt3; wt12=wt4; wt13=wt5; wt14=wt6; wt15=wt7;

		if deltaI=1 or deltaC=1 then do;*p11 only required;
			gauss_total = 0;
			do i=1 to 15;
				x_scaled = ((xx[i]+1)/2)*(year-lag_year)+lag_year;
				wt_scaled = ((year-lag_year)/2)*wt[i];
				H1 = wt_scaled*exp( lambda121*(x_scaled<=pwc_knot1) + lambda122*(pwc_knot1<x_scaled<=pwc_knot2) + lambda123*(pwc_knot2<x_scaled) + gamma121*X1cont + gamma122*X3cat
									+ alpha12*(beta0 + beta1*x_scaled + b0 + b1*x_scaled + beta2*X1cont) );
				gauss_total + H1;
			end;
			*p11 = exp(-gauss_total);
			log_p11 = -gauss_total;
			log_lik = log_p11;
		end;

		else if deltaI=2 or deltaC=2 or deltaC=6 then do;*p12 only required;
			gauss_total=0;
			do i=1 to 15;
				gaussH1=0;
				gaussH2=0;
				x_scaled = ((xx[i]+1)/2)*(year-lag_year)+lag_year;
				wt_scaled = ((year-lag_year)/2)*wt[i];
				q12 = exp(lambda121*(x_scaled<=pwc_knot1) + lambda122*(pwc_knot1<x_scaled<=pwc_knot2) + lambda123*(pwc_knot2<x_scaled) + gamma121*X1cont + gamma122*X3cat
									+ alpha12*(beta0 + beta1*x_scaled + b0 + b1*x_scaled + beta2*X1cont) );	
				do j=1 to 15;
					x2_scaled = ((xx[j]+1)/2)*(x_scaled - lag_year) + lag_year;
					wt2_scaled = ((x_scaled - lag_year)/2)*wt[j];
					H1 = wt2_scaled*exp( lambda121*(x2_scaled<=pwc_knot1) + lambda122*(pwc_knot1<x2_scaled<=pwc_knot2) + lambda123*(pwc_knot2<x2_scaled) + gamma121*X1cont + gamma122*X3cat
										+ alpha12*(beta0 + beta1*x2_scaled + b0 + b1*x2_scaled + beta2*X1cont) );
					x3_scaled = ((xx[j]+1)/2)*(year - x_scaled) + x_scaled;
					wt3_scaled = ((year - x_scaled)/2)*wt[j];
					H2 = wt3_scaled*exp( lambda231*(x3_scaled<=pwc_knot1) + lambda232*(pwc_knot1<x3_scaled<=pwc_knot2) + lambda233*(pwc_knot2<x3_scaled) + gamma231*X1cont /*+ gamma232*X3cat*/
										+ alpha23*(beta0 + beta1*x3_scaled + b0 + b1*x3_scaled + beta2*X1cont));
					gaussH1+H1;
					gaussH2+H2;
					piece = exp(-gaussH1)*q12*exp(-gaussH2);
				end;
				weighted_piece = wt_scaled*piece;
				gauss_total + weighted_piece;
			end;
			integral_approximation = gauss_total;
			p12 = integral_approximation;
			log_p12 = log(p12);

			if deltaI=2 or deltaC=2 then do;
				log_lik = log_p12;
			end;
			else if deltaC=6 then do;
				log_lik = log_p12 + lambda231*(year<=pwc_knot1) + lambda232*(pwc_knot1<year<=pwc_knot2) + lambda233*(pwc_knot2<year) + gamma231*X1cont /*+ gamma232*X3cat*/
								+ alpha23*(beta0 + beta1*year + b0 + b1*year + beta2*X1cont);
			end;
		end;

		else if deltaI=3 or deltaC=3 or deltaC=5 or deltaC=7 then do;*p22 only required;
			gauss_total = 0;
			do i=1 to 15;
				x_scaled = ((xx[i]+1)/2)*(year-lag_year)+lag_year;
				wt_scaled = ((year-lag_year)/2)*wt[i];
				H2 = wt_scaled*exp( lambda231*(x_scaled<=pwc_knot1) + lambda232*(pwc_knot1<x_scaled<=pwc_knot2) + lambda233*(pwc_knot2<x_scaled)+ gamma231*X1cont /*+ gamma232*X3cat*/
									+ alpha23*(beta0 + beta1*x_scaled + b0 + b1*x_scaled + beta2*X1cont) );
				gauss_total + H2;
			end;
			*p22 = exp(-gauss_total);
			log_p22 = -gauss_total;

			if deltaI=3 or deltaC=3 or deltaC=5 then do;
				log_lik = log_p22;
			end;
			else if deltaC=7 then do;
				log_lik = log_p22 + lambda231*(year<=pwc_knot1) + lambda232*(pwc_knot1<year<=pwc_knot2) + lambda233*(pwc_knot2<year) + gamma231*X1cont /*+ gamma232*X3cat*/
							+ alpha23*(beta0 + beta1*year + b0 + b1*year + beta2*X1cont);
			end;
		end;

		else if deltaC=4 then do;*p11 and p12 required;
			gauss_total_p11 = 0;
			do i=1 to 15;
				x_scaled = ((xx[i]+1)/2)*(year-lag_year)+lag_year;
				wt_scaled = ((year-lag_year)/2)*wt[i];
				H1 = wt_scaled*exp( lambda121*(x_scaled<=pwc_knot1) + lambda122*(pwc_knot1<x_scaled<=pwc_knot2) + lambda123*(pwc_knot2<x_scaled) + gamma121*X1cont + gamma122*X3cat
									+ alpha12*(beta0 + beta1*x_scaled + b0 + b1*x_scaled + beta2*X1cont) );
				gauss_total_p11 + H1;
			end;
			p11 = exp(-gauss_total_p11);

			gauss_total_p12=0;
			do i=1 to 15;
				gaussH1=0;
				gaussH2=0;
				x_scaled = ((xx[i]+1)/2)*(year-lag_year)+lag_year;
				wt_scaled = ((year-lag_year)/2)*wt[i];
				q12 = exp(lambda121*(x_scaled<=pwc_knot1) + lambda122*(pwc_knot1<x_scaled<=pwc_knot2) + lambda123*(pwc_knot2<x_scaled) + gamma121*X1cont + gamma122*X3cat
									+ alpha12*(beta0 + beta1*x_scaled + b0 + b1*x_scaled + beta2*X1cont) );	
				do j=1 to 15;
					x2_scaled = ((xx[j]+1)/2)*(x_scaled - lag_year) + lag_year;
					wt2_scaled = ((x_scaled - lag_year)/2)*wt[j];
					H1 = wt2_scaled*exp( lambda121*(x2_scaled<=pwc_knot1) + lambda122*(pwc_knot1<x2_scaled<=pwc_knot2) + lambda123*(pwc_knot2<x2_scaled) + gamma121*X1cont + gamma122*X3cat
										+ alpha12*(beta0 + beta1*x2_scaled + b0 + b1*x2_scaled + beta2*X1cont) );
					x3_scaled = ((xx[j]+1)/2)*(year - x_scaled) + x_scaled;
					wt3_scaled = ((year - x_scaled)/2)*wt[j];
					H2 = wt3_scaled*exp( lambda231*(x3_scaled<=pwc_knot1) + lambda232*(pwc_knot1<x3_scaled<=pwc_knot2) + lambda233*(pwc_knot2<x3_scaled) + gamma231*X1cont /*+ gamma232*X3cat*/
										+ alpha23*(beta0 + beta1*x3_scaled + b0 + b1*x3_scaled + beta2*X1cont));
					gaussH1+H1;
					gaussH2+H2;
					piece = exp(-gaussH1)*q12*exp(-gaussH2);
				end;
				weighted_piece = wt_scaled*piece;
				gauss_total_p12 + weighted_piece;
			end;
			integral_approximation = gauss_total_p12;
			p12 = integral_approximation;
			log_lik = log(p11 + p12);
		end;

	end;

	else if outcome_number=2 then do;
		mean = beta0 + beta1*year + b0 + b1*year + beta2*X1cont;
		dens = -0.5*log(2*3.14159265358) - log(sigma) - 0.5*(y-mean)**2/(sigma**2);
		log_lik=dens;
	end;

	dummy=1;

	model dummy ~ general(log_lik);
	random b0 b1 ~ normal([0,0],[d11*d11,d21,d22*d22]) subject = id;

run;
