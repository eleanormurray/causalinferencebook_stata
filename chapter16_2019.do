/***************************************************************
Stata code for Causal Inference: What If by Miguel Hernan & Jamie Robins
Date: 10/10/2019
Author: Eleanor Murray 
For errors contact: ejmurray@bu.edu
***************************************************************/

/***************************************************************
PROGRAM 16.1
Estimating the average causal effect using the standard IV estimator
via the calculation of sample averages
Data from NHEFS
Section 16.2
***************************************************************/

clear
use "nhefs.dta"

summarize price82

/* ignore subjects with missing outcome or missing instrument for simplicity*/
foreach var of varlist wt82 price82 {
drop if `var'==.
}


/*Create categorical instrument*/
gen byte highprice  = (price82 > 1.5 & price82 < .)

/*Calculate P[A|Z=Z]*/
tab highprice qsmk, row

/*Calculate P[Y|Z=z]*/
ttest wt82_71, by(highprice)

/*Final IV estimate, OPTION 1: Hand calculations*/
/*Numerator: num = E[Y|Z=1] - E[Y|Z=0] = 2.686 - 2.536 = 0.150*/
/*Denominator: denom = P[A=1|Z=1] - P[A=1|Z=0] = 0.258 - 0.195 = 0.063 */ 
/*IV estimator: E[Ya=1] - E[Ya=0] = (E[Y|Z=1]-E[Y|Z=0])/(P[A=1|Z=1]-P[A=1|Z=0]) = 0.150/0.063 = 2.397*/
display 2.686 - 2.536
display 0.258 - 0.195 
display 0.150/0.063

/*OPTION 2 2: automated calculation of instrument*/
*Calculate P[A=1|Z=z], for each value of the instrument, and store in a matrix*
quietly summarize qsmk if (highprice==0)
matrix input pa = (`r(mean)')
quietly summarize qsmk if (highprice==1)
matrix pa = (pa ,`r(mean)')
matrix list pa
*Calculate P[Y|Z=z], for each value of the instrument, and store in a second matrix*
quietly summarize wt82_71 if (highprice==0)
matrix input ey = (`r(mean)')
quietly summarize wt82_71 if (highprice==1)
matrix ey = (ey ,`r(mean)')
matrix list ey
*Using Stata's built-in matrix manipulation feature (Mata), calculate numerator, denominator and IV estimator*
*Numerator: num = E[Y|Z=1] - E[Y|Z=0]*mata
*Denominator: denom = P[A=1|Z=1] - P[A=1|Z=0]*
*IV estimator: iv_est = IV estimate of E[Ya=1] - E[Ya=0] *
mata
pa = st_matrix("pa")
ey = st_matrix("ey")
num = ey[1,2] - ey[1,1] 
denom = pa[1,2] - pa[1,1]
iv_est = num / denom 
num
denom
iv_est
end


/***************************************************************
PROGRAM 16.2
Estimating the average causal effect using the standard IV estimator
via two-stage-least-squares regression
Data from NHEFS
Section 16.2
***************************************************************/

/*ivregress fits the model in two stages: */
/*first model: qsmk = highprice*/
/*second model: wt82_71 = predicted_qsmk*/
ivregress 2sls wt82_71 (qsmk = highprice)

/***************************************************************************
PROGRAM 16.3
Estimating the average causal effect using the standard IV estimator 
via an additive marginal structural model
Data from NHEFS
Checking one possible value of psi.
See Chapter 14 for program that checks several values 
and computes 95% confidence intervals   
Section 16.2
***************************************************************************/

gen psi = 2.396
gen hspi = wt82_71 -psi*qsmk

logit highprice hspi


/***************************************************************************
PROGRAM 16.4
Estimating the average causal effect using the standard IV estimator
based on alternative proposed instruments
Data from NHEFS
Section 16.5
***************************************************************************/

/*Instrument cut-point: 1.6*/
replace highprice = .
replace highprice = (price82 >1.6 & price82 < .)

ivregress 2sls wt82_71 (qsmk = highprice)


/*Instrument cut-point: 1.7*/
replace highprice = .
replace highprice = (price82 >1.7 & price82 < .)

ivregress 2sls wt82_71 (qsmk = highprice)


/*Instrument cut-point: 1.8*/
replace highprice = .
replace highprice = (price82 >1.8 & price82 < .)

ivregress 2sls wt82_71 (qsmk = highprice)


/*Instrument cut-point: 1.9*/
replace highprice = .
replace highprice = (price82 >1.9 & price82 < .)

ivregress 2sls wt82_71 (qsmk = highprice)



/***************************************************************************
PROGRAM 16.5
Estimating the average causal effect using the standard IV estimator
conditional on baseline covariates
Data from NHEFS
Section 16.5
***************************************************************************/

replace highprice = .
replace highprice = (price82 >1.5 & price82 < .)

ivregress 2sls wt82_71 sex race c.age c.smokeintensity c.smokeyrs i.exercise i.active c.wt7 (qsmk = highprice)
