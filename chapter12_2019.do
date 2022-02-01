/***************************************************************
Stata code for Causal Inference: What If by Miguel Hernan & Jamie Robins
Date: 10/10/2019
Author: Eleanor Murray 
For errors contact: ejmurray@bu.edu
***************************************************************/

/***************************************************************
PROGRAM 12.1
Descriptive statistics from NHEFS data (Table 12.1)
***************************************************************/

clear
use "nhefs.dta"

/*Provisionally ignore subjects with missing values for follow-up weight*/
/*Sample size after exclusion: N = 1566*/
drop if wt82==.

/* Calculate mean weight change in those with and without smoking cessation*/
label define qsmk 0 "No smoking cessation" 1 "Smoking cessation"
label values qsmk qsmk
by qsmk, sort: egen years = mean(age) if age < . 
label var years "Age, years"
by qsmk, sort: egen male = mean(100 * (sex==0)) if sex < . 
label var male "Men, %"
by qsmk, sort: egen white = mean(100 * (race==0)) if race < . 
label var white "White, %"
by qsmk, sort: egen university = mean(100 * (education == 5)) if education < .
label var university "University, %"
by qsmk, sort: egen kg = mean(wt71) if wt71 < .
label var kg "Weight, kg"
by qsmk, sort: egen cigs = mean(smokeintensity) if smokeintensity < . 
label var cigs "Cigarettes/day"
by qsmk, sort: egen meansmkyrs = mean(smokeyrs) if smokeyrs < .
label var meansmkyrs "Years smoking"
by qsmk, sort: egen noexer = mean(100 * (exercise == 2)) if exercise < . 
label var noexer "Little/no exercise"
by qsmk, sort: egen inactive = mean(100 * (active==2)) if active < . 
label var inactive "Inactive daily life"

/*Output table*/
foreach var of varlist years male  white university kg cigs meansmkyrs noexer inactive {
tabdisp qsmk, cell(`var') format(%3.1f)
}



/***************************************************************
PROGRAM 12.2
Estimating IP weights for Section 12.2
Data from NHEFS
***************************************************************/

/*Fit a logistic model for the IP weights*/ 
logit qsmk sex race c.age##c.age ib(last).education c.smokeintensity##c.smokeintensity ///
c.smokeyrs##c.smokeyrs ib(last).exercise ib(last).active c.wt71##c.wt71 

/*Output predicted conditional probability of quitting smoking for each individual*/
predict p_qsmk, pr

/*Generate nonstabilized weights as P(A=1|covariates) if A = 1 and 1-P(A=1|covariates) if A = 0*/
gen w=.
replace w=1/p_qsmk if qsmk==1
replace w=1/(1-p_qsmk) if qsmk==0
/*Check the mean of the weights; we expect it to be close to 2.0*/
summarize w

/*Fit marginal structural model in the pseudopopulation*/
/*Weights assigned using pweight = w*/
/*Robust standard errors using cluster() option where 'seqn' is the ID variable*/
regress wt82_71 qsmk [pweight=w], cluster(seqn) 



/***************************************************************
PROGRAM 12.3
Estimating stabilized IP weights for Section 12.3
Data from NHEFS
***************************************************************/

/*Fit a logistic model for the denominator of the IP weights and predict the conditional probability of smoking*/ 
logit qsmk sex race c.age##c.age ib(last).education c.smokeintensity##c.smokeintensity ///
c.smokeyrs##c.smokeyrs ib(last).exercise ib(last).active c.wt71##c.wt71  
predict pd_qsmk, pr

/*Fit a logistic model for the numerator of ip weights and predict Pr(A=1) */ 
logit qsmk 
predict pn_qsmk, pr

/*Generate stabilized weights as f(A)/f(A|L)*/
gen sw_a=.
replace sw_a=pn_qsmk/pd_qsmk if qsmk==1
replace sw_a=(1-pn_qsmk)/(1-pd_qsmk) if qsmk==0

/*Check distribution of the stabilized weights*/
summarize sw_a

/*Fit marginal structural model in the pseudopopulation*/
regress wt82_71 qsmk [pweight=sw_a], cluster(seqn) 


/**********************************************************
FINE POINT 12.2
Checking positivity
**********************************************************/

/*Check for missing values within strata of covariates, for example: */
tab  age qsmk if race==0 & sex==1 & wt82!=.
tab  age qsmk if race==1 & sex==1 & wt82!=.



/***************************************************************
PROGRAM 12.4
Estimating the parameters of a marginal structural mean model
with a continuous treatment Data from NHEFS
Section 12.4
***************************************************************/
drop sw_a

/*Analysis restricted to subjects reporting <=25 cig/day at baseline: N = 1162*/
keep if  smokeintensity <=25

/*Fit a linear model for the denominator of the IP weights and calculate the mean expected smoking intensity*/ 
regress smkintensity82_71 sex race c.age##c.age ib(last).education c.smokeintensity##c.smokeintensity ///
c.smokeyrs##c.smokeyrs ib(last).exercise ib(last).active c.wt71##c.wt71
quietly predict p_den

/*Generate the denisty of the denomiator expectation using the mean expected smoking intensity and the residuals, assuming a normal distribution*/
/*Note: The regress command in STATA saves the root mean squared error for the immediate regression as e(rmse), thus there is no need to calculate it again. */
gen dens_den = normalden(smkintensity82_71, p_den, e(rmse))


/*Fit a linear model for the numerator of ip weights, calculate the mean expected value, and generate the density*/
quietly regress smkintensity82_71
quietly predict p_num
gen dens_num = normalden( smkintensity82_71, p_num, e(rmse))

/*Generate the final stabilized weights from the estimated numerator and denominator, and check the weights distribution*/
gen sw_a=dens_num/dens_den
summarize sw_a

/*Fit a marginal structural model in the pseudopopulation*/
regress wt82_71  c.smkintensity82_71##c.smkintensity82_71 [pweight=sw_a], cluster(seqn)


/*Output the estimated mean Y value when smoke intensity is unchanged from baseline to 1982 */
lincom _b[_cons]

/*Output the estimated mean Y value when smoke intensity increases by 20 from baseline to 1982*/
lincom _b[_cons] + 20*_b[smkintensity82_71 ] +400*_b[c.smkintensity82_71#c.smkintensity82_71]


/***************************************************************
PROGRAM 12.5
Estimating the parameters of a marginal structural logistic model
Data from NHEFS
Section 12.4
***************************************************************/

clear
use "nhefs.dta"

/*Provisionally ignore subjects with missing values for follow-up weight*/
/*Sample size after exclusion: N = 1566*/
drop if wt82==.

/*Estimate the stabilized weights for quitting smoking as in PROGRAM 12.3*/
/*Fit a logistic model for the denominator of the IP weights and predict the conditional probability of smoking*/ 
logit qsmk sex race c.age##c.age ib(last).education c.smokeintensity##c.smokeintensity ///
c.smokeyrs##c.smokeyrs ib(last).exercise ib(last).active c.wt71##c.wt71  
predict pd_qsmk, pr
/*Fit a logistic model for the numerator of ip weights and predict Pr(A=1) */ 
logit qsmk 
predict pn_qsmk, pr
/*Generate stabilized weights as f(A)/f(A|L)*/
gen sw_a=.
replace sw_a=pn_qsmk/pd_qsmk if qsmk==1
replace sw_a=(1-pn_qsmk)/(1-pd_qsmk) if qsmk==0
summarize sw_a

/*Fit marginal structural model in the pseudopopulation*/
/*NOTE: Stata has two commands for logistic regression, logit and logistic*/
/*Using logistic allows us to output the odds ratios directly*/
/*We can also output odds ratios from the logit command using the or option (default logit output is regression coefficients*/
logistic death qsmk [pweight=sw_a], cluster(seqn) 


***************************************************************
*PROGRAM 12.6
*Assessing effect modification by sex using a marginal structural mean model
*Data from NHEFS
*Section 12.5
***************************************************************/
drop pd_qsmk pn_qsmk sw_a

/*Check distribution of sex*/
tab sex

/*Fit logistc model for the denominator of IP weights, as in PROGRAM 12.3 */
logit qsmk sex race c.age##c.age ib(last).education c.smokeintensity##c.smokeintensity ///
c.smokeyrs##c.smokeyrs ib(last).exercise ib(last).active c.wt71##c.wt71 
predict pd_qsmk, pr

/*Fit logistic model for the numerator of IP weights, no including sex */
logit qsmk sex
predict pn_qsmk, pr

/*Generate IP weights as before*/
gen sw_a=.
replace sw_a=pn_qsmk/pd_qsmk if qsmk==1
replace sw_a=(1-pn_qsmk)/(1-pd_qsmk) if qsmk==0

summarize sw_a

/*Fit marginal structural model in the pseudopopulation, including interaction term between quitting smoking and sex*/
regress wt82_71 qsmk##sex [pw=sw_a], cluster(seqn)



/***************************************************************
PROGRAM 12.7
Estimating IP weights to adjust for selection bias due to censoring
Data from NHEFS
Section 12.6
***************************************************************/

clear
use "nhefs.dta"

/*Analysis including all individuals regardless of missing wt82 status: N=1629*/
/*Generate censoring indicator: C = 1 if wt82 missing*/
gen byte cens = (wt82 == .)

/*Check distribution of censoring by quitting smoking and baseline weight*/
tab cens qsmk, column
bys cens: summarize wt71

/*Fit logistic regression model for the  denominator of IP weight for A*/
logit qsmk sex race c.age##c.age ib(last).education c.smokeintensity##c.smokeintensity ///
c.smokeyrs##c.smokeyrs ib(last).exercise ib(last).active c.wt71##c.wt71 
predict pd_qsmk, pr

/*Fit logistic regression model for the  numerator of IP weights for A*/
logit qsmk
predict pn_qsmk, pr

/*Fit logistic regression model for the  denominator of IP weights for C, including quitting smoking*/
logit cens qsmk sex race c.age##c.age ib(last).education c.smokeintensity##c.smokeintensity ///
c.smokeyrs##c.smokeyrs ib(last).exercise ib(last).active c.wt71##c.wt71 
predict pd_cens, pr

/*Fit logistic regression model for the  numerator of IP weights for C, including quitting smoking */
logit cens qsmk
predict pn_cens, pr

/*Generate the stabilized weights for A (sw_a)*/
gen sw_a=.
replace sw_a=pn_qsmk/pd_qsmk if qsmk==1
replace sw_a=(1-pn_qsmk)/(1-pd_qsmk) if qsmk==0

/*Generate the stabilized weights for C (sw_c)*/
/*NOTE: the conditional probability estimates generated by our logistic models for C represent the conditional probability of being censored (C=1)*/
/*We want weights for the conditional probability of bing uncensored, Pr(C=0|A,L)*/
gen sw_c=.
replace sw_c=(1-pn_cens)/(1-pd_cens) if cens==0

/*Generate the final stabilized weights and check distribution*/
gen sw=sw_a*sw_c
summarize sw

/*Fit marginal structural model in the pseudopopulation*/
regress wt82_71 qsmk [pw=sw], cluster(seqn)
