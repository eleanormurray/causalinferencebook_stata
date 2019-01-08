/***************************************************************
Stata code for Causal Inference by Miguel Hernan & Jamie Robins
Date: 1/8/2019
Author: Eleanor Murray 
For errors contact: emurray@mail.harvard.edu
***************************************************************/

/*********************
*Data preprocessing***
**********************/

clear
use "nhefs.dta"
gen byte cens = (wt82 == .)


/***************************************************************
# PROGRAM 14.1
# Ranks of extreme observations
# Data from NHEFS
***************************************************************/

/*Ranking of extreme observations*/
extremes wt82_71 seqn

/*Estimate unstabilized censoring weights for use in g-estimation models*/
glm cens qsmk sex race c.age##c.age ib(last).education c.smokeintensity##c.smokeintensity c.smokeyrs##c.smokeyrs ib(last).exercise ib(last).active c.wt71##c.wt71
predict pr_cens
gen w_cens = 1/(1-pr_cens)
summarize w_cens

/*Analyses restricted to N=1566*/
drop if wt82 == .
summarize wt82_71

/***************************************************************
# PROGRAM 14.2
# G-estimation of a 1-parameter structural nested mean model
# Brute force search
# Data from NHEFS
# Author: Ehsan Karim ehsan@stat.ubc.ca
***************************************************************/

/*Generate test value of Psi = 3.446*/
gen psi = 3.446

/*Generate H(Psi) for each individual using test value of Psi and their own values of weight change and smoking status*/
gen Hpsi = wt82_71 - psi * qsmk 

/*Fit a model for smoking status, given confounders and H(Psi) value, with censoring weights and display H(Psi) coefficient*/
logit qsmk sex race c.age##c.age ib(last).education c.smokeintensity##c.smokeintensity c.smokeyrs##c.smokeyrs ib(last).exercise ib(last).active c.wt71##c.wt71 Hpsi [pw = w_cens], cluster(seqn)
di _b[Hpsi]

drop psi Hpsi

/***************************************************************
# G-estimation: Checking multiple possible values of psi
***************************************************************/
drop psi Hpsi

local seq_start = 2
local seq_end = 5
local seq_by = 0.1
local seq_len = (`seq_end'-`seq_start')/`seq_by' + 1
matrix results = J(`seq_len',3,0)
matrix list results
gen psi = .
gen Hpsi = .
local j = 0
forvalues i =  `seq_start'(`seq_by')`seq_end' {
	local j = `j'+1
	replace psi = `i'
	replace Hpsi = wt82_71 - psi * qsmk 
	quietly logit qsmk sex race c.age##c.age ib(last).education c.smokeintensity##c.smokeintensity c.smokeyrs##c.smokeyrs ib(last).exercise ib(last).active c.wt71##c.wt71 Hpsi [pw = w_cens], cluster(seqn)
	local b = _b[Hpsi]
	di "coeff `b' is generated from psi `i'"
	matrix results[`j',1]= `i'
	matrix results[`j',2]= `b'
	matrix results[`j',3]= abs(`b')
}
matrix colnames results = "psi" "B(Hpsi)" "AbsB(Hpsi)"
mat li results 

mata
res = st_matrix("results")
for(i=1; i<= rows(res); i++) { 
if (res[i,3] == colmin(res[,3])) res[i,1]	
}
end
* Setting seq_by = 0.01 will yield the result 3.46

/***************************************************************
# PROGRAM 14.3
# G-estimation for 2-parameter structural nested mean model
# Closed form estimator
# Data from NHEFS
# Author: Ehsan Karim ehsan@stat.ubc.ca
***************************************************************/


/***************************************************************
# G-estimation: Closed form estimator linear mean models  
***************************************************************/

logit qsmk sex race c.age##c.age ib(last).education c.smokeintensity##c.smokeintensity c.smokeyrs##c.smokeyrs ib(last).exercise ib(last).active c.wt71##c.wt71 [pw = w_cens], cluster(seqn)  
predict pr_qsmk
summarize pr_qsmk

ssc inst tomata
putmata *, replace
mata: diff = qsmk - pr_qsmk
mata: num = w_cens :* wt82_71 :* diff
mata: den = sum(w_cens :* qsmk :* diff)
mata: sum( num/den )

/***************************************************************
# G-estimation: Closed form estimator for 2-parameter model
***************************************************************/
mata
diff = qsmk - pr_qsmk
diff2 = w_cens :* diff

lhs = J(2,2, 0)
lhs[1,1] = sum( qsmk :* diff2)
lhs[1,2] = sum( qsmk :* smokeintensity  :* diff2 )
lhs[2,1] = sum( qsmk :* smokeintensity :* diff2)
lhs[2,2] = sum( qsmk :* smokeintensity :* smokeintensity :* diff2 )
                                                                
rhs = J(2,1,0)
rhs[1] = sum(wt82_71 :* diff2 )
rhs[2] = sum(wt82_71 :* smokeintensity :* diff2 )

psi = (lusolve(lhs, rhs))'
psi
psi = (invsym(lhs'lhs)*lhs'rhs)'
psi
end



