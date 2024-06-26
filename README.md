# _Causal Inference_ Stata Code

This repo contains Stata code for the book Causal Inference: What If, by Miguel Hernán and James Robins [book site](https://www.hsph.harvard.edu/miguel-hernan/causal-inference-book/). 

The code here corresponds to the SAS programs found at the book site.

## Usage notes
1. For all bootstrapped confidence interval code, the corresponding point estimate code needs to be run before the bootstrap code. This is because the bootstrap code only calculates the values in the bootstrapped data, but as output shows the point estimate and intervals. If you receive the error message "observe not found" this means you have not run the point estimate code.
2. Chapter 17, program 17.5 is not included in the current code set due to the complexity of the code. I am working on finalizing the code for program 17.5 -- anyone interested in helping with this code can contact me for the draft code. This note will be removed once program 17.5 is added.

## Data

The data can be obtained from the [book site](https://www.hsph.harvard.edu/miguel-hernan/causal-inference-book/). The do files all assume that the Stata version of the data (.dta) has been saved in the same directory as the do files. To automatically set the working directory, open Stata by directly double-clicking on the do files. 

If you open the do files from within Stata you will need to set the working directory manually: before running the command: use "nhefs.dta" add your path (eg. use "C:\Users\Downloads\nhefs.dta")

## Author

Ellie Murray
