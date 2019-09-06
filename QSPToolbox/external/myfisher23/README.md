# myfisher23
Fisher's Exact Probability Test on 2x3 matrix.<br/>
Fisher's exact test of 2x3 contingency tables permits calculation of
precise probabilities in situation where, as a consequence of small cell
frequencies, the much more rapid normal approximation and chi-square
calculations are liable to be inaccurate. The Fisher's exact test involves
the computations of several factorials to obtain the probability of the
observed and each of the more extreme tables. Factorials growth quickly,
so it's necessary use logarithms of factorials. This computations is very
easy in Matlab because x!=gamma(x+1) and log(x!)=gammaln(x+1). This
function is fully vectorized to speed up the computation.

Syntax: 	p=myfisher23(x)
     
    Inputs:
          X - 2x3 data matrix 
    Outputs:
          - Three p-values

          Created by Giuseppe Cardillo
          giuseppe.cardillo-edta@poste.it

To cite this file, this would be an appropriate format:
Cardillo G. (2007) MyFisher23: a very compact routine for Fisher's exact
test on 2x3 matrix
http://www.mathworks.com/matlabcentral/fileexchange/15399
