A Bayesian formalism for the generalised <br/> Lomb-Scargle periodogram
=======================================================================

The BGLS tool calculates the Bayesian Generalized Lomb-Scargle periodogram as described in Mortier, Faria et al. (2014). It is written in Python (tested on Python 2.7).

The code contains the definition of the algorithm, takes as input arrays with a time series, a dataset and errors on those data, and returns arrays with sampled periods and the periodogram values at those periods.

In order to run, it requires the following python packages:

    * numpy (http://www.numpy.org/)
    * mpmath (http://mpmath.org/)

More information can be found in Mortier, Faria et al. (2014)
"" LINK TO PAPER TO BE INSERTED ""
