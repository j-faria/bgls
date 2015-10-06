A Bayesian formalism for the generalised <br/> Lomb-Scargle periodogram
=======================================================================

The BGLS tool calculates the Bayesian Generalized Lomb-Scargle periodogram as described in Mortier et al. (2014). It is written in Python (tested on Python 2.7).

The code contains the definition of the algorithm, takes as input arrays with a time series, a dataset and errors on those data, and returns arrays with sampled periods and the periodogram values at those periods.

In order to run, it requires the following python packages:

    * numpy (http://www.numpy.org/)
    * mpmath (http://mpmath.org/)

More information can be found in the paper

![Mortier, A., Faria, J. P., Correia, C. M., Santerne, A., Santos, N. C. 2015, A&A, 573, A101](http://adsabs.harvard.edu/abs/2015A%26A...573A.101M)

If you use the code provided here, please cite the above paper. We welcome feedback via the [issues](https://github.com/j-faria/bgls/issues).

The code is distributed under the MIT license. See the LICENSE file.
