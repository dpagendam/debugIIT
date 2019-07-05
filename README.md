## debugIIT

## Stochastic modelling of Wolbachia Incompatible Insect Technique (IIT) programs for control of mosquito populations.
**Authors**: Dan Pagendam (CSIRO)

**Contributors**: Nigel Snoad (Verily), Nigel Beebe (University of Queensland), Brendan Trewin (CSIRO)

debugIIT provides a set of functions for stochastic modelling of population dynamics using continuous-time Markov chains (birth-death processes).  The package is useful for stochastic simulation of possible population trajectories, which when aggregated, can be used to construct probabilistic predictive distributions for IIT program outcomes such as the likelihood of driving a population to extinction. The simulation code is largely written in C++ and makes use of the Rcpp package.

The package was developed at [CSIRO](http://www.csiro.au), Australia as part of collaboration with Verily's [Debug Project](https://debugproject.com/), targeting the mosquito _Aedes aegypti_.


### Package installation
First, clone the package from this repository using 

``` git clone https://github.com/dpagendam/debugIIT ```

You will be prompted to enter your GitHub username and password to complete the clone.

To install the package then type

```	
	R CMD build debugIIT
	R CMD INSTALL debugIIT_1.0.tar.gz
```



### Using this package

Please see the example scripts in the Example Scripts folder for examples of using the model under three different release scenarios: (i) constant releases (same numbers released each release day); (ii) adaptive releases, where Wolbachia release numbers are dynamically altered to maintain a costant overflooding level; and (iii) crude adaptive releases, where the release numbers are reduced once the population reaches 50% of its initial size and then again once the population reaches 10% of its initial size.
