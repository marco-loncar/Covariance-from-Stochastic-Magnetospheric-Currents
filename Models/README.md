# Contents of *Models* folder

As described in the higher level README file, this folder contains Jovian magnetic field reconstructions derived from the Juno data within the **Dataset** folder. The two model files contain the spherical harmonic model parameters (or Gauss coefficients) organised as follows: 

* *L\_int* and *L\_ext*

> The truncation degrees used for the internal and external fields

* *lambda* 

> The damping parameter $\lambda$ (units of $\text{nT}^2$)

* *beta* 

> The scaling parameter for stochastic magnetospheric currents $\beta$ (units of $\text{A}$ $\text{m}^{-1/2}$)

* Rows of model parameters in the canonical ordering, with columns indicating 

    * *g*, *h*, *G* or *H*
    
    > Even/odd (*g*/*h*) and internal/external (*g*/*G*) associated coefficients

    * *\#* *\#*
    
    > Degree *l* and order *m* of the model parameter respectively 
    
    * *\#\#\#\#\#\#.\#\#\#\#\#\#* 
    
    > Value of model parameter (units of $\text{nT}$) 
