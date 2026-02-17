# Contents of *Dataset* folder

As described in the higher level README file, this folder contains processed Juno magnetometer data. The data is split into individual orbit files, each with nine columns which are organised as follows: 

* *Decimal day (\#\#\#\#)* 

> Time of measurement given in decimal days of the year \#\#\#\#

* *r / R\_J*

> Radial position of Juno satellite from Jupiter's centre (units of Jovian radii $R_J$)

* *theta / radians*

> Colatitudinal position of Juno satellite in System III spherical coordinates (units of radians)

* *phi / radians* 

> Azimuthal position of Juno satellite in System III spherical coordinates (units of radians)

* *B\_r / nT* 

> Radial magnetic field strength $B_r$ measured by Juno (units of nT)

* *B\_t / nT* 

> Colatitudinal magnetic field strength $B_{\theta}$ measured by Juno (units of nT)

* *B\_p / nT* 

> Azimuthal magnetic field strength $B_{\varphi}$ measured by Juno (units of nT)

* *Error / nT* 

> Uncertainty in the Juno measurement, a combination of attitude, quantisation and baseline uncertainties (units of nT)

* *Range*

> Current dynamic measurement range of the Juno magnetometer (between 4 and 6 for the regions of interest here)
