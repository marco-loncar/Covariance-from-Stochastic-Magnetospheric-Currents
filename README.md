# The Impact of Magnetospheric Currents on Planetary Magnetic Field Reconstruction

This repository contains data, programs and field models used in the work 'The Impact of Magnetospheric Currents on Jovian Magnetic Field Reconstruction' (Loncar & Jackson, 2026). Namely, this includes the augmented Juno dataset used in field reconstructions, the machinery used to calculate covariance between measurements as a result of stochastic magnetospheric currents and the resulting spherical harmonic field models. 

* The data processing of publically available Juno data (found on the NASA Planetary Data System (Connerney, 2024)) is outlined in appendix A of Loncar & Jackson (2026) 
* The mathematical background behind the evaluation of covariance functions is similarly contained in the main body and appendix B of the same paper 
* The resulting models are presented and discussed throughout this paper also

The dataset, code and models listed here were assembled by Marco Loncar and any questions regarding their use or development will be happily answered upon request (marco.loncar@eaps.ethz.ch).

The contents of this repository are: 
    
* Folder entitled **Dataset** containing Juno magnetic field data files **Orbit\#\#.txt** used in field reconstructions (where \#\# refers to the orbit number)

    * Spherical component data $B_r$, $B_{\theta}$, $B_{\varphi}$
    * Up to orbit 70 of Juno, omitting orbits 2 and 19
    * Restricted to within $2R_J$ of Jupiter's centre
    * Despiked and offset at dynamic range changes (described in Loncar & Jackson, 2026)
    * Downsampled to $30s$ intervals
    * Organised into columns of measurement: time, position, value, uncertainty and instrument range

* Program entitled **FindingCovariance.py** 

    * Contains Python functions to find the covariance between magnetic field measurements as a result of stochastic magnetospheric currents, $C_{ij}$ as described in Loncar & Jackson (2026)
    * The variable naming conventions here match those used in the paper apart from using $R_P$ instead of $R_J$ (to offer a more general *planetary* description as opposed to a *Jupiter*-specifc case)

* Folder entitled **Models** containing two spherical harmonic, Jovian magnetic field reconstructions 

    * Includes a reference model with no consideration of stochastic currents **JREF.txt** and a model that accounts for stochastic currents **JSTO.txt**
    * The text files first note: spherical harmonic truncation degrees $L_{\text{int}}$ and $L_{\text{ext}}$, damping parameter $\lambda / \text{nT}^2$ and current correlation scaling parameter $\beta / \text{A m}^{-\frac{1}{2}}$
    * Model parameters (or Gauss coefficients) $g_l^m$, $h_l^m$, $G_l^m$, $H_l^m$ are listed in canonical order (second and third columns corresponding to $l$ and $m$)

For further information regarding the use of this data, code or models in the construction of covariance matrices, figures or other, please feel free to reach out. 
