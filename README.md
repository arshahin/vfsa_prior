# vfsa_prior

This package demonstrate how one can incorporate the priori information into VFSA optimization algorithm. 

Simply run this code "vfsa_prior_v02.m" and it will find all the dependencies.

Contribution are welcome to include other types of prior distribution functions (PDF) and other objective fucntions.
Converting the codes to other programming languages are also welcome.

The current objective function is a 1D finction with a global minimum at (0,0).

There are 7 types of prior PDF are included and you can select by setting prior_type to one the follwing options between 1 to 7. 

%%%%%%%% 1 for A unimodal Prior with low uncertainty :  sharp Gaussian N(0,1)
%%%%%%%% 2 for A unimodal Prior with high uncertainty: broad Gaussian N(0,3)
%%%%%%%% 3 for A deviated unimodal Prior simulating : broad Gaussian N(2,3)
%%%%%%%% 4 for A equiprobable bimodal Prior : N(0,1)+N(-4,1)
%%%%%%%% 5 for A bimodal Prior simulating high probability on target : N(0,1)+0.5*N(-4,1)
%%%%%%%% 6 for A bimodal Prior simulating low probability on target : 0.5*N(0,1)+N(-4,1)
%%%%%%%% 7 for An erroneous unimodal Prior  : N(6,1)


You can also customize the way that prior information will be added to VFSA duriung iteration. Currently the follwoing options are included and they can be selected by setting weight_type as follows:
%%%%%%%% 1 for pure VFSA (no prior is applied)
%%%%%%%% 2 for Polynomial increase (with increasing the iteration number weights are increasing with a 4th order Polynomial)
%%%%%%%% 3 linear incease 
%%%%%%%% 4 for exponential
%%%%%%%% 5 fixed weigth 0.50 


Please make sure to check the PDF file provided to give you an idea of different visualizations that you can reproduce.



