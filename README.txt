    Title: Model-robust designs for nonlinear quantile regression.

    Authors: Selvakkadunko Selvaratnam, Linglong Kong, and Douglas P. Wiens (University of Alberta))

 The codes are written in MATLAB. We have two main functions: (i) finalmainadaptive and (ii) nonlinearsequential. ALso (iii) compare
    
 (i) finalmainadaptive.m: The function finalmainadaptive will generate results for the adaptive designs that are described in Section 3.2. These results are a table for the average losses, a plot for the adaptive designs, and a plot for the comparisons between true and estimated scale functions. We consider the two types of initial designs: (a) 2-point design, (b) the equispaced design. So, we have to choose one method and deactivate the other method in the function finalmainadaptive.

 (ii) nonlinearsequential.m: The function nonlinearsequential will produce results for the sequential local designs that are described in Section 3.1. These results are tables for losses, and a plot for the sequential local designs.

 (iii) compare.m: Prepares the final two figures in the article, using the values in 'tables.pdf', which is also in the .zip file.