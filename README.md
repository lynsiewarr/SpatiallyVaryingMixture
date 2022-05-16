# spatiallyvaryingmixture
Data and Code associated with the paper Distributional Validation of Precipitation Data Products with Spatially Varying Mixture Models (2022 in JABES).

SVMixFixed_SameComp.R contains the code to fit the spatially varying mixture model with component distributions the same across different data products.

SVAphro.R contains an example implementation of the function in SVMixFixed_SameComp.R .

AMCMCUpdate.R contains a function for adaptive MCMC that is utilized in SVMixFixed_SameComp.R .

aphro.Rdata , era5.Rdata , merra2.Rdata , and trmm.Rdata are the data files used in the paper.
