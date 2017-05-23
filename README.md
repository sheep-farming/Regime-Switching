# README

Regime Switching
*An UZH Asset Management: Advanced Portfolio Management FS17 Project*

## Author

Kun Yu       kun.yu@uzh.ch 16-704-389
Luoyi Zou    luoyi.zou@uzh.ch 16-743-536

## File List

### Data Files 

##### Historical Return Data
- data.xls
 - The original data file with data retrieved from **BenchmarkReturns_v2.xlsx** *(provided by Dr. Rohner)* and Bloomberg
- dataCSV.csv
 - Same set of data, **data.xls** cleaned for the purpose of being imported to MATLAB
- new.mat
 - Same set of data stored in MATLAB-friendly .mat format

##### Completed Estimation Results
- est.mat
 - It takes tens of minutes to run the **mainEstimation.m** script to estimate all parameters, so we stored the resulting variables that can be simply loaded by **backTesting.m** to process the remaining calculation & run backtest.
 - It is also possible to run the **mainEstimation.m** to go through from the beginning of the process.

### MATLAB Scripts & Functions
- mainEstimation.m
 - The entry point of the project. Estimate the P, Q and other regime-related parameters using the MS_Regress-MATLAB package.
- getBetas.m
 - Called by **mainEstimation.m** to calculate the beta coefficients of assets for every observation windows.
- backTesting.m
 - The script that calculates asset-and-regime-specific expected returns and volatilities. It also computes the optimal weights of assets across months, and calculate the overall portfolio return of the strategy.
- plotting.m
 - Called by **backTesting.m** by the end of the script. Used to visualise the data estimated and calculated.

## Acknowledgements

- Used in **mainEstimation.m**
 - **Marcelo Perlin** - *MS_Regress - The MATLAB Package for Markov Regime Switching Models* - [SSRN](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=1714016) - [GitHub Page](https://github.com/msperlin/MS_Regress-Matlab)


