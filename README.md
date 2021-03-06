
# True and Error Models

## R program to analyze True and Error Models - TEMAP2.R
These instructions will help you run the program presented in Birnbaum, M. H., & Quispe-Torreblanca, E. G. (2018). 
TEMAP2: Computer program to analyze True and Error Models. Published in [Judgment and Decision Making](
http://www.sjdm.org/journal/18/18507/jdm18507.pdf)


The Excel file can be downloaded from the [Judgment and Decision Making website](
http://www.sjdm.org/journal/vol13.5.html)

The Excel file used for input to the program contains three worksheets.

The “READ ME” sheet contains information on how to organize the data and how to specify the inputs. 

The “Inputs” worksheet contains values that can be adjusted by the user to request that the program analyze one or more
of TE-4, TE-2, or TE-1, to specify either G or χ2, to allow parameters to be free or fixed, to request Monte Carlo simulations, and to request bootstrapping analysis. 

The notation in the program is a bit different from the notation used in the paper for the example of the Allais paradox: The notations 0 and 1 are used to denote first or second responses, respectively, which in this example are choice of the “risky (R)” or “safe (S)” alternatives. The notations, a_00, a_01, a_10, and a_11 in the program correspond to pRR′, pRS′, pSR′, and pSS′, respectively. To set up EU, fix the values of a_01 and a_10 to 0.

The “participant responses” worksheet contains the data: response frequencies (counts) for the 16 possible response patterns. 
The first column is the case number. The first case of Example.xlsx contains the data of Table 2. 
Several cases can be analyzed in the same computer run. Each line represents a different case, which may represent aggregated data for a group of participants (for gTET) or for an individual (iTET).

The R-program is documented by many comments, statements on a line following #. 
The section beginning with line 2330 is used to create the output of the program, and this section should be easiest to modify by those familiar with R. Comments include suggestions for revising this section to access additional information generated by the program. The example file, Example.xlsx, has been configured to request only 100 samples to illustrate the program. Once the program is running properly, the 100 on the “Inputs” worksheet of Example.xlsx can be changed to a higher value for better accuracy of Monte Carlo simulation and bootstrapping. 

Sample results in the paper are based on the setup in Example.xlsx, for the first case (the data in Table 2), except using 10000 instead of 100.

Some examples of outputs of the program can be seen [here](https://github.com/edikaQT/TE-model/blob/master/TEMAP_example.md).
 
## True and Error Models for testing transitivity violations - TE_analysis.R

These instructions will help you run the program presented in Birnbaum, M. H., Navarro-Martinez, D., Ungemach, C., Stewart, N. & Quispe-Torreblanca, E. G. (2016). Risky decision making: Testing for violations of transitivity predicted by an editing mechanism. 
Published in [Judgment and Decision Making](
http://journal.sjdm.org/15/15615a/jdm15615a.pdf)

The Excel file can be downloaded from the [Judgment and Decision Making website](
http://journal.sjdm.org/vol11.1.html)

To run the program, select the program, TE_analysis.R, as the source file (from the File menu) or type the command as follows (with the appropriate path to the program) in the Rconsole:
```
source("<folder>TE_analysis.R")
```
The program generates 21 new files: 15 figures for the first listed case (11 distributions of the bootstrapped estimates of the TE model, and 4 distributions of the test statistic
for Monte Carlo samples, for the TE model and Independence model by either conservative or re-fit method).
It also creates 2 comma separated values (CSV) files (one
for each model) with the parameter estimates for each case along with its index of fit and conventional p-value. Two additional
CSV files are created for the simulated Monte Carlo p-values by each of the two methods for each model. A CSV
file with bootstrapped parameter estimates for the TE model is also saved, as well as a text file containing other output. 

It takes about 15 min. on a MacBook Pro to run the program with example.txt (6 cases with 1000 Monte Carlo simulations
and 1000 Bootstrapped samples for each case). One can set the number of samples to 10,000 for greater accuracy,
once the program is running correctly

