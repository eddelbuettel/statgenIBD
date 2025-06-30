# statgenIBD 1.0.9

-   Functions no longer rely on soft-deprecated ggplot2 functions.

# statgenIBD 1.0.8

-   `readRabbit` now reads input files produced by the open source version of RABBIT, available at <https://github.com/Biometris/RABBIT>

# statgenIBD 1.0.7

-   Plotting didn't work correctly for populations where there was a heterozygous probability in the IBD output. This is now corrected for plotType "allGeno", "meanProbs", and "totalCoverage".
-   It is now possible to create a plot of plotType "singleGeno" when the position of the markers in the map have large values, i.e. when they are in base pairs.
-   A small bug in the concatenation of `calcIBD` objects with an identical number of genotypes is fixed.

# statgenIBD 1.0.6

-   A bug in the `plot` function is fixed. The "singleGeno" plots now show the correct probabilities.
-   Functions no longer rely on soft-deprecated ggplot2 functions.

# statgenIBD 1.0.5

-   Two new functions for writing (`writeIBDs`) and reading IBD (`readIBDs`) probabilities to a custom format have been added. The format is explained in a new vignette `vignette("IBDFileFormat", package = "statgenIBD")`.
-   Two new plots are added for objects of class `IBDprob`, a plot showing the mean coverage of each parent across the population (`plotType = "meanProbs"`) and a plot showing the total coverage of each parent (`plotType = "totalCoverage"`).
-   Parent names in the output of `calcIBD` no longer have a leading "p".
-   The function `readRabbit`, for reading output from the RABBIT software for computing IBDs for complex crosses, is moved here from `statgenMPP`.

# statgenIBD 1.0.4

-   Increased speed of IBD calculations by improving an inefficient part of the algorithm.

# statgenIBD 1.0.3

-   It is now possible to plot the pedigree of an object of class `IBDprob` by specifying `plotType = "pedigree"`.
-   A small bug in the calculation of probabilities for populations C3 and C4 is fixed.

# statgenIBD 1.0.2

-   It is now possible to specify chromosomes in character format, i.e. 1a, 1b, 1c as is common in certain crops.
-   The plot function for calculated IBD probabilities now has an extra option `plotType` that can have values "singleGeno", for the original plot, and "allGeno" for a new plot showing probabilities for all genotypes in a single plot.
-   A bug that in some cases could lead to misnaming of the markers when combining multiple imputed populations is fixed.

# statgenIBD 1.0.1

-   A bug where written Flapjack files could not be read is fixed.

# statgenIBD 1.0.0

-   Initial CRAN version
