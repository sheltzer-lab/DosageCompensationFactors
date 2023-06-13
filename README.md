# DosageCompensationFactors

## Requirements

* R (>= 4.3.0)

## Usage

All scripts are to be called from this directory (location of this `README.md` file) as the working directory.

R-Scripts can be executed with an IDE of your choice (RStudio, PyCharm, ...) or by executing 
```bash
RScript ./path/to/script.R
```

### Package installation

`Renv` is used in this project to keep track of all required packages.
Install the `renv` package using te following command in an R console:

```r
install.packages("renv")
```

To restore the environment and install all necessary packages execute the following command in the R console
from the directory where `renv.lock` is located at:

```r
renv::restore()
```
