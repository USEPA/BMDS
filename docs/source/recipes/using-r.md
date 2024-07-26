# Using R

You can use use RStudio to execute `pybmds` using the [reticulate](https://rstudio.github.io/reticulate/) package; follow the standard [installation](../installation.md) guide and install within a Python environment. After a successful installation, load the reticulate package and activate the Python environment with `pybmds` installed:

```R
library(reticulate)

# Use one of the commands below, depending on environment:
# ... if using Python, the path to the virtual environment python executable
use_python('~/dev/bmds-desktop/bin/python')
# ... if using Anaconda, the environment name
use_condaenv('bmds-desktop')
```

Load `pybmds` and then execute a logistic model fit:

```R
# import the pybmds package and show version information
pybmds <- import('pybmds')
pybmds$citation()

# fit a logistic model to a dataset
ds = pybmds$DichotomousDataset(
    doses = c(0, 10, 50, 150, 300),
    ns = c(25, 25, 24, 24, 24),
    incidences = c(0, 3, 7, 11, 15)
)
model = pybmds$models$dichotomous$Logistic(
    dataset=ds,
    settings=dict(bmr=0.10)
)
model$execute()
```

After execution, you can view model results or generate a plot:

```R
paste0("BMD: ", model$results$bmd)
plot(
    model$results$plotting$dr_x,
    model$results$plotting$dr_y,
    type='l'
)
```


To show a plot generated in Python in RStudio:

```R
showPyFig <- function(fig){
    tmp <- tempfile(fileext='.png')
    fig$savefig(tmp)
    grid::grid.raster(png::readPNG(tmp))
}

fig = model$plot()
showPyFig(fig)
```

