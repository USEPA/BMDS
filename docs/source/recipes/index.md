# pybmds Recipes

The examples below describe some key capabilities of `pybmds`.

## Modeling

```{toctree}
:maxdepth: 2

dichotomous
dichotomous_ma
multitumor
continuous
nested_dichotomous
batch
```

## Data Manipulation

```{toctree}
:maxdepth: 2

preparing-datasets
custom-excel-exports
using-r
```

## Jupyter Notebooks

Recipes are written using [Jupyter](https://jupyter.org/) notebooks. There are many excellent tutorials online describing how to use Jupyter notebooks and this is beyond the scope of `pybmds`. However, a few tips are described below for new Jupyter users.

### Rendering figures

The plotting functionality of `pybmds` uses the [matplotlib](https://matplotlib.org/) package to generate figures. To render matplotlib figures within a Jupyter notebook, add this "magic function" to a cell at the top of the notebook, and execute the cell to enable inline plot rendering:

```python
%matplotlib inline
```

After running this cell, matplotlib figures will appear inline after execution, if it is the final output in a cell. For more information, see this [summary](https://stackoverflow.com/q/43027980/906385) from Stack Overflow.

### Displaying extra output

By default, notebook cells only display the last output of the code that is executed in them, and anything that is printed to standard out (for example anything using the `print()` function). You can also display content at any time in a cell by using the [display](https://ipython.readthedocs.io/en/stable/api/generated/IPython.display.html#IPython.display.display) method. This will allow you to show plots or text in the middle of loops, or show multiple outputs from a single cell. We utilize this function in some recipes.

```python
from IPython.display import display

display(...)
```

Anything that is passed to this function will be shown after the cell is executed.
