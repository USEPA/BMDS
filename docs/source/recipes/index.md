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

## Using Jupyter

Recipes are written using [Jupyter](https://jupyter.org/) notebooks for running data-science code in Python. There are many excellent tutorials online describing how to use Jupyter notebooks and this is beyond the scope of `pybmds`. Below are a few tips that may be helpful when copying text from the recipes.

### Rendering figures

The build in plotting functionality of `pybmds` uses [matplotlib](https://matplotlib.org/). To render matplotlib figures within a Jupyter notebooks, add a this line to the top of the notebook[^1]:

```python
%matplotlib inline
```

After running this cell, matplotlib figures will appear inline after a cell is executed, if it is the final output in a cell.

[^1]: A Stack Overflow [summary](https://stackoverflow.com/q/43027980/906385) of the purpose of `%matplotlib inline`

### Display extra cell output

By default, notebook cells only display the last output of the code that is executed in them, and anything that you `print` to the screen. You can also `display`[^2] content at any time in a cell by importing the display method. This will allow you to show plots or text in the middle of loops, for example:

```python
from IPython.display import display

display(...)
```

[^2]: For more information on displaying outputs; see the IPython [documentation](https://ipython.readthedocs.io/en/stable/api/generated/IPython.display.html#IPython.display.display).

