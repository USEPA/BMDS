<p align="center" style="margin-top: 40px; margin-bottom: 40px;">
  <img src="_static/img/pybmds.png" width="210px">
</p>

**pybmds** is a Python package for executing the U.S. EPA Benchmark Dose Modeling Software (BMDS). The package includes dose-response models for multiple types of dose-response data, including dichotomous, continuous, nested dichotomous, and tumor (as well as multitumor modeling). The **pybmds** library is designed for those familiar with basic scripting techniques to conduct dose-response modeling, however it can be extended using the **bmds-ui** library to run BMDS via a graphical user interface (GUI) known as BMDS Desktop. The **bmds-ui** library is also used for [BMDS Online](https://bmdsonline.epa.gov); BMDS Desktop allows for users to use the same methods, but offline and locally on their desktop.

This guide describes the use and installation of both **pybmds** and **BMDS Desktop**.

:::{admonition} NOTE - this will not work until we officially release!
:class: error

The content below below will not work until we are cleared to release the the software on [pypi](https://pypi.org/); please follow the [EPA installation guide](./installation.md#internal-epa) until we can release formally on the internet.  This message will be removed prior to official release.

---

To install **BMDS Desktop** (which includes pybmds):

```bash
pip install bmds-ui
```

To install just **pybmds**:

```bash
pip install bmds
```
:::

:::{tip}

We recommend using the [installation guide](./installation.md) if you plan on seriously using pybmds or BMDS Desktop for work on multiple projects; the guide makes it easier to setup multiple versions, and describes in detail how upgrade or uninstall existing versions, and describes possible issues and solutions that may arise during installation.
:::

**Highlights:**

* Dose response modeling for multiple dataset types
* Batch execution for multiple datasets or option set configuration
* Reporting in standard Microsoft Excel and Microsoft Word reports

A detailed guide for installation is described

# Quickstart

```python
import pybmds

# create a dataset
dataset = pybmds.DichotomousDataset(
    doses=[0, 10, 50, 150, 400],
    ns=[25, 25, 24, 24, 24],
    incidences=[0, 3, 7, 11, 15],
)

# create a BMD session
session = pybmds.Session(dataset=dataset)

# add all default models
session.add_default_models()

# execute the session
session.execute()

# show a summary figure
session.plot(colorize=True)
```


# Contents

```{eval-rst}
.. toctree::
   :maxdepth: 2

   installation
   desktop
   quickstart
   reference/index
   recipes/index
   history
```

```{eval-rst}
.. toctree::
   :caption: Links
   :maxdepth: 2

   US EPA BMDS <https://epa.gov/bmds>
   pybmds Github <https://github.com/USEPA/bmds>
   bmds-ui Github <https://github.com/USEPA/bmds-ui>
```
