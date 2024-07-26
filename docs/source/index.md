<p align="center" style="margin-top: 40px; margin-bottom: 40px;">
  <img src="_static/img/pybmds.png" width="210px" style="margin-right: 10px">
  <img src="_static/img/bmds-desktop-logo.png" width="210px">
</p>

**pybmds** is a Python package for executing the U.S. EPA Benchmark Dose Modeling Software (BMDS). The package includes dose-response models for multiple types of dose-response data, including dichotomous, continuous, nested dichotomous, and tumor (as well as multitumor modeling). The pybmds library is designed for those familiar with scripting to conduct dose-response modeling.

**BMDS Desktop** is an application for running dose-response modeling using a Graphical User Interface (GUI) locally. The **bmds-ui** Python package installs BMDS Desktop; it can also be used to deploy a copy of the [BMDS Online](https://bmdsonline.epa.gov) application.

:::{note}

This software/application has been approved for release by the U.S. Environmental Protection Agency (USEPA). Although the software has been subjected to rigorous review, the USEPA reserves the right to update the software as needed pursuant to further analysis and review. No warranty, expressed or implied, is made by the USEPA or the U.S. Government as to the functionality of the software and related material nor shall the fact of release constitute any such warranty. Furthermore, the software is released on condition that neither the USEPA nor the U.S. Government shall be held liable for any damages resulting from its authorized or unauthorized use.

:::

This guide describes the use and installation of both **pybmds** and **BMDS Desktop**.

:::{admonition} NOTE - this will not work until we officially release!
:class: caution

The content below below will not work until we are cleared for public release; please follow the detailed [installation guide](./installation.md). This message will be removed prior to official release.

---

To install **BMDS Desktop** (which includes `pybmds`):

```bash
pip install bmds-ui
```
:::

Please follow the [installation guide](./installation.md) if you plan on seriously using `pybmds` or BMDS Desktop for work on multiple projects; the guide makes it easier to setup multiple versions, and describes in detail how upgrade or uninstall existing versions, and describes possible issues and solutions that may arise during installation.

**Highlights:**

* Dose response modeling for multiple dataset types (continuous, dichotomous, nested dichotomous, multitumor)
* Plotting and summary tables capabilities
* Model recommendation logic
* Reporting in standard Microsoft Excel and Microsoft Word reports
* Batch execution for multiple datasets or option set configuration

The BMDS Desktop main interface for dose response analyses:

```{figure} _static/img/bmds-desktop.jpg
:alt: Screenshot of BMDS Desktop Application
```

Output showing a Hill model fit and summary table:

```{figure} _static/img/bmds-output.jpg
:alt: Screenshot of BMDS Desktop Dose Response Output
```

Example use the `pybmds` software to run a dose-response session:

```python
import pybmds

dataset = pybmds.DichotomousDataset(
   doses=[0, 10, 50, 150, 400],
   ns=[25, 25, 24, 24, 24],
   incidences=[0, 3, 7, 11, 15],
)
session = pybmds.Session(dataset=dataset)
session.add_default_models()
session.execute()
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
   Github pybmds <https://github.com/USEPA/bmds>
   Github bmds-ui <https://github.com/USEPA/bmds-ui>
```
