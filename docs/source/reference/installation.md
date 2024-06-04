# Installation

Lake layers of an onion, there are a few components in the BMDS software suite:

1. `bmdscore` is the innermost layer which includes the dose-response models and fitting algorithms
2. `pybmds` is a wrapper around `bmdscore` that makes the methods easier to call and generates standardized analyses and reports
3. `BMDS Desktop` wraps `pybmds` and provides a user interface to setup analyses and visualize results

```{figure} ../_static/img/diagram.png
:alt: Overview of how BMDS Python ecosystem is organized

Conceptual diagram of how the bmds software packages are organized. The `pybmds` library embeds a copy of `bmdscore`, and has all the necessary functionality to execute BMDS modeling simulations and generate reports by writing Python code. The `BMDS Desktop` package is an optional additional package that allows you to run BMDS modeling simulations using a user-interface in your browser, with no coding required. Even though a web browser is used, all modeling occurs locally on your desktop computer.
```

If you're familiar with Python and packages (or R packages), then installing these packages should be should be familiar. If you're not, see the guides below.

(epa-guide)=
## Quick Start (EPA Guide)

:::{tip}
This portion of the guide is recommended as we develop pre-release versions. Recommended for EPA staff. Please check the FAQ below with any questions, or reach out to [Andy Shapiro](mailto:shapiro.andy@epa.gov) for any feedback.  If anything is missing from the guide or anything is confusing, please let us know!
:::

After installing Python, open up your terminal on your computer. Ensure that Python is on your path:

```bash
python -V
```

This should return a Python version (eg., `Python 3.12.3`).  Note that the Python version must be 3.11 or higher.

Next, while on the EPA VPN, install the packages:

```bash
python -m pip install pybmds --index-url https://gitlab.epa.gov/api/v4/projects/1508/packages/pypi/simple
python -m pip install bmds-ui --index-url https://gitlab.epa.gov/api/v4/projects/1508/packages/pypi/simple
```

This should install `pybmds` (the BMDS execution engine) and `bmds-ui` (the BMDS Desktop User Interface) along with it's related dependencies. Now, you're ready to use the application.

Start the application:

```bash
bmds-desktop
```

Check the version:

```bash
bmds-desktop --version
```

## Python Novice

Here's a simple guide for installation, where we make a few decisions for you for how to install to simplify the installation. If you run into issues, please see the FAQ below, or feel free to contact us.  The only requirement is Python on you computer.

First, install Python. Ideally, use the most recent version available; you'll need at least Python 3.11 (released in 2022). You can install it from [python.org](https://python.org), or [Anaconda](https://www.anaconda.com/), or whatever your organization prefers. During installation, ensure it is added to your path (this is an option on some Windows installers).

After installing, open a terminal on your computer. On Windows, use Windows Terminal (it should be built-in to recent Windows versions). On Mac, you can use the built-in Terminal. After starting, check if Python is available:

```bash
python -V
```

This should return a Python version (eg., `Python 3.12.3`). If Python wasn't found, follow the FAQ below to add to your path, and then continue after adding to your path.

Like many other programming languages (R, JavaScript), Python has a package installer (the <u>P</u>ackage <u>I</u>nstaller for <u>P</u>ython, or [pip](https://pypi.org/project/pip/)) for short, which makes it possible to install 3rd party packages (like `pybmds` or `bmds-desktop`). Pip is installed when you install Python automatically.

To install, run the installer command. This will install the latest version available and compatible with your computer:

```bash
python -m pip install bmds-ui
```

:::{error}
This command will not work until EPA releases the software publicly. See the [guide above](epa-guide) for EPA staff. However, we'd love feedback on the guide as these instructions will be public when we release.
:::

After install, you can start BMDS Desktop:

```bash
bmds-desktop
```

You can also check the version of BMDS Desktop installed:

```bash
bmds-desktop --version
```

### Installing multiple version of BMDS Desktop

It is possible to install multiple version of BMDS Desktop on your computer if you'd like to run different versions for different projects. For experienced Python users, you install different versions in different [virtual environments ](https://docs.python.org/3/tutorial/venv.html). For new Python users, we recommend using another package for managing environments, [pipx](https://pipx.pypa.io/). This would allow you to have multiple versions of BMDS Desktop on your computer at the same time:

First, install pipx:

```bash
python -m pip install pipx
```

Next, install a particular version of BMDS Desktop:

```bash
pipx install --suffix=-24.1 bmds-ui==24.1
pipx --list
```

After installation, you can start a particular version:

```bash
bmds-desktop-24.1 --version
bmds-desktop-24.1
```

## Familiar with Python

If you're familiar with Python, installation should be straightforward. You'll need Python 3.11 or higher. You can use any version of Python - python.org, anaconda, etc. Using [pip](https://pypi.org/project/pip/), or your preferred Python package manager, install the packages.

```bash
python -m pip install -U bmds-ui
```

We recommend using virtual environments since that allows you to have multiple versions of the software installed on your computer. The guide above describes using pipx for virtual environment management, but this is optional. You do not need to install `bmds-desktop` if you only wish to run `pybmds` without the web-based user interface.

## Frequently Asked Questions (FAQ)

One of the toughest parts of using Python is getting packages installed. When things don't follow the guide, these frequently asked questions (FAQ) may be helpful. There are also excellent resources online for installing Python package if this guide is insufficient, or contact us for help!

### How do I Add Python to my Path (Windows)?

First, you'll need to find where Python was installed on your computer. Open up Windows Explorer and search for `python.exe`. Copy the complete path of where it was installed (up to the directory, not including the executable). Next, we'll add that location to your default Path.

On the Start Menu, search for "environment" and then click the "Edit the system environment variables". This should open the "System Properties" dialog. Select "Environment Variables", and under the "User Variables" section, there should be a variable named Path. Edit the Path, and add the folder that contains the Python path for the Python installation above.

After adding the location, restart your terminal. If you type the command `python -V` the version of Python should now appear! You are ready to continue installation.

### How do I Add Python to my Path (Mac)?

First, you'll need to find where Python was installed on your computer. Open Finder, and then search for `python` on your mac; you may want to change the "Kind" to "Executable".  Next, create or edit a special file named `.zprofile` at the root of your home folder. You'll want to add a line at the bottom of the file:

```bash
export PATH="$PATH:/path/to/python/"
```

Where `/path/to/python` is the python directory which contains the Python executable.

### Can I use RStudio or R?

Yes, you should be able to; we're using standard Python packaging standards and recent versions of RStudio has nice integration with Python. Please refer to guides for how to setup Python and install a Python package using [RStudio Desktop](https://posit.co/download/rstudio-desktop/) and/or [reticulate](https://rstudio.github.io/reticulate/).

### Installing from the EPA website

The guide above uses the Python packaging index ([pypi](https://pypi.org/)) to install the BMDS packages, which automatically picks the appropriate file to install for your operating system and Python version and installs it. However, you may wish to download these packages from the U.S EPA website directly. From the U.S. EPA BMDS website go to the [downloads page](https://www.epa.gov/bmds/download-bmds). Then, you'll need to download the correct files for your operating system, as well as the Python version you have installed. For example, if you have Python 3.12 installed and are using Windows, download the following files:

- `pybmds-24.1-cp312-cp312-win_amd64.whl`
- `bmds-desktop-24.1-cp312-cp312-win_amd64.whl`

Next, start up your terminal. Ensure you have Python installed and available on your Path. It may be simplest to copy the text below and edit them to the specifics of where items are on your computer. Then, paste the commands in, one by one:

```bash
# change directory to where you downloaded the files (or move the files to your current path)
cd "path-to-downloads"
# install pipx, just for you as a user, not globally on your computer
python -m pip install --user -U pipx
# install the files, using the complete filenames for the files on your computer
pipx install --preinstall pybmds-24.1-cp312-cp312-win_amd64.whl --suffix=-24.1 bmds-ui-24.1-cp312-cp312-win_amd64.whl
# they should now appear on your terminal
pipx --list
```

### Updating BMDS Desktop

To upgrade to the latest version:

```bash
python -m pip install --upgrade bmds-ui
```

### Uninstalling BMDS Desktop

If you install BMDS Desktop using `pip`, you should be able to uninstall the same way:

```bash
python -m pip uninstall bmds-desktop
```

If you followed the guide which used `pipx` for installing multiple version of BMDS Desktop on your computer at the same time, you can uninstall using `pipx`:

```bash
# show
pipx list

# uninstall your environment with the correct prefix (see the output from pip list)
pipx uninstall bmds-ui-24.1

# you can even uninstall pipx after removing all your environments
python -m pip uninstall -U pipx
```

At this point, you can uninstall Python if you'd like.
