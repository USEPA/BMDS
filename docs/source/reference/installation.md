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

## Python Novice

Here's a simple guide for installation, where we make a few decisions for you for how to install to simplify the installation. If you run into issues, please see the FAQ below, or feel free to contact us.  The only requirement is Python on you computer.

First, install Python. Ideally, use the most recent version available; you'll need at least Python 3.11 (released in 2022). You can install it from [python.org](https://python.org), or [Anaconda](https://www.anaconda.com/), or whatever your organization prefers. During installation, ensure it is added to your path.

After installing, open a terminal on your computer. On Windows, use Windows Terminal (it should be built-in to recent Windows versions). On Mac, you can use the built-in Terminal. After starting, check if Python is available:

```bash
python -V
```

If Python wasn't found, follow the FAQ below to add to your path, and then continue after adding to your path.

Like many other programming languages (R, JavaScript), Python has a built-in package installer ([pip](https://pypi.org/project/pip/)), which makes it possible to install 3rd party packages (like `pybmds` or `bmds-desktop`). Generally, that's as easy as installing `python -m pip install <package-name>`. However, we recommend installing in another more advanced tool, [pipx](https://pipx.pypa.io/), so that you can have multiple versions of BMDS Desktop on your computer at the same time:

```bash
# install the latest version of the installer
python -m pip install -U pip
# install pipx, just for you as a user, not globally on your computer
python -m pip install --user -U pipx
# now install a particular version of pybmds and bmds desktop
pipx install --suffix=-24.1 --preinstall pybmds==24.1 bmds-desktop=24.1
# they should now appear on your terminal
pipx --list
```

After successful installation, you should be able to start the application:

```bash
bmds_desktop-24.1 --version
bmds_desktop-24.1
```

## Familiar with Python

If you're familiar with Python, installation should be straightforward. You'll need Python 3.11 or higher. You can use any version of Python - python.org, anaconda, etc. Using [pip](https://pypi.org/project/pip/), or your preferred Python package manager, install the packages.

```bash
python -m pip install -U pip
python -m pip install -U pybmds bmds-desktop
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
# install the latest version of the installer
python -m pip install -U pip
# install pipx, just for you as a user, not globally on your computer
python -m pip install --user -U pipx
# install the files, using the complete filenames for the files on your computer
pipx install --preinstall pybmds-24.1-cp312-cp312-win_amd64.whl --suffix=-24.1 bmds-desktop-24.1-cp312-cp312-win_amd64.whl
# they should now appear on your terminal
pipx --list
```

### Uninstalling BMDS Desktop

Assuming you followed the guide describe above which puts BMDS in an isolated environment using `pipx`:

```bash
# show
pipx list

# uninstall your environment with the correct prefix (see the output from pip list)
pipx uninstall bmds_desktop-24.1

# you can even uninstall pipx after removing all your environments
python -m pip uninstall -U pipx
```

At this point, you can uninstall Python if you'd like.
