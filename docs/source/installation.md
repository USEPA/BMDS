# Installation

Generally, we recommend installing `bmds-ui` on your computer, which makes it possible to run both `BMDS Desktop` and `pybmds`. The guides below describe installation.

:::{note}
This guide follows pretty standard installation processes for packages in the Python ecosystem. If you have a colleague or friend that uses Python, they can likely help! If you do not, please [contact us](https://ecomments.epa.gov/bmds)!
:::

## Part 1 - install Python

First, install Python. Ideally, use the most recent version available; the minimum support version is version 3.11 (released in 2022). You can install it from [python.org](https://python.org)[^disclaimer], [Anaconda](https://www.anaconda.com/)[^disclaimer], or anywhere else you or your organization prefers. In general, using an Anaconda style distribution may make things simpler as environment management is integrated into the software (described [below](part-2)).

During installation, ensure it is added to your path (this is an option on some Windows installers).

After installing, open a terminal, and check that Python is available from your terminal:

```bash
python --version
```

:::{admonition} Which terminal should I use?
:class: info

If you installed Anaconda or one or other `conda` applications you can start a terminal from your start menu with that environment active.

If you installed Python directly and am on Windows, use Windows Terminal (it should be built-in to recent Windows versions). The terminal allows you to run different shells; use the "Command Prompt" instead of "Powershell".

If you installed Python direclty am am on a Mac, you can use the built int Terminal, or any other application you're comfortable with.
:::




(part-2)=
## Part 2 - create an environment

Python has a concept of virtual environments; is your using


::::{tab-set}

:::{tab-item} Label1
Content 1

```python
def foo():
    raise ValueError(   )
```
:::

:::{tab-item} Label2
Content 2
:::

::::

## Part 3 - install BMDS Desktop and pybmds


(internal-epa)=
### For internal EPA users

:::{tip}
This portion of the guide should be used as as we develop pre-release versions. Recommended for EPA staff. Please check the FAQ below with any questions, or reach out to [Andy Shapiro](mailto:shapiro.andy@epa.gov) for any feedback.  If anything is missing from the guide or anything is confusing, please let us know! This will be removed when we can officially release.
:::

After installing Python, open up your terminal on your computer. Ensure that Python is on your path:

```bash
python --version
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

This should return a Python version (eg., `Python 3.12.4`). If Python wasn't found, follow the FAQ below to add to your path, and then continue after adding to your path.

Like many other programming languages (R, JavaScript), Python has a package installer (the <u>P</u>ackage <u>I</u>nstaller for <u>P</u>ython, or [pip](https://pypi.org/project/pip/)) for short, which makes it possible to install 3rd party packages (like `pybmds` or `bmds-desktop`). Pip is installed when you install Python automatically.

To install, run the installer command. This will install the latest version available and compatible with your computer:

```bash
python -m pip install bmds-ui
```

:::{error}
This command will not work until EPA releases the software publicly. See the [guide above](internal-epa) for EPA staff. However, we'd love feedback on the guide as these instructions will be public when we release.
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

- [Python basics](python-basics)
- Installation Issues
    - [Adding Python to my Path (Windows)](path-windows)?
    - [Adding Python to my Path (Mac)](path-mac)?
    - [I got an error: "Syntax Error"](syntax-error)
    - [I got an error: `Running Scripts is Disabled on your System` (Windows)](psh-no-scripts)
    - [I got an error: `bmds-desktop: path not found`](path-not-found)
    - [Upgrading BMDS Desktop](upgrading-bmds-desktop)
    - [Uninstalling BMDS Desktop](uninstalling-bmds-desktop)
- Python environments
    - [Creating an environment](create-venv)
    - [Activating an environment](activate-venv)
    - [Deleting an environment](delete-venv)
- Writing pybmds code
    - [Using RStudio](using-rstudio)

One of the toughest parts of using Python is getting packages installed. When things don't follow the guide, these frequently asked questions (FAQ) may be helpful. There are also excellent resources online for installing Python package if this guide is insufficient, or contact us for help!

(syntax-error)=
#### I got an error: "Syntax Error"

(psh-no-scripts)=
#### I got an error: `Running Scripts is Disabled on your System` (Windows)

(path-not-found)=
#### I got an error: `bmds-desktop: path not found`

Variations may include: 'bmds-desktop' is not recognized as an internal or external command,
operable program or batch file.

To fix, active your Python virtual environment! See the guide below.

(python-basics)=
### Python basics

**Python** is a programming language commonly used for many different software products, from data science to websites. The `pybmds` and `bmds-ui` packages are written in the Python language and you must hav ea version of Python installed on your computer to use these packages.

A Python **package** is essentially a folder of code that you can download from the internet (or other places) that allows you to extend the Python programming language to perform a specific operation. In the context of BMDS, we have developed two packages, `pybmds`, which allows you to execute dose-response models using python, and `bmds-ui` which builds a user-interface for setting up these models. The `BMDS Desktop` application, as well as [BMDS Online](https://bmdsonline.epa.gov) use the `bmds-ui` package.

The **pip** tool which is commonly built into Python allows a user to install Python packages from the internet. Packages are most commonly stored on the [pypi.org](https://pypi.org) website.

A Python **virtual environment** makes it easy to setup multiple Python projects on your computer at the same time. For example, one day in the future, you may need multiple different versions of **pybmds** and **bmds-ui** installed if you have different projects running different versions of the software (perhaps one project was started today, one project 3 years in the past). By using a virtual environment, you can have multiple versions of these software on your computer at the same time, making it easy to switch between projects.

### Installation issues

(path-windows)=
#### Adding Python to my Path (Windows)

First, you'll need to find where Python was installed on your computer. Open up Windows Explorer and search for `python.exe`. Copy the complete path of where it was installed (up to the directory, not including the executable). Next, we'll add that location to your default Path.

On the Start Menu, search for "environment" and then click the "Edit the system environment variables". This should open the "System Properties" dialog. Select "Environment Variables", and under the "User Variables" section, there should be a variable named Path. Edit the Path, and add the folder that contains the Python path for the Python installation above.

After adding the location, restart your terminal. If you type the command `python --version` the version of Python should now appear! You are ready to continue installation.

(path-mac)=
#### Adding Python to my Path (Mac)

First, you'll need to find where Python was installed on your computer. Open Finder, and then search for `python` on your mac; you may want to change the "Kind" to "Executable".  Next, create or edit a special file named `.zprofile` at the root of your home folder. You'll want to add a line at the bottom of the file:

```bash
export PATH="$PATH:/path/to/python/"
```

Where `/path/to/python` is the python directory which contains the Python executable.

(using-pybmds)=
#### Recommended environments for using `pybmds`

If you'd prefer to not write any code, then install `bmds-ui` and start BMDS Desktop.

Otherwise, you have many options for running the `pybmds` software. The package is a standard Python package which means if you're familiar with Python, you can generally use any approach for writing code you're familiar with. However, if you're new to the Python ecosystem, we generally recommend using [jupyterlab](https://jupyter.org/)[^disclaimer], Microsoft [Visual Studio Code](https://code.visualstudio.com/)[^disclaimer], or Posit [Positron](https://github.com/posit-dev/positron)[^disclaimer] for writing and executing code.



(upgrading-bmds-desktop)=
#### Upgrading BMDS Desktop

To upgrade to the latest version:

```bash
python -m pip install --upgrade bmds-ui
```

(uninstalling-bmds-desktop)=
#### Uninstalling BMDS Desktop

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

### Python environments

(create-venv)=
#### Create an environment

(activate-venv)=
#### Activating an environment

(delete-venv)=
#### Deleting an environment

### Writing pybmds code

(using-rstudio)=
#### Using RStudio

We have had success executing the `pybmds` library using RStudio as a developer environment.  The `pybmds` and `bmds-ui` packages are installed using standard installation methods described above and RStudio works nicely with Python environments. Please refer to guides for how to setup Python and install a Python package using [RStudio Desktop](https://posit.co/download/rstudio-desktop/), and the guide above for installing `pybmds` and/or `bmds-ui` into that environment. You can also call python code directly from R using the [reticulate](https://rstudio.github.io/reticulate/) package.

[^disclaimer]: Mention of or referral to commercial products or services, and/or links to non-EPA sites does not imply official EPA endorsement of or responsibility for the opinions, ideas, data, or products presented at those locations, or guarantee the validity of the information provided.
