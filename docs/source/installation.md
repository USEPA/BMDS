# Installation

We recommend installing `bmds-ui` on your computer; this makes it possible to run both **BMDS Desktop** for a user interface and **pybmds** for scripting. This guide describes the basics from how to install python, to setting up an environment, to installing the specific packages.


:::{note}
This guide documents a standard installation processes for packages in the Python ecosystem. If you have a colleague or friend that uses Python, they can likely help! If you do not and need support, please [contact us](https://ecomments.epa.gov/bmds).
:::

(part-1)=
## Part 1 - Install Python

First, install Python. Ideally, use the most recent version available; the minimum support version is version 3.11 (released in 2022). You can install it from [python.org](https://python.org)[^disclaimer], [Anaconda](https://www.anaconda.com/)[^disclaimer], or anywhere else you or your organization prefers. In general, using an Anaconda style distribution may make things simpler as environment management is more tightly integrated into the software (described [below](part-2)).

During installation, ensure Python is added to your path (this is an option on some Windows installers).

(open-terminal)=
After installing, open a terminal.

:::{admonition} Which terminal should I use?
:class: info

If you installed Anaconda or one or other `conda` applications, you can start a terminal from your start menu with that environment active.

If you installed Python directly on Windows, use Windows Terminal (it should be built-in to recent Windows versions). The terminal allows you to run different shells; use the "Command Prompt" instead of "Powershell".

If you installed Python directly on a Mac, you can use the built-in Terminal, or any other application you're comfortable with.
:::

In the terminal, confirm that Python is available:

```bash
python --version
```

This should return a version, for example `Python 3.12.4`.  If you see a python version after typing this command, you're ready for the next step!  Otherwise, check the [FAQ](faq) for possible solutions.


(part-2)=
## Part 2 - Create an Environment

With Python installed and available in your terminal, you can now install BMDS Desktop and `pybmds`. However, we recommend creating a [virtual environment](https://docs.python.org/3/tutorial/venv.html#creating-virtual-environments), and installing inside the environment. Virtual environments are essentially copies of Python, but each one can have different packages with different versions, and you can have multiple environments on the same computer.

Using virtual environments instead of installing in the "root" Python environment is advantageous for a few reasons:

1. With some computers and enterprise setups, it may be not be possible to install packages in the "root" Python environment without administrative rights
2. If packages are installed, they may not easily be found by the "root" Python
3. When using virtual environments, you can install multiple versions of the BMDS Desktop software in case you'd like to use different versions for different projects

Creating python virtual environments are simple, but depend on if you installed Python or Anaconda:

::::{tab-set}

:::{tab-item} Anaconda

There is one decision you'll need to make when creating an environment:

1. The environment name. You can use any name for the environment; in the example below we use the name `bmds-desktop`. If you plan on installing multiple versions of BMDS Desktop simultaneously, you may want to add the version to the name, for example, `bmds-desktop-24-1`. You cannot rename an environment after you create it (but you can delete or make another one).

An example setup is below:

```bash
# create an environment named `bmds-desktop`
# with Python 3.12 in the environment
conda create --name bmds-desktop python=3.12
```

This creates an environment (by default in a path in your home directory).  Anaconda handles finding environments on your filesystem for you (this simplifies the process versus Python below).

:::

:::{tab-item} Python (Windows)

There are two decisions you'll need to make when creating a virtual environment:

1. The location of the environment on your computer. You will want to put this environment in a place you in your home folder, but ideally a location that is not set up with cloud syncing software such as OneDrive. In the example below, we make a new `dev` folder in your home directory, but you use other locations.
2. The environment name. You can use any name for the environment; in the example below we use the name `bmds-desktop`. If you plan on installing multiple versions of BMDS Desktop simultaneously, you may want to add the version to the name, for example, `bmds-desktop-24-1`. You cannot rename an environment after you create it (but you can delete or make another one).

An example setup is below. First, navigate to the folder we want to create the environment in (where `USERNAME` is your username):

```batch
cd C:\Users\USERNAME
mkdir dev
cd dev
```

And then create the environment:

```batch
python -m venv bmds-desktop
```

The instructions above create a virtual environment here: `C:\Users\USERNAME\dev\bmds-desktop`.

:::

:::{tab-item} Python (Mac/Linux)

There are two decisions you'll need to make when creating a virtual environment:

1. The location of the environment on your computer. You will want to put this environment in a place you in your home folder, but ideally a location that is not set up with cloud syncing software such as OneDrive. In the example below, we make a new `dev` folder in your home directory, but you use other locations.
2. The environment name. You can use any name for the environment; in the example below we use the name `bmds-desktop`. If you plan on installing multiple versions of BMDS Desktop simultaneously, you may want to add the version to the name, for example, `bmds-desktop-24-1`. You cannot rename an environment after you create it (but you can delete or make another one).

An example setup is below:

```bash
# create a dev folder in your home directory
cd ~
mkdir dev
cd dev

# create a virtual environment named `bmds-desktop`
python -m venv bmds-desktop
```

The instructions above create a virtual environment here: `~/dev/bmds-desktop`.

:::

::::

**Additional References:**

- [Anaconda Environments](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)
- [Python Virtual Environments](https://packaging.python.org/en/latest/tutorials/installing-packages/#creating-virtual-environments)

(activate-venv)=
### Activating an environment

After creating an environment, you'll need to activate the environment. Activating the environment means that instead of looking at the global python environment, you look at the contents within the environment.

::::{tab-set}

:::{tab-item} Anaconda

You'll need to remember the environment name you created and then activate it:

```bash
conda activate bmds-desktop
```

If successful, you should see the environment name in front of your terminal in parenthesis, e.g.:

```bash
(bmds-desktop) $
```

:::


:::{tab-item} Python (Windows)

You'll need change directories to where you created the virtual environment, and then activate it:

```batch
cd C:\Users\USERNAME\dev

bmds-desktop\Scripts\activate
```

If successful, you should see the environment name in front of your terminal in parenthesis, e.g.:

```batch
(bmds-desktop) C:\Users\USERNAME\dev>
```

:::

:::{tab-item} Python (Mac/Linux)

You'll need change directories to where you created the virtual environment, and then activate it:

```bash
# change directories to the location you created the environment
cd ~/dev
# activate the environment
bmds-desktop/bin/activate
```

If successful, you should see the environment name in front of your terminal in parenthesis, e.g.:

```bash
(bmds-desktop) $
```

:::

::::

:::{tip}

You'll need to activate the environment **every time** you want to run BMDS Desktop. It may be helpful to keep these instructions handy, or even create a document on your computer with specific instructions for how to access with your particular setup.

:::


(part-3)=
## Part 3 - Install BMDS Desktop and pybmds

With Python installed and a virtual environment created, we're ready to install BMDS Desktop and `pybmds`.  You'll need to activate [your environment](activate-venv), and then install the packages using Python's package installer, [pip](https://pip.pypa.io/), which is included with Python.

(internal-epa)=
:::{admonition} NOTE - for internal EPA testing
:class: tip

This portion of the guide should be used as as we develop pre-release versions. Recommended for EPA staff. Please check the FAQ below with any questions, or reach out to [Andy Shapiro](mailto:shapiro.andy@epa.gov) for any feedback.  If anything is missing from the guide or anything is confusing, please let us know! This will be removed when we can officially release.

---

While connected to the EPA VPN install the packages:

```bash
python -m pip install pybmds --index-url https://gitlab.epa.gov/api/v4/projects/1508/packages/pypi/simple
python -m pip install bmds-ui --index-url https://gitlab.epa.gov/api/v4/projects/1508/packages/pypi/simple
```

:::

:::{admonition} NOTE - this will not work until we officially release!
:class: caution

The content below below will not work until we are cleared for public release; please follow the [EPA installation guide](internal-epa). This message will be removed prior to official release.

---

To install **BMDS Desktop** (which includes `pybmds`):

```bash
python -m pip install bmds-ui
```

You can install a specific version of `bmds-ui` by specifying the version number, for example:

```bash
python -m pip install bmds-ui==24.1
```

If no version is specified, it will install the latest.

:::

After installation, we're ready to use.


(part-4)=
## Part 4 - Using BMDS Desktop and pybmds

To use BMDS Desktop, [open a terminal](open-terminal) and then [activate](activate-venv) an environment. Then, you can start the application the application:

```bash
bmds-desktop
```

If this works, you're good to go! However, if you'd rather not use the command line, now might be a good time to create a [desktop shortcut](desktop-shortcut) instead. After the shortcut is created, you may not need to use the terminal for future work.

(desktop-shortcut)=
### The BMDS Desktop Manager

Starting the application may be difficult for users who do not frequently use the terminal. Therefore, we've created the BMDS Desktop Manager that allows you to start BMDS Desktop (or upgrade it) by double clicking an icon on your Desktop. To create the shortcut:

```bash
bmds-desktop --create-shortcut
```

This creates the file wherever you currently are located in your terminal. You can move this file where you'd like using Windows Explorer or Mac Finder. We recommend creating a shortcut and putting the shortcut on your Desktop for easy access; you can even rename the shortcut "BMDS Desktop Manager" if you'd prefer.

To start the BMDS Desktop, double-click the BMDS Desktop Manager icon.

(faq)=
## Frequently Asked Questions (FAQ)

- [Python Fundamentals](python-fundamentals)
- [Installation Issues](installation-issues)
    - [Adding Python to my Path (Windows)](path-windows)
    - [Adding Python to my Path (Mac)](path-mac)
    - [I got an error: "Syntax Error"](syntax-error)
    - [I got an error: `Running Scripts is Disabled on your System` (Windows)](psh-no-scripts)
    - [I got an error: `bmds-desktop: path not found`](path-not-found)
- [Upgrading BMDS Desktop](upgrading-bmds-desktop)
- [Installing Multiple Versions of BMDS Desktop](multiple-versions)
- [Uninstalling BMDS Desktop](uninstalling-bmds-desktop)
- [Writing `pybmds` Code](writing-pybmds-code)
- [Writing `pybmds` Code in R](using-r)

(python-fundamentals)=
### Python Fundamentals

**Python** is a programming language commonly used for many different software products, from data science to websites. The `pybmds` and `bmds-ui` packages are written in the Python language and you must hav ea version of Python installed on your computer to use these packages.

A Python **package** is essentially a folder of code that you can download that allows you to extend the Python programming language to perform a specific operation. In the context of BMDS, we have developed two packages, `pybmds`, which allows you to execute dose-response models using python, and `bmds-ui` which builds a user-interface for dose-response modeling. The `BMDS Desktop` application, as well as [BMDS Online](https://bmdsonline.epa.gov) use the `bmds-ui` package.

The **pip** tool which is commonly built into Python allows a user to install Python packages from the internet. Packages are most commonly stored on the [pypi.org](https://pypi.org) website.

A Python **virtual environment** makes it easy to setup multiple Python projects on your computer at the same time. By using a virtual environment, you can have multiple versions of these software on your computer at the same time, making it easy to switch between projects. Environments also make it easier for your computer to find Python packages associated with a particular application.

(installation-issues)=
### Installation issues

(path-windows)=
#### Adding Python to my Path (Windows)

Variations of this error may include: `python` is not recognized as an internal or external command, operable program or batch file.

First, you'll need to find where Python was installed on your computer. Open up Windows Explorer and search for `python.exe`. Copy the complete path of where it was installed (up to the directory, not including the executable). Next, we'll add that location to your default Path.

On the Start Menu, search for "environment" and then click the "Edit the system environment variables". This should open the "System Properties" dialog. Select "Environment Variables", and under the "User Variables" section, there should be a variable named Path. Edit the Path, and add the folder that contains the Python path for the Python installation above.  You should be able to update "User Variables" without administrative rights.

After adding the location, restart your terminal. If you type the command `python --version` the version of Python should now appear! You are ready to continue installation.

(path-mac)=
#### Adding Python to my Path (Mac)

First, you'll need to find where Python was installed on your computer. Open Finder, and then search for `python` on your mac; you may want to change the "Kind" to "Executable".  Next, create or edit a special file named `.zprofile` at the root of your home folder. You'll want to add a line at the bottom of the file:

```bash
export PATH="$PATH:/path/to/python/"
```

Where `/path/to/python` is the python directory which contains the Python executable.

(syntax-error)=
#### I got an error: "Syntax Error"

You may see a Syntax error message if you've started a session of Python to interactively type commands. From this view, you could can write Python scripts and execute functions in `pybmds` or the standard library, but it's not where you'll want to for the installation guide; you'll need to be in a terminal shell in a location where python can be executed.

To exit and return to your terminal, type `exit()` and press return.

(psh-no-scripts)=
#### I got an error: `Running Scripts is Disabled on your System` (Windows)

Windows has two built-in command line shells - Command Prompt (cmd) and PowerShell. We recommend using the Command Prompt for installation, not PowerShell.  If you see this error, you're running PowerShell. Follow the guide on [opening terminals](open-terminal) above.

(path-not-found)=
#### I got an error: `bmds-desktop: path not found`

Variations may include: 'bmds-desktop' is not recognized as an internal or external command,
operable program or batch file.

To fix, [activate](activate-venv) your Python virtual environment after successfully install BMDS Desktop.

(upgrading-bmds-desktop)=
### Upgrading BMDS Desktop

After [activating your environment](activate-venv), you can upgrade to the latest version:

:::{admonition} NOTE - this will not work until we officially release!
:class: caution

The content below below will not work until we are cleared for public release; please follow the [EPA installation guide](internal-epa). This message will be removed prior to official release.

---

```bash
python -m pip install --upgrade bmds-ui
```
:::

:::{admonition} NOTE - for internal EPA testing
:class: tip

Use this guide instead until we are cleared for release.

---

While connected to the EPA VPN, install the packages:

```bash
python -m pip install --upgrade bmds-ui --index-url https://gitlab.epa.gov/api/v4/projects/1508/packages/pypi/simple
```

:::

(multiple-versions)=
### Installing Multiple Versions of BMDS Desktop

It is possible to install multiple version of BMDS Desktop on your computer if you'd like to run different versions for different projects. Follow the [environment guide](part-2) above, and create a new environment for each version.  Then, active the environment. You cannot have multiple versions installed in the same environment.

(uninstalling-bmds-desktop)=
### Uninstalling BMDS Desktop

After [activating your environment](activate-venv), you can uninstall:

```bash
python -m pip uninstall bmds-desktop
```

You can remove the environment if you'd like as well.

::::{tab-set}

:::{tab-item} Anaconda

```bash
conda remove -n bmds-desktop --all
```

Where `bmds-desktop` is the name of your environment.

:::

:::{tab-item} Python

Open the path to your virtual environment using Windows Explorer and delete the folder.

:::

::::

(writing-pybmds-code)=
### Writing pybmds Code

You have many options for running the `pybmds` software. The package is a standard Python package which means if you're familiar with Python, you can generally use any approach for writing code you're familiar with. However, if you're new to the Python ecosystem, we generally recommend using [jupyterlab](https://jupyter.org/)[^disclaimer], Microsoft [Visual Studio Code](https://code.visualstudio.com/)[^disclaimer], or Posit [Positron](https://github.com/posit-dev/positron)[^disclaimer] for writing and executing code.

If you'd prefer to not write any code, then install `bmds-ui` and start BMDS Desktop.

(using-r)=
### Writing `pybmds` Code in R

You can use use RStudio to execute `pybmds` using the [reticulate](https://rstudio.github.io/reticulate/) package; you'll still need to follow the guide above to install the software on your computer in a Python environment. See the [using R](recipes/using-r.md) recipe for an example.

[^disclaimer]: Mention of or referral to commercial products or services, and/or links to non-EPA sites does not imply official EPA endorsement of or responsibility for the opinions, ideas, data, or products presented at those locations, or guarantee the validity of the information provided.
