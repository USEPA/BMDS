# Installation

We recommend installing `bmds-ui` on your computer; this makes it possible to run both **BMDS Desktop** for a user interface and **pybmds** for scripting.

This guide describes the basics for a user new to Python, from how to install Python, to setting up an environment, to installing the specific packages. A [simplified installer script](simplified-installer) is also available which reduces the amount of scripting required for installation.

:::{note}
This guide documents a standard installation process for packages in the Python ecosystem. If you know someone who uses Python, they can likely help! If you need support, please [contact us](https://ecomments.epa.gov/bmds).
:::

(part-1)=
## Part 1 - Install Python

If possible, use the most recent version available; the minimum support version is version 3.11 (released in 2022).

You can install Python from [python.org](https://python.org)[^disclaimer], [Anaconda](https://www.anaconda.com/)[^disclaimer], or anywhere else you or your organization prefers. In general, using an Anaconda style distribution may make things simpler as environment management is more tightly integrated into the software (described [below](part-2)).

During installation, ensure Python is added to your path (this is an option on some Windows installers).

(open-terminal)=
After installing, open a terminal.

:::{admonition} Which terminal should I use?
:class: info

On Windows:
* If you installed Anaconda or other applications which have a `conda` command (e.g., Miniconda, Miniforge, etc.), in the Start Menu, search for the Anaconda Prompt - it may be called something like "Anaconda prompt" or "Miniconda prompt".
* If you installed Python from python.org, in the Start menu search for and run the "Command Prompt". Make sure to use the "Command Prompt" instead of "Powershell".

On Mac/Linux:
* You can use the built-in Terminal, or any other application you're comfortable with.
:::

In the terminal, confirm that Python is available:

```bash
python --version
```

This should return a version, for example `Python 3.12.5`.  If you see a python version after typing this command, you're ready for the next step!  Otherwise, check the [FAQ](faq) for possible solutions.


(part-2)=
## Part 2 - Create an Environment

We recommend creating a [virtual environment](https://docs.python.org/3/tutorial/venv.html#creating-virtual-environments) for BMDS Desktop and pybmds.

Virtual environments are essentially copies of Python, but you can have multiple environments on the same computer. In addition, each environment can contain different packages with different versions of those packages.

Using virtual environments instead of installing in the "root" Python environment has several advantages:

1. With some computers and enterprise setups, it may be not be possible to install packages in the "root" Python environment without administrative rights.
2. If packages are installed, they may not easily be found by the "root" Python.
3. When using virtual environments, you can install multiple versions of the BMDS Desktop software in case you'd like to use different versions for different projects.

Python virtual environments are simple to create, but the steps are specific depending on whether you installed Python or Anaconda:

::::{tab-set}

:::{tab-item} Anaconda

There is one decision you'll need to make when creating an environment:

1. The environment name. You can use any name for the environment; in the example below we use the name `bmds-desktop`. If you plan on installing multiple versions of BMDS Desktop simultaneously, you may want to add the version to the name, for example, `bmds-desktop-24-1`. You cannot rename an environment after you create it (but you can delete or create a new one).

An example setup is below:

```bash
# create an environment named `bmds-desktop`
# with Python 3.12 in the environment
conda create --name bmds-desktop python=3.12
```

This creates an environment (by default in a path in your home directory).  Anaconda manages switching or moving among environments.

:::

:::{tab-item} Python (Windows)

There are two decisions you'll need to make when creating a virtual environment:

1. The location of the environment on your computer. You will want to put this environment in your home folder, but ideally not in folders managed by cloud-syncing software such as OneDrive or Dropbox. Environments create many small files that do not need to be backed up, and backing up will hinder the performance of your cloud sync application and may slow down your computer. In the example below, we make a new `dev` folder in your home directory, but you can use other locations.
2. The environment name. You can use any name for the environment; in the example below we use the name `bmds-desktop`. However, if you plan on installing multiple versions of BMDS Desktop simultaneously, you may want to add the version to the name, for example, `bmds-desktop-24-1`. You cannot rename an environment after you create it (but you can delete or create another one).

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

1. The location of the environment on your computer. You will want to put this environment in your home folder, but ideally not in folders managed by cloud-syncing software such as OneDrive or Dropbox. Environments create many small files that do not need to be backed up, and backing up will hinder the performance of your cloud sync application and may slow down your computer. In the example below, we make a new `dev` folder in your home directory, but you can use other locations.
2. The environment name. You can use any name for the environment; in the example below we use the name `bmds-desktop`. However, if you plan on installing multiple versions of BMDS Desktop simultaneously, you may want to add the version to the name, for example, `bmds-desktop-24-1`. You cannot rename an environment after you create it (but you can delete or create another one)

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

After creating an environment, you'll need to activate the environment. Activating the environment means that, instead of looking at the global python environment and anything that may be installed there, your system looks at the contents within the environment you specified.

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

You'll need to activate the environment **every time** you want to run BMDS Desktop from your terminal. It may be helpful to keep these instructions handy, or even create a document on your computer with specific instructions for how to access your particular setup.

:::


(part-3)=
## Part 3 - Install BMDS Desktop and pybmds

With Python installed and a virtual environment created, you're ready to install BMDS Desktop and `pybmds`.  You'll need to [activate your environment](activate-venv), and then install the packages using Python's package installer, [pip](https://pip.pypa.io/), which is included with Python.

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

To install **BMDS Desktop**, which includes `pybmds`:

```bash
python -m pip install bmds-ui
```

You can install a specific version of `bmds-ui` by specifying the version number, for example:

```bash
python -m pip install bmds-ui==24.1
```

If no version is specified, `pip` will install the latest version.

:::

After installation, you're ready to use BMDS Desktop and `pybmds`.

(desktop-shortcut)=
### The BMDS Desktop Manager

Starting the application may be difficult for users who do not frequently use the terminal. Therefore, we've created the BMDS Desktop Manager that allows you to start BMDS Desktop (or upgrade it) by double clicking an icon on your Desktop. To create the shortcut:

```bash
bmds-desktop --create-shortcut
```

This command creates the shortcut file at the location in the terminal window.

The following confirmation message is displayed:

    BMDS Desktop Manager Created:
    -----------------------------
    C:\Users\USER\dev\bmds-desktop-bmds-desktop-manager.bat

    Opening this file will start BMDS Desktop.
    You can move this file or create a shortcut to it.

    Would you like to open the folder to view "bmds-desktop-manager.bat"? (y/n)

Typing `y` and pressing Enter will open a new Finder or File Manager window displaying the folder containing the batch file.

We recommend creating a shortcut to the desktop manager file and putting the shortcut on your Desktop for easy access; you can even rename the shortcut "BMDS Desktop Manager" if you'd prefer.

To start the BMDS Desktop, double-click the BMDS Desktop Manager shortcut icon.

(part-4)=
## Part 4 - Starting BMDS Desktop

To use BMDS Desktop, [open a terminal](open-terminal) and [activate](activate-venv) an environment. Then start the application:

```bash
bmds-desktop
```

:::{tip}

If you'd rather not use the command line, consider creating a shortcut to the [BMDS Desktop Manager](desktop-shortcut) instead. After the shortcut is created, you may not need to use the terminal for future work, and you can double-click the shortcut file to start the application.

:::

The `pybmds` package is available in this environment, so you can start a Python interpreter in this environment and you can begin scripting.

### A complete example

The complete startup after installation can be summarized below, assuming you followed the guide above.


::::{tab-set}

:::{tab-item} Anaconda

Assuming you installed a conda version of Python on Windows, using the defaults in the guide above, the entire sequence would look like this:


1. In the Windows Start Menu Search field, enter "Anaconda Prompt" and start the application
2. In the Prompt, enter

```batch
conda activate bmds-desktop
```

You should see the following prompt:

```batch
(bmds-desktop) C:\Users\USERNAME>
```

3. At the prompt, enter:

```batch
bmds-desktop
```

:::

:::{tab-item} Python

Assuming you installed a python.org version of Python on Windows, using the defaults in the guide above, the entire sequence would look like this:

1. In the Windows Taskbar's Search field, enter "Command Prompt"
2. In the Command Prompt's terminal window, enter (changing `USERNAME` to your username):

```batch
cd C:\Users\USERNAME\dev
bmds-desktop\Scripts\activate
```

You should see the following prompt:

```batch
(bmds-desktop) C:\Users\USERNAME\dev>
```

3. At the prompt, enter:

```batch
bmds-desktop
```

:::

::::

(simplified-installer)=
## Simplified installer

The simplified installer (experimental) may simplify the installation process by automating the installation process and using reasonable defaults.  If you wish to have multiple versions of BMDS Desktop installed at the same time on your computer, the simplified installer script may be too simple and you may need to follow the detailed guide.

1. Install Python and open your terminal, following [Part 1](part-1) of the guide below. Install the most recent version available from [https://python.org](https://python.org).
2. Download the <a href="_static/install-bmds-desktop.py" download>installation script</a>. Move the downloaded file to the same location that is open in your terminal. Then, run the command:
    ```
    python install-bmds-desktop.py
    ```
3. The script will install BMDS Desktop, and then will create a BMDS Desktop Manager file. Create a shortcut to the [BMDS Desktop Manager](desktop-shortcut) so you can open BMDS Desktop in the future or update to a more recent version.

After running the installation script and creating a shortcut to the BMDS Desktop Manager, you shouldn't need to open your terminal in the future to update or run BMDS Desktop.

To uninstall, open your terminal and navigate to the installation script as described above, but run the command below:

```
python install-bmds-desktop.py --uninstall
```

(faq)=
## Frequently Asked Questions (FAQ)

- [Python Fundamentals](python-fundamentals)
- [Installation Issues](installation-issues)
    - [Adding Python to my Path (Windows)](path-windows)
    - [Adding Python to my Path (Mac)](path-mac)
    - [I got an error: `Syntax Error`](syntax-error)
    - [I got an error: `Running Scripts is Disabled on your System` (Windows)](psh-no-scripts)
    - [I got an error: `bmds-desktop: path not found`](path-not-found)
- [Upgrading BMDS Desktop](upgrading-bmds-desktop)
- [Installing Multiple Versions of BMDS Desktop](multiple-versions)
- [Uninstalling BMDS Desktop](uninstalling-bmds-desktop)
- [Writing `pybmds` Code](writing-pybmds-code)
- [Writing `pybmds` Code in R](using-r)

(python-fundamentals)=
### Python Fundamentals

**Python** is a programming language commonly used for many different software products, from data science to websites. The `pybmds` and `bmds-ui` packages are written in the Python language and you must have a version of Python installed on your computer to use these packages.

A Python **package** is essentially a folder of code you download that enables you to extend the Python programming language to perform a specific predefined operations in that package. In the context of BMDS, we have developed two packages:

* `pybmds`, which allows you to execute dose-response models in Python scripts, and
* `bmds-ui` which builds a user-interface for dose-response modeling. The BMDS Desktop application and [BMDS Online](https://bmdsonline.epa.gov) use the same package, configured differently.

The **pip** tool which is built into Python, allows a user to install Python packages from the internet. Packages are available to all Python users on [pypi.org](https://pypi.org) (this is Python's equivalent to R's [CRAN](https://cran.r-project.org/)).

A Python **virtual environment** makes it easy to set up multiple Python projects on your computer at the same time, making it easy to switch between projects. Environments also make it easier for your computer to find Python packages associated with a particular application.

(installation-issues)=
### Installation issues

(path-windows)=
#### Adding Python to my Path (Windows)

Variations of this error may include: `python is not recognized as an internal or external command, operable program or batch file.`

To troubleshoot this error (instructions assume Windows 11):

1. Locate where Python was installed on your computer. Open Windows Explorer and search for `python.exe`.
2. Copy the complete path of where `python.exe` is found (up to the directory only; do not include the filename). You will add this folder location to your default Path.
3.	In the Start Menu, search for "environment" and then click the "Edit environment variables for your account" item. The system should open the "Environment Variables" dialog. This will enable you to update the Path variable *without* needing administrative rights.
    ```{figure} _static/img/windows-edit-enviro-1.jpg
    :alt: Screenshot of how to find the "Edit environment variables for your account" icon
    ```
4.	Under the "User variables for USERNAME" section, locate the variable named Path.
    ```{figure} _static/img/windows-edit-enviro-2.jpg
    :alt: Screenshot of User Environment Variables
    ```
5.	Select the Path variable and then the Edit button to display the "Edit environment variable" dialog.
6.	In the "Edit environment variable" dialog, select New to create a new line in the display. Click inside the new line to select it and paste the folder path you copied from Step 2 above.
    ```{figure} _static/img/windows-edit-enviro-3.jpg
    :alt: Screenshot of adding a new directory to the Path User environment variable
    ```
7.	Shut down the Command Prompt window and restart it.
8.	At the reopened Command Prompt window, enter the command `python --version` and the Python version number should now appear! Continue with the installation.


(path-mac)=
#### Adding Python to my Path (Mac)

1.	To find where Python was installed on your macOS, open Finder, and search for "python". In the search bar, you may want to change the "Kind" to "Executable" to limit results.
2.	With the file displayed in the search results, highlight the file and press Command + I to open the Get Info window.
3.	In the Get Info window, triple-click the file path beside "Where" to select it. Then press Command (⌘)+C to copy the file or folder path.
4.	In a new Finder window, navigate to the root of your home folder under Users. You will need to add a path to a special file named ".zprofile". So-called *dot-files* are typically hidden in the standard Finder view. Press the keyboard shortcut Shift + Command + . (period) to toggle the display of any dot-files in the folder.
5.	Use TextEdit to open the .zprofile file if it exists or use TextEdit to create a new text file and save it as .zprofile in the root of your home folder.
6.	In either case, add the following line to the bottom of the .zprofile file:

    ```bash
    export PATH="$PATH:/path/to/python/"
    ```

    where /path/to/python is the python directory containing the Python executable. You can paste the path you copied in Step 3.

7.	Restart your terminal and enter the command `python --version` and the Python version number should now appear! Continue with the installation.


(syntax-error)=
#### I got an error: `Syntax Error`

You may see a Syntax Error message if you've started the Python interpreter and then typed installation commands not written in the Python language.

From within the Python interpreter, you can write Python scripts using `pybmds`. However, it's it's not where you want to be for the installation guide. For installing BMDS Desktop, you need to be in a terminal in a location where Python can be executed. The good news is, if you're inside a Python interpreter, you're in the correct place for installation, but you'll need to exit the interpreter.

To exit the interactive Python interpreter and return to your terminal, type `exit()` and press Enter.

(psh-no-scripts)=
#### I got an error: `Running Scripts is Disabled on your System` (Windows)

Windows has two built-in command line shells - Command Prompt and PowerShell. We recommend using the Command Prompt, not PowerShell.

If you see this error, you're likely running PowerShell.

Follow the guide on [opening terminals](open-terminal) above.

(path-not-found)=
#### I got an error: `bmds-desktop: path not found`

Variations may include: `'bmds-desktop' is not recognized as an internal or external command,
operable program or batch file`.

To fix, [activate](activate-venv) your Python virtual environment after successfully installing BMDS Desktop.

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

It is possible to install multiple versions of BMDS Desktop on your computer. You can even run different version of BMDS Desktop for different projects. For example, you may start an assessment using one version of BMDS Desktop and continue using that version for consistency, but may want to use a newer version for a different assessment. You can install the same version of BMDS Desktop in multiple separate environments, with each dedicated to a specific assessment project. That way, all datasets and results are clustered in their own isolated, reproducible environment.

Follow the [environment guide](part-2) above and create a new environment for each version.  Then, activate the environment.

You cannot install multiple versions of BMDS Desktop in the same environment.

(uninstalling-bmds-desktop)=
### Uninstalling BMDS Desktop

After [activating your environment](activate-venv), you can uninstall BMDS Desktop with the following command:

```bash
python -m pip uninstall bmds-desktop
```

You can also remove the environment:

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

You have many options for running the `pybmds` software. The package is a standard Python package which means if you're familiar with Python, you can generally use any approach for writing code you're familiar with.

However, if you're new to the Python ecosystem, we generally recommend using [jupyterlab](https://jupyter.org/)[^disclaimer], Microsoft [Visual Studio Code](https://code.visualstudio.com/)[^disclaimer], or Posit [Positron](https://github.com/posit-dev/positron)[^disclaimer] for writing and executing code.

If you'd prefer to not write any code, then install `bmds-ui` and start BMDS Desktop.

(using-r)=
### Writing `pybmds` Code in R

You can use use RStudio to execute `pybmds` using the [reticulate](https://rstudio.github.io/reticulate/) package; you'll still need to follow the guide above to install the software on your computer in a Python environment.

See the [using R](recipes/using-r.md) recipe for an example.

[^disclaimer]: Mention of or referral to commercial products or services, and/or links to non-EPA sites does not imply official EPA endorsement of or responsibility for the opinions, ideas, data, or products presented at those locations, or guarantee the validity of the information provided.
