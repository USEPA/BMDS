# Installation

We recommend installing `bmds-ui` on your computer; this makes it possible to run both **BMDS Desktop** for a user interface and **pybmds** for scripting.

This guide describes the basics for a user new to Python, from how to install Python, to setting up an environment, to installing the specific packages.

:::{note}
This guide documents a standard installation process for packages in the Python ecosystem. If you know someone who uses Python, they can likely help! If you need support, please go to the [BMDS web site](https://www.epa.gov/bmds) and select the Contact Us link.
:::

(part-1)=
## Part 1 - Install Python

You can install Python from [python.org](https://www.python.org/downloads/)[^disclaimer], [Anaconda](https://www.anaconda.com/download)[^disclaimer], or anywhere else you or your organization prefers. In general, using an Anaconda style distribution (e.g., Anaconda, Miniconda, Miniforge, etc.) may make installation simpler as environment management is more tightly integrated into the software (described [below](part-2)).

* If installing Python from [python.org](https://www.python.org/downloads/)[^disclaimer], make sure that Python is added to your path (this is an option on Windows installers). If possible, use the most recent version available; the minimum supported version is 3.11.0 (released in 2022).
* If installing a Conda distribution, make sure it is a relatively recent version (released within the last 3 years).

(open-terminal)=
After installing, open a terminal to check that either `python` or `conda` is available from your terminal.  Select the appropriate guide below, depending on your method of installing.

::::{tab-set}
:::{tab-item} Anaconda

If you installed Anaconda style distribution (e.g., Anaconda, Miniconda, Miniforge, etc.), in the Start Menu, search for the Anaconda Prompt - it may be called something like "Anaconda Prompt", "Miniconda Prompt", or "Miniforge Prompt". Select this item to open the terminal.

In the terminal, confirm that `conda` is available:

```bash
conda --version
```

This should return a conda version, for example `conda 24.7.1`. Conda can install different versions of Python, and you should be able to install a version that is compatible with `pybmds` and BMDS Desktop as long as the conda installation is not too old.

If you see a conda version after typing this command, you’re ready for the next step! Otherwise, check the [FAQ](faq) for possible solutions.

:::
:::{tab-item} Python

**On Windows:**

* In the Start menu, search for and run the "Command Prompt". Make sure to use the "Command Prompt" instead of "Powershell".

**On Mac/Linux:**

* You can use the built-in Terminal, or any other application you're comfortable with.

In the terminal, confirm that `python` is available:

```bash
python --version
```

This should return a Python version, for example `Python 3.13.0`. The minimum supported version of Python for `pybmds` and BMDS Desktop is Python 3.11.0; anything more recent than this also should work.

If you see a Python version after typing this command, you’re ready for the next step! Otherwise, check the [FAQ](faq) for possible solutions.

:::
::::

(part-2)=
## Part 2 - Create an Environment

We recommend creating a [virtual environment](https://docs.python.org/3/tutorial/venv.html#creating-virtual-environments) for BMDS Desktop and pybmds.

Virtual environments are essentially copies of Python, and you can have multiple environments on the same computer. In addition, each environment can contain different packages with different versions of those packages.

Using virtual environments instead of installing in the "root" Python environment has several advantages:

1. With some computers and enterprise setups, it may not be possible to install packages in the "root" Python environment without administrative rights.
2. If packages are installed, they may not easily be found by the "root" Python.
3. When using virtual environments, you can install multiple versions of the BMDS Desktop software in case you'd like to use different versions for different projects.

Python virtual environments are simple to create, but the steps are specific depending on whether you installed Python or Anaconda:

::::{tab-set}

:::{tab-item} Anaconda

There is one decision you'll need to make when creating an environment: **the environment name**.

You can use any name for the environment; in the example below we use the name `bmds-desktop`. If you plan on installing multiple versions of BMDS Desktop simultaneously, you may want to add the version to the name, for example, `bmds-desktop-24-1`. You cannot rename an environment after you create it (but you can delete or create a new one).

An example setup is below:

```bash
# create an environment named `bmds-desktop`
# with Python 3.13 in the environment
conda create --name bmds-desktop python=3.13
```

This creates an environment (by default in a path in your home directory).  Anaconda manages switching or moving among environments.

:::

:::{tab-item} Python (Windows)

There are two decisions you'll need to make when creating a virtual environment:

1. **The location of the environment on your computer.** You will want to put this environment in your home folder, but ideally not in folders managed by cloud-syncing software such as OneDrive or Dropbox. Environments create many small files that do not need to be backed up; backing them up may hinder the performance of your cloud sync application and may slow down your computer. In the example below, we make a new `dev` folder in the home directory, but you can use other locations.
2. **The environment name.** You can use any name for the environment; in the example below we use the name `bmds-desktop`. However, if you plan on installing multiple versions of BMDS Desktop simultaneously, you may want to add the version to the name, for example, `bmds-desktop-24-1`. You cannot rename an environment after you create it (but you can delete or create another one).

An example setup is below. First, navigate to the folder where we want to create the environment (where `USERNAME` is your username):

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

1. **The location of the environment on your computer.** You will want to put this environment in your home folder, but ideally not in folders managed by cloud-syncing software such as OneDrive or Dropbox. Environments create many small files that do not need to be backed up; backing them up may hinder the performance of your cloud sync application and may slow down your computer. In the example below, we make a new `dev` folder in the home directory, but you can use other locations.
2. **The environment name.** You can use any name for the environment; in the example below we use the name `bmds-desktop`. However, if you plan on installing multiple versions of BMDS Desktop simultaneously, you may want to add the version to the name, for example, `bmds-desktop-24-1`. You cannot rename an environment after you create it (but you can delete or create another one)

An example setup is below:

```bash
cd ~  # navigate to your home directory
mkdir dev
cd dev
```

And then create the environment:

```bash
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

After creating an environment, you'll need to activate the environment. Activating the environment means that, instead of looking at the global Python environment and anything that may be installed there, your system looks at the contents within the environment you specified.

::::{tab-set}

:::{tab-item} Anaconda

You'll need to remember the environment name you created and then activate it:

```bash
conda activate bmds-desktop
```

If successful, you should see the environment name in front of your terminal prompt in parenthesis, e.g.:

```bash
(bmds-desktop) $
```

:::


:::{tab-item} Python (Windows)

Change the directory to where you created the virtual environment, and then activate it:

```batch
cd C:\Users\USERNAME\dev

bmds-desktop\Scripts\activate
```

If successful, you should see the environment name in front of your terminal prompt in parenthesis, e.g.:

```batch
(bmds-desktop) C:\Users\USERNAME\dev>
```

:::

:::{tab-item} Python (Mac/Linux)

Change the directory to where you created the virtual environment, and then activate it:

```bash
cd ~/dev

bmds-desktop/bin/activate
```

If successful, you should see the environment name in front of your terminal prompt in parenthesis, e.g.:

```bash
(bmds-desktop) ~/dev/bmds-desktop $
```

:::

::::

:::{tip}

You'll need to activate the environment **every time** you want to run BMDS Desktop from your terminal. It may be helpful to keep these instructions handy, or even create a document on your computer with specific instructions for how to access your particular setup.

:::


(part-3)=
## Part 3 - Install BMDS Desktop and pybmds

With Python installed and a virtual environment created, you're ready to install BMDS Desktop and `pybmds`.  You'll need to [activate your environment](activate-venv), and then install the packages using Python's package installer, [pip](https://pip.pypa.io/), which is included with Python.

To install **BMDS Desktop**, which includes `pybmds`:

```bash
python -m pip install bmds-ui
```

You can install a specific version of `bmds-ui` by specifying the version number, for example:

```bash
python -m pip install bmds-ui==24.1
```

If no version is specified, `pip` will install the latest version.

After installation, you're ready to use BMDS Desktop and `pybmds`.

(desktop-shortcut)=
### Create the BMDS Desktop Manager Shortcut

Starting BMDS Desktop may be difficult for users who do not frequently use the terminal. Therefore, we've created the BMDS Desktop Manager that enables you to start or upgrade BMDS Desktop by double clicking an icon on your Desktop.

Follow the prior steps to install BMDS Desktop. With the environment active, create the shortcut:

```bash
bmds-desktop --create-shortcut
```

This command creates the shortcut file at the location in the terminal window.

The following confirmation message is displayed:

    BMDS Desktop Manager Created:
    -----------------------------
    C:\Users\USERNAME\dev\bmds-desktop-bmds-desktop-manager.bat

    Opening this file will start BMDS Desktop.
    You can move this file or create a shortcut to it.

    Would you like to open the folder to view "bmds-desktop-manager.bat"? (y/n)

Typing `y` and pressing Enter will open a new Finder or File Manager window displaying the folder containing the batch file.

We recommend creating a shortcut to the desktop manager file and putting the shortcut on your Desktop for easy access; you can even rename the shortcut "BMDS Desktop Manager" if you'd prefer.

To start the BMDS Desktop, double-click the BMDS Desktop Manager shortcut icon.

(part-4)=
## Part 4 - Starting BMDS Desktop

It is recommended to start BMDS Desktop using the BMDS Desktop Shortcut Manager, described in the prior [section](desktop-shortcut).

Alternatively, you can start BMDS Desktop from your terminal. To use BMDS Desktop, [open a terminal](open-terminal), and [activate](activate-venv) an environment. Then start the application:

```bash
bmds-desktop
```

:::{tip}

If you'd rather not use the command line, consider creating a shortcut to the [BMDS Desktop Manager](desktop-shortcut) instead. After the shortcut is created, you may not need to use the terminal for future work, and you can double-click the shortcut file to start the application.

:::

The `pybmds` package is available in this environment, so you can start a Python interpreter in this environment and you can begin scripting.

### A complete example

The complete startup after installation on Windows can be summarized below, assuming you followed the guide above.

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

A Python **package** is essentially a folder of code you download that enables you to extend the Python programming language to perform specific predefined operations in that package. In the context of BMDS, we have developed two packages:

* `pybmds`, which allows you to execute dose-response models in Python scripts, and
* `bmds-ui`, which builds a user-interface for dose-response modeling. The BMDS Desktop application and [BMDS Online](https://bmdsonline.epa.gov) use the same package, but configured differently.

The **pip** tool built into Python, allows a user to install Python packages from the internet. Packages are available to all Python users on [pypi.org](https://pypi.org) (this is Python's equivalent to R's [CRAN](https://cran.r-project.org/)).

A Python **virtual environment** makes it easy to set up multiple Python projects on your computer at the same time, making it easy to switch between projects. Environments also make it easier for your computer to find Python packages associated with a particular application.

(installation-issues)=
### Installation issues

(path-windows)=
#### Adding Python to my Path (Windows)

Variations of this error may include: `python is not recognized as an internal or external command, operable program or batch file.`

To troubleshoot this error (instructions assume Windows 11):

1. Locate where Python was installed on your computer. Open Windows Explorer and search for `python.exe`.
2. Copy the complete path of where `python.exe` is found (up to the directory only; do not include the filename). You will add this folder location to your default Path.
3.	In the Start Menu, search for "environment" and then click the "Edit environment variables for your account" item (as shown in the following screening). The system should open the "Environment Variables" dialog. This will enable you to update the Path variable *without* needing administrative rights.
    ```{figure} _static/img/windows-edit-enviro-1.jpg
    :alt: Screenshot of how to find the "Edit environment variables for your account" icon
    ```
4.	Under the "User variables for USERNAME" section, locate the variable named Path.
    ```{figure} _static/img/windows-edit-enviro-2.jpg
    :alt: Screenshot of User Environment Variables
    ```
5.	Select the Path variable and then the Edit button to display the "Edit environment variable" dialog.
6.	In the "Edit environment variable" dialog (shown in the following screenshot), select New to create a new line in the display. Click inside the new line to select it and paste the folder path you copied from Step 2 above. Select OK to save your changes and close the dialog.
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
    where `/path/to/python` is the directory containing the Python executable. You can paste the path you copied in Step 3.

7.	Restart your terminal and enter the command `python --version` and the Python version number should now appear! Continue with the installation.


(syntax-error)=
#### I got an error: `Syntax Error`

You may see a Syntax Error message if you've started the Python interpreter and then typed installation commands not written in the Python language.

From within the Python interpreter, you can write Python scripts using `pybmds`. However, the Python interpreter is not where you want to be for the installation guide. For installing BMDS Desktop, you need to be in a terminal in a location where Python can be executed. The good news is, if you're inside a Python interpreter, you're in the correct place for installation, but you'll need to exit the interpreter.

To exit the interactive Python interpreter and return to your terminal, type `exit()` and press Enter.

(psh-no-scripts)=
#### I got an error: `Running Scripts is Disabled on your System` (Windows)

Windows has two built-in command line shells - Command Prompt and PowerShell. We recommend using the Command Prompt, not PowerShell.

If you see this error, you're likely running PowerShell.

Follow the guide on [opening terminals](open-terminal) above.

(path-not-found)=
#### I got an error: `bmds-desktop: path not found`

Variations on this error may include: `'bmds-desktop' is not recognized as an internal or external command,
operable program or batch file`.

To fix, [activate](activate-venv) your Python virtual environment after successfully installing BMDS Desktop.

(upgrading-bmds-desktop)=
### Upgrading BMDS Desktop

After [activating your environment](activate-venv), you can upgrade to the latest version:

```bash
python -m pip install --upgrade bmds-ui
```

(multiple-versions)=
### Installing Multiple Versions of BMDS Desktop

It is possible to install multiple versions of BMDS Desktop on your computer. You can even run different versions of BMDS Desktop for different projects.

For example, you may start an assessment using one version of BMDS Desktop and continue using that version for consistency, but may want to use a newer version for a different assessment. Or, you can install the same version of BMDS Desktop in multiple separate environments, with each dedicated to a specific assessment project. That way, all datasets and results are clustered in their own isolated, reproducible environment.

Follow the [environment guide](part-2) above and create a new environment for each version.  Then, activate the environment.

However, note the following constraint: you **cannot** install multiple versions of BMDS Desktop in the same environment. An environment can hold only one version of BMDS Desktop at a time, but you can create multiple environments on your computer.

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

You have many options for running the `pybmds` software. The package is a standard Python package. If you're familiar with Python, you can use any approach for writing code you may already be familiar with.

If you are new to the Python ecosystem, you can use one of the following commonly used development environments - [jupyterlab](https://jupyter.org/)[^disclaimer], Microsoft [Visual Studio Code](https://code.visualstudio.com/)[^disclaimer], or Posit [Positron](https://github.com/posit-dev/positron)[^disclaimer].

If you'd prefer to not write any code, then install `bmds-ui` and start BMDS Desktop.

(using-r)=
### Writing `pybmds` Code in R

You can use use RStudio to execute `pybmds` using the [reticulate](https://rstudio.github.io/reticulate/) package; you'll still need to follow the guide above to install the software on your computer in a Python environment.

See the [using R](recipes/using-r.md) recipe for an example.

[^disclaimer]: Mention of or referral to commercial products or services, and/or links to non-EPA sites does not imply official EPA endorsement of or responsibility for the opinions, ideas, data, or products presented at those locations, or guarantee the validity of the information provided.
