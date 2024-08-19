# BMDS Desktop

BMDS Desktop is a graphical user interface to execute dose-response modeling on your computer. It allows you to execute analyses fully offline and store your data in a database file that can be shared with others.

BMDS Desktop is identical to [BMDS Online](https://bmdsonline.epa.gov), with a few additional features:

* Analyses (dose response analyses) and data storage is fully offline
* Database files (projects) are single files containing all analyses
* Within a project, analyses can be labelled and organized

Follow the [installation](installation.md) guide for to install the software. Make sure to create the [BMDS Desktop Manager](./installation.md#the-bmds-desktop-manager) shortcut. To start BMDS Desktop, double-click the shortcut and then enter option `1` to start the application:

```{figure} _static/img/bmds-desktop-manager.jpg
:alt: BMDS Desktop Manager

Start page for the BMDS Desktop Manager. Type a number to execute the command, or type q to exit the manager.  From the manager page you can start the BMDS Desktop Interface, or update BMDS Desktop, or get diagnostic version information for troubleshooting.
```

If you are a more experienced developer or prefer to start from your terminal directly, [activate](./installation.md/#activating-an-environment) your environment and then run the command:

```bash
bmds-desktop
```

## BMDS Desktop Startup Interface

The BMDS Desktop Startup Interface is the gateway to create a BMDS project and start the application.

```{figure} _static/img/desktop-startup.jpg
:alt: Screenshot of BMDS Desktop Startup

BMDS Desktop Startup Interface. On startup, you see a list existing projects that have been run with BMDS Desktop; you can create new projects from this screen as well. Navigate the interface using your keyboard or mouse.
```

Each project in BMDS Desktop contains all of that project's analyses stored in a single file. You can create a single project and store all of your analyses in a single file, or multiple projects - one project per chemical, for example.

### Project Creation and Management

```{figure} _static/img/create-db.jpg
:alt: Screenshot of BMDS Desktop Project Creation

BMDS Desktop Project Creation. Create a new project by specifying a path and a database filename. The path must already exist on your computer; you can copy and paste the path from a file manager such as Windows Explorer.
```

:::{important}
Creating a new project creates an accompanying database file at the with the path and filename specified if it does not already exist. Database files should have the `.db` extension. BMDS Desktop also creates other files with different extensions in that directory- `.db-shm` and `.db-wal`. **Do not delete those files.** They allow multiple users to work with the same project concurrently.
:::

You can also update a database's location, in th e event that you moved the database to a new directory. Updating the database refreshes the path/location so BMDS Desktop can find it again; the updating process does not change the file's contents.

Deleting a project from the BMDS Desktop Startup Interface deletes its entry in the list of recent databases, but does not delete the database itself. To fully delete the database files, navigate to that project in your system's file manager and manually delete.

## BMDS Desktop Application

After at least one project has been created, press the "Start"
 button to run the project. Starting a project will open a new BMDS Desktop tab in your default browser. For a new project, initial startup may take up to a minute for the browser tab to appear.

You can run only one BMDS Desktop project at a time.

```{figure} _static/img/bmds-desktop.jpg
:alt: Screenshot of BMDS Desktop Application

BMDS Desktop Application. The main home page of the BMDS Desktop application. From BMDS Desktop you can create and execute dose-response analyses, and generate summary reports. You can also star and label analyses, and search and filter by the analysis name or other information.
```

### Modeling Dose Response Data

Please see the BMDS User Guide for more information on dose-response modeling and execution at [https://www.epa.gov/bmds](https://www.epa.gov/bmds).
