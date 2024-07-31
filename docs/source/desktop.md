# BMDS Desktop

BMDS Desktop is a graphical user interface to execute dose-response modeling on your computer. It allows you to execute analyses fully offline and store your data in a database file that can be shared with others.

BMDS Desktop is identical to [BMDS Online](https://bmdsonline.epa.gov), with a few additional features:

* Analyses and data storage is fully offline
* Database files are single files containing all analyses
* Within a project, analyses can be labelled and organized

After following the [installation](installation.md) guide, start from your terminal:

```bash
bmds-desktop
```

:::{tip}
Make sure to [activate](./installation.md/#activating-an-environment) the environment prior to running this command, otherwise it won't be able to find the command.
:::

## BMDS Desktop Startup Interface

The BMDS Desktop Startup Interface is the gateway to create a BMDS project and start the application.

```{figure} _static/img/desktop-startup.jpg
:alt: Screenshot of BMDS Desktop Startup

BMDS Desktop Startup Interface. A list of projects show existing projects which have been run with this version of BMDS Desktop; you can create new projects as well. You can use your keyboard or mouse to navigate the interface.
```

Each project in BMDS Desktop contains all analyses in that project stored in a single database file. You can create a single project and store all of your analyses in a single file, or multiple projects, one per chemical analysis for example.

### Project Creation and Management

```{figure} _static/img/create-db.jpg
:alt: Screenshot of BMDS Desktop Project Creation

BMDS Desktop Project Creation. Create a new project by specifying a database filename. It must be in a path that already exists on your folder; you can copy-paste the path from a file explorer as well.
```

:::{important}
When you create a new project, it creates a database file, for example, `bmds-database.db`, in the location specified. It also makes a few other files with different extension in that directory, `.db-shm` and `.db-wal`, do not delete those files. They allow multiple people to work with the same project concurrently.
:::

You can also update a database as well. Note that updating the database filename or path doesn't actually change the file; it just allows you select the same file in case it moved.

If you delete a project, it will delete the project from the startup interface, but it will not delete the database project itself. You can delete that file manually in a file explorer, should you wish.

## BMDS Desktop Application

After at least one project has been created, you're can run the Desktop Application creating a project, you can "Start" the project, which will create a new tab in your browser window and open BMDS Desktop. It may take a a few seconds to startup the first time being used. You can only run one BMDS Desktop project at a time.

```{figure} _static/img/bmds-desktop.jpg
:alt: Screenshot of BMDS Desktop Application

BMDS Desktop Application. The main home page of the BMDS Desktop application. From BMDS Desktop you can create and execute dose response analyses, as well as generate summary reports. You can also star and label analyses and search and filter by the analysis name as well as other information.
```

### Modeling Dose Response Data

Please see the BMDS User Guide for more information on dose-response modeling and execution at https://www.epa.gov/bmds.
