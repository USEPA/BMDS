# BMDS Desktop

BMDS Desktop is a graphical user interface to execute dose-response modeling on your computer. It allows you to execute analyses fully offline and store your data in a database file that can be shared with others.

BMDS Desktop is identical to [BMDS Online](https://bmdsonline.epa.gov), with a few additional features:

* Analyses (dose response analyses) and data storage are fully offline
* Database files (projects) are single files containing all analyses
* Within a project, analyses can be labelled and organized

Follow the [installation](installation.md) guide to install the software. Make sure to create the [BMDS Desktop Manager](./installation.md#create-the-bmds-desktop-manager-shortcut) shortcut. To start BMDS Desktop, double-click the shortcut and then enter option `1` to start the application:

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

BMDS Desktop Startup Interface. On startup, you see a list of existing projects that have been run with BMDS Desktop. You can create new projects from this screen as well. Navigate the interface using your keyboard or mouse.
```

Each project in BMDS Desktop contains all of that project's analyses stored in a single file. You can create a single project and store all of your analyses in a single file, or multiple projects - one project per chemical, for example.

### Project Creation and Management

```{figure} _static/img/create-db.jpg
:alt: Screenshot of BMDS Desktop Project Creation

BMDS Desktop Project Creation. Create a new project by specifying a path and a database filename. The path must already exist on your computer; you can copy and paste the path from a file manager such as Windows Explorer.
```

:::{important}
Creating a new project creates an accompanying database file at the path and filename specified if it does not already exist. Database files should have the `.db` extension. BMDS Desktop also creates other files with different extensions in that directory, such as `.db-shm` and `.db-wal`. **Do not delete those files -** they allow multiple users to work with the same project concurrently.
:::

You can also update a database's location, in the event that you moved the database to a new directory. Updating the database refreshes the path/location so BMDS Desktop can find it again; the updating process does not change the file's contents.

Deleting a project from the BMDS Desktop Startup Interface deletes its entry in the list of recent databases, but does not delete the database itself. To fully delete the database files, navigate to that project in your system's file manager and manually delete.

## BMDS Desktop Application

After at least one project has been created, select the "Start" button to run the project. Starting a project will open a new BMDS Desktop tab in your default browser. For a new project, initial startup may take up to a minute before the browser tab appears.

You can run only one BMDS Desktop project at a time.

```{figure} _static/img/bmds-desktop.jpg
:alt: Screenshot of BMDS Desktop Application

BMDS Desktop Application. The main home page of the BMDS Desktop application. From BMDS Desktop you can create and execute dose-response analyses, and generate summary reports. 
```
BMDS Desktop allows users to organize all the analyses within a single project by applying "Collection" tags to the analyses and starring analyses to indicate which ones are most important.  Both starring and organizing into Collections allows users to filter their analyses as they see fit.  

For example, using the Collection dropdown menu on the main page, users can create Collections that allow them to organize analyses by organ system, study name, sex, strain, or any descriptive tag they require.

```{figure} _static/img/bmds-desktop-collection-dropdown.jpg
:alt: Screenshot of BMDS Desktop Main Page with Collections drop down menu displayed with different tag options

The Collections drop down menu on the BMDS Desktop Main Page.
```

Once a tag is created within the Collections drop down menu on the Main Page, it can be applied to an analysis by clicking on an analysis and then applying the tag using the Actions button.


```{figure} _static/img/bmds-desktop-collections-actions-button.jpg
:alt: Screenshot of BMDS Desktop Settings tab for an analysis with the Collection field within the Actions drop down menu displayed

Collection tags can be applied to individual analyses from the Actions drop down menu.
```
Once a Collections tag has been applied to an analysis, that tag will display next to the analysis name on the Main Page.

```{figure} _static/img/bmds-desktop-collection-tag-applied.jpg
:alt: Screenshot of BMDS Desktop Main Page with Collections tags applied to all analyses

Collection tags display next to analyses once applied.
```

Analyses can also be starred using the Starred toggle button.  Saved analyses can be searched by Analysis name or sorted by Collections tags or whether they are starred.


### Modeling Dose Response Data

Please see the BMDS User Guide for more information on dose-response modeling and execution at [https://www.epa.gov/bmds](https://www.epa.gov/bmds).
