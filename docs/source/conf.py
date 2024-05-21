project = "pybmds"
copyright = "Public Domain"
author = "U.S. EPA"

extensions = ["myst_nb"]
souce_suffix = {
    ".md": "markdown",
    ".rst": "restructuredtext",
}
templates_path = ["_templates"]

exclude_patterns = []

# HTML settings
html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]
html_css_files = ["css/style.css"]
html_sidebars = {"**": ["globaltoc.html"]}
html_theme_options = {
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/USEPA/bmds-private",
            "icon": "fa-brands fa-github",
            "type": "fontawesome",
        },
        {
            "name": "U.S. Environmental Protection Agency",
            "url": "https://epa.gov/bmds",
            "icon": "_static/img/epa_logo.png",
            "type": "local",
        },
    ],
    "use_edit_page_button": False,
    "show_toc_level": 2,
    "navbar_end": [
        "theme-switcher.html",
        "navbar-icon-links.html",
    ],
}

# Latex / PDF settings
latex_elements = {
    "printindex": "",
    "sphinxsetup": "hmargin={0.9in,0.9in}, vmargin={0.9in,0.9in}, marginpar=1.0in",
    "papersize": "letterpaper",
    "pointsize": "10pt",
    "figure_align": "htbp",
}
