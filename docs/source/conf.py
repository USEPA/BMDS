project = "pybmds"
copyright = "MIT License"

extensions = ["myst_nb", "sphinx_design"]
souce_suffix = {
    ".md": "markdown",
    ".rst": "restructuredtext",
}
templates_path = ["_templates"]

exclude_patterns = []

myst_links_external_new_tab = True
myst_enable_extensions = [
    "attrs_inline",
    "colon_fence",
    "dollarmath",
    "fieldlist",
    "html_image",
    "smartquotes",
]

# HTML settings
html_theme = "furo"
html_static_path = ["_static"]
html_css_files = ["css/style.css"]
html_theme_options = {
    "sidebar_hide_name": True,
    "light_logo": "img/pybmds.png",
    "dark_logo": "img/pybmds.png",
    "source_repository": "https://github.com/USEPA/bmds",
    "source_branch": "main",
    "source_directory": "docs/",
}

# Latex / PDF settings
latex_elements = {
    "printindex": "",
    "sphinxsetup": "hmargin={0.9in,0.9in}, vmargin={0.9in,0.9in}, marginpar=1.0in",
    "papersize": "letterpaper",
    "pointsize": "10pt",
    "figure_align": "htbp",
}
