from matplotlib.colors import to_rgba

PLOT_FIGSIZE = (8, 5)
DPI = 100
PLOT_MARGINS = 0.05
DATASET_POINT_FORMAT = dict(ms=7, fmt="o", c="k", capsize=3, lw=1, zorder=50)
DATASET_INDIVIDUAL_FORMAT = dict(s=35, color=to_rgba("#ffffff", 0.5), edgecolors="black")
LEGEND_OPTS = dict(loc="best", fontsize=9, frameon=True, facecolor="white", markerscale=0.5)
LINE_FORMAT = dict(c="#6470C0", lw=3, zorder=100)
AXLINE_FORMAT = dict(c="#b8a800", lw=3, zorder=100)
INDIVIDUAL_MODEL_COLORS = ["#6e40aa", "#e7298a", "#1b9e77", "#b8a800", "#666666"]
INDIVIDUAL_LINE_STYLES = ["solid", "dotted", "dashed", "dashdot"]
CDF_SOLID = dict(c="#b8a800", lw=3, zorder=50)
CDF_DASHED = dict(c="#b8a800", lw=3, ls=(1, (2, 1)), zorder=50)
BMD_LABEL_FORMAT = dict(size=9)
BMD_LINE_FORMAT = dict(
    c="#6470C0",
    markeredgecolor="white",
    markeredgewidth=2,
    fmt="d",
    ecolor=to_rgba("#6470C0", 0.7),
    ms=12,
    elinewidth=7,
    zorder=150,
)
