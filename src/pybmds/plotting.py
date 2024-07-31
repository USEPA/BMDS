import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.colors import to_rgba

mpl.use("Agg")  # prevent matplotlib framework issues

__all__ = []

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
FAILURE_MESSAGE_FORMAT = dict(
    style="italic",
    weight="bold",
    bbox={"facecolor": "red", "alpha": 0.35, "pad": 10},
    horizontalalignment="center",
    verticalalignment="center",
)


def create_empty_figure(figsize: tuple[float, float] | None = None):
    plt.style.use("seaborn-v0_8-darkgrid")
    mpl.rcParams.update({"font.size": 10})
    fig, ax = plt.subplots(figsize=figsize or PLOT_FIGSIZE, dpi=DPI)
    return fig


def close_figure(fig):
    plt.close(fig)


def add_bmr_lines(
    ax, bmd: float, bmd_y: float, bmdl: float, bmdu: float, axlines: bool = False, **kw
):
    if bmd <= 0:
        return
    if axlines:
        xmin, xmax = ax.get_xlim()
        xrange = xmax - xmin
        ymin, ymax = ax.get_ylim()
        yrange = ymax - ymin
        ax.axhline(
            y=bmd_y,
            xmin=0,
            xmax=(bmd - xmin) / xrange,
            # label="BMD/BMDL",
            **AXLINE_FORMAT,
        )
        ax.axvline(
            x=bmd,
            ymin=0,
            ymax=(bmd_y - ymin) / yrange,
            # label="BMD/BMDL",
            **AXLINE_FORMAT,
        )
        if bmdl > 0:
            ax.axvline(
                x=bmdl,
                ymin=0,
                ymax=(bmd_y - ymin) / yrange,
                label="BMD/BMDL",
                **AXLINE_FORMAT,
            )
    else:
        styles = {**BMD_LINE_FORMAT, **kw}
        lower = 0 if bmdl < 0 else max(0, bmd - bmdl)
        upper = 0 if bmdu < 0 else max(0, bmdu - bmd)
        ax.errorbar(bmd, bmd_y, label="BMDL-BMD-BMDU", xerr=[[lower], [upper]], **styles)


def improve_bmd_diamond_legend(ax: Axes):
    # manually improve legend icon for BMD diamond
    legend = ax.get_legend()
    rows = legend._legend_handle_box.get_children()[0].get_children()
    for row in rows:
        text = row.get_children()[1].properties()["text"]
        if text == "BMDL-BMD-BMDU":
            area = row.get_children()[0]

            # reduce line height for BMD diamond
            area.get_children()[0].set_linewidth(5)

            # expand line width for BMD diamond
            paths = area.get_children()[0].get_paths()[0]
            paths.vertices[0][0] -= 3
            paths.vertices[1][0] += 3

            # increase diamond size; reduce white edge
            area.get_children()[2].set_markeredgewidth(1)
            area.get_children()[2].set_markersize(8)
