import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from .constants import AXLINE_FORMAT, BMD_LINE_FORMAT, DPI, PLOT_FIGSIZE


def create_empty_figure(
    rows: int = 1, cols: int = 1, figsize: tuple[float, float] | None = None
) -> Figure:
    plt.style.use("seaborn-v0_8-darkgrid")
    mpl.rcParams.update({"font.size": 10})
    fig, _ = plt.subplots(rows, cols, figsize=figsize or PLOT_FIGSIZE, dpi=DPI)
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
            label="BMD/BMDL",
            **AXLINE_FORMAT,
        )
        ax.axvline(
            x=bmd,
            ymin=0,
            ymax=(bmd_y - ymin) / yrange,
            **AXLINE_FORMAT,
        )
        if bmdl > 0:
            ax.axvline(
                x=bmdl,
                ymin=0,
                ymax=(bmd_y - ymin) / yrange,
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
