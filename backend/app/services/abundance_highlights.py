"""Selection and annotations for rank-abundance plots, including portable exports."""
import re

import numpy as np


def abundance_selection(frame, search=None, top_high=0, top_low=0):
    """Select plotted rows, independently within each condition; never change ranks."""
    if (search is None or len(search) == 0) and top_high <= 0 and top_low <= 0:
        return frame["Rank"].notna() & False

    def normalize(value):
        return re.sub(r"_(human|mouse)$", "", str(value).strip().lower())

    entries = [search] if isinstance(search, str) else ([] if search is None else search)
    terms = {normalize(term) for entry in entries for term in str(entry).split()}
    selected = frame["Feature"].map(lambda value: normalize(value) in terms or any(
        normalize(token) in terms for token in re.split(r"[;,]", str(value))))
    finite = np.isfinite(frame["log10Intensity"])
    for _, group in frame.loc[finite].groupby("Condition", observed=True, sort=False):
        ordered = group.sort_values("Rank", kind="stable")
        selected.loc[ordered.head(max(0, int(top_high))).index] = True
        selected.loc[ordered.tail(max(0, int(top_low))).index] = True
    return selected & finite


def abundance_highlight_style(advanced):
    settings = advanced or {}
    return {
        "width": max(.1, float(settings.get("abundanceHighlightWidth") or 1)),
        "color": settings.get("abundanceHighlightColor") or "#333333",
        "legend": settings.get("abundanceHighlightLegend") is not False,
        "name": settings.get("abundanceHighlightLegendLabel") or "Highlighted",
        "label_size": settings.get("abundanceLabelSize") or 7,
    }


def abundance_static_highlights(fig, ax, frame, colors, names, search, top_high, top_low,
                                outline, labels, advanced=None, legend=True):
    selected = frame.loc[abundance_selection(frame, search, top_high, top_low)]
    style = abundance_highlight_style(advanced)
    if outline and not selected.empty:
        sizes = {}
        for condition in selected["Condition"].unique():
            base = next((collection for collection in ax.collections if collection.get_label() == names[condition]), None)
            diameter = float(np.sqrt(base.get_sizes()[0])) if base is not None and len(base.get_sizes()) else np.sqrt(5)
            sizes[condition] = diameter ** 2
        ax.scatter(selected["Rank"], selected["log10Intensity"], facecolors="none",
                   edgecolors=[style["color"] or colors[cond] for cond in selected["Condition"]],
                   s=[sizes[cond] for cond in selected["Condition"]], linewidths=style["width"],
                   zorder=4, label=style["name"] if style["legend"] and legend else "_nolegend_")
    if not labels or selected.empty:
        return
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    legend_artist = ax.get_legend()
    placed = [legend_artist.get_window_extent(renderer)] if legend_artist is not None and legend_artist.get_visible() else []
    multiple = frame["Condition"].nunique() > 1
    for row in selected.itertuples():
        label = str(row.Feature) + (f" ({names[row.Condition]})" if multiple else "")
        color = colors[row.Condition]
        for dx, dy in ((6, 6), (6, -14), (-6, 6), (-6, -14), (6, 26), (6, -34), (-6, 26), (-6, -34), (6, 46), (-6, 46), (6, -54), (-6, -54), (-6, 66)):
            text = ax.annotate(label, (row.Rank, row.log10Intensity), xytext=(dx, dy),
                               textcoords="offset points", color=color, fontsize=style["label_size"],
                               ha="left" if dx > 0 else "right", zorder=5,
                               bbox={"facecolor": ax.get_facecolor(), "edgecolor": "none", "alpha": .85, "pad": 0},
                               arrowprops={"arrowstyle": "-", "color": color, "linewidth": .5})
            text.update_positions(renderer)
            text.update_bbox_position_size(renderer)
            box = text.get_bbox_patch().get_window_extent(renderer).expanded(1.04, 1.35)
            if (ax.bbox.contains(box.x0, box.y0) and ax.bbox.contains(box.x1, box.y1)
                    and not any(box.overlaps(other) for other in placed)):
                placed.append(box)
                break
            text.remove()


def abundance_interactive_highlights(fig, frame, colors, names, search, top_high, top_low,
                                     outline, labels, advanced=None):
    """Keep rings in their condition legend group; interactive labels use hover."""
    labels = False
    selected = frame.loc[abundance_selection(frame, search, top_high, top_low)]
    style = abundance_highlight_style(advanced)
    multiple = frame["Condition"].nunique() > 1
    for condition, group in selected.groupby("Condition", observed=True, sort=False):
        if not outline and not labels:
            continue
        name = names[condition]
        base = next((trace for trace in fig.data if trace.legendgroup == name), None)
        diameter = base.marker.size if base is not None and base.marker.size is not None else 6
        size = float(np.max(diameter))
        fig.add_scatter(
            x=group["Rank"], y=group["log10Intensity"],
            mode="markers+text" if outline and labels else "markers" if outline else "text",
            marker={"size": size, "color": style["color"] or colors[condition], "symbol": "circle-open", "line": {"width": style["width"]}},
            text=[str(feature) + (f" ({name})" if multiple else "") for feature in group["Feature"]] if labels else None,
            textposition="top center", textfont={"color": colors[condition], "size": style["label_size"] * 96 / 72},
            customdata=[[str(feature), str(name)] for feature in group["Feature"]],
            hovertemplate="%{customdata[0]}<br>Condition: %{customdata[1]}<br>Rank: %{x}<br>log10 intensity: %{y:.4f}<extra></extra>",
            name=f"{name} ({style['name']})", legendgroup=name,
            showlegend=False,
        )

    if outline and style["legend"] and not selected.empty:
        fig.add_scatter(x=[None], y=[None], mode="markers", hoverinfo="skip",
            marker={"size": 6, "symbol": "circle-open", "color": style["color"], "line": {"width": style["width"]}},
            name=style["name"], legendgroup="_copylot_highlight_legend", showlegend=True)
