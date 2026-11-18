from __future__ import annotations

import io
import math
from contextvars import ContextVar
from threading import RLock

plot_output_format: ContextVar[str] = ContextVar("plot_output_format", default="png")
# False while a request asks for a download in another format: that render is a
# file to save, not the figure on screen, so it gets no notebook export of its own.
notebook_capture_allowed: ContextVar[bool] = ContextVar("notebook_capture_allowed", default=True)
# Why this render offers no notebook export. The request middleware puts a
# fresh dict here; the route runs in a worker thread with a copy of the context,
# so the reason is written into the shared dict rather than set on the variable.
notebook_unavailable: ContextVar[dict | None] = ContextVar("notebook_unavailable", default=None)


def note_notebook_unavailable(reason) -> None:
    """Say, under the plot, why it has no notebook export (GUI only)."""
    holder = notebook_unavailable.get()
    text = str(reason).strip()
    if holder is not None and text and "reason" not in holder:
        holder["reason"] = text


plot_notices: ContextVar[list | None] = ContextVar("plot_notices", default=None)


def note_plot_notice(text: str) -> None:
    """Say something about this plot under it (GUI only, never in the figure or report)."""
    holder = plot_notices.get()
    text = str(text).strip()
    if holder is not None and text and text not in holder:
        holder.append(text)

# Formats matplotlib writes as vectors rather than pixels.
_VECTOR_FORMATS = frozenset({"pdf", "svg"})
_SAVE_LOCK = RLock()
# One figure is built at a time. The server answers requests on several
# threads, and matplotlib is not thread-safe: its math-text parser (used for
# every tick label) keeps parser state, so two figures laid out at once could
# fail with an internal error ("pop from empty list") or corrupt each other.
# Re-entrant, and the same lock as the save, so a render that saves holds it
# once.
RENDER_LOCK = _SAVE_LOCK


class PlotBytes(bytes):
    """Image payload with a notice kept outside the exported artwork."""

    def __new__(cls, content: bytes, notice: str = ""):
        result = super().__new__(cls, content)
        result.notice = notice
        return result

# Upper bounds for a figure that is allowed to grow past the requested size.
# They mirror MAX_PLOT_WIDTH_CM / MAX_PLOT_HEIGHT_CM / MAX_PLOT_PIXELS in
# plot_limits, so an auto-grown canvas can never end up larger than a canvas the
# user is allowed to ask for. Duplicated rather than imported because
# plot_limits pulls in Starlette and this module is used off the HTTP path too.
MAX_AUTO_FIT_CM = 60.0
MAX_AUTO_FIT_PIXELS = 24_000_000
# Breathing room added around the measured label extent, in inches.
AUTO_FIT_MARGIN_INCHES = 0.12


def _tick_label_space(fig, renderer) -> tuple[float, float]:
    """Extra inches needed to separate adjacent tick labels, without hiding any."""
    extra = [0.0, 0.0]
    for ax in getattr(fig, "axes", []):
        for dimension, labels in enumerate((ax.get_xticklabels(), ax.get_yticklabels())):
            boxes = [label.get_window_extent(renderer) for label in labels
                     if label.get_visible() and label.get_text()]
            boxes.sort(key=lambda box: box.x0 if dimension == 0 else box.y0)
            ratio = 1.0
            for first, second in zip(boxes, boxes[1:]):
                first_center = (first.x0 + first.x1) / 2 if dimension == 0 else (first.y0 + first.y1) / 2
                second_center = (second.x0 + second.x1) / 2 if dimension == 0 else (second.y0 + second.y1) / 2
                distance = second_center - first_center
                if distance <= 0:
                    continue
                extent = (first.width + second.width) / 2 if dimension == 0 else (first.height + second.height) / 2
                ratio = max(ratio, (extent + fig.dpi * 0.02) / distance)
            span = ax.bbox.width if dimension == 0 else ax.bbox.height
            extra[dimension] = max(extra[dimension], span * (ratio - 1) / fig.dpi)
    return tuple(extra)


def _fit_crowded_legend(fig, dpi: int) -> None:
    """Reflow a single crowded axes legend without shrinking its text or data."""
    legends = [(ax, ax.get_legend()) for ax in getattr(fig, "axes", [])
               if ax.get_legend() is not None]
    if len(legends) != 1 or getattr(fig, "_copylot_legend_bottom", None) is not None:
        return
    ax, old = legends[0]
    if getattr(old, "_copylot_manual_position", False):
        return
    renderer = fig.canvas.get_renderer()
    width, height = fig.get_size_inches()
    bounds = old.get_window_extent(renderer)
    if len(old.get_texts()) < 8 or bounds.height / fig.dpi < height * 0.55:
        return
    texts = list(old.get_texts())
    handles = old.legend_handles
    # Measure candidate columns, including handles and inter-column padding.
    # Rebuilding through the public API preserves handler-generated artists.
    best = 1
    for columns in range(1, len(texts) + 1):
        candidate = ax.legend(handles, [t.get_text() for t in texts],
                              ncols=columns, prop=texts[0].get_fontproperties(),
                              title=old.get_title().get_text(),
                              title_fontproperties=old.get_title().get_fontproperties(),
                              loc="lower center", bbox_to_anchor=(0.5, 0.02),
                              bbox_transform=fig.transFigure)
        if candidate.get_window_extent(renderer).width / fig.dpi > width - 0.3:
            break
        best = columns
    legend = ax.legend(handles, [t.get_text() for t in texts], ncols=best,
                       prop=texts[0].get_fontproperties(),
                       title=old.get_title().get_text(),
                       title_fontproperties=old.get_title().get_fontproperties(),
                       loc="lower center", bbox_to_anchor=(0.5, 0.02),
                       bbox_transform=fig.transFigure, frameon=old.get_frame_on(),
                       markerscale=1, handlelength=old.handlelength,
                       handletextpad=old.handletextpad, columnspacing=old.columnspacing,
                       labelspacing=old.labelspacing, borderpad=old.borderpad)
    for original, replacement in zip(texts, legend.get_texts()):
        replacement.set_fontproperties(original.get_fontproperties())
        replacement.set_color(original.get_color())
    legend.get_title().set_color(old.get_title().get_color())
    legend.get_frame().update_from(old.get_frame())
    legend_height = legend.get_window_extent(renderer).height / fig.dpi + 0.35
    # Keep at least the usual 65% of the requested height for the axes, with
    # another inch for title/ticks. Grow only when the legend needs that room.
    target_height = max(height, height * 0.65 + 1 + legend_height)
    target_height = min(target_height, MAX_AUTO_FIT_CM / 2.54,
                        MAX_AUTO_FIT_PIXELS / max(72, dpi)**2 / width)
    fig.set_size_inches(width, max(height, target_height))
    fig._copylot_legend_bottom = legend_height
    # The legend is in figure coordinates. Reserve its measured space rather
    # than letting tight_layout treat it as an axes decoration.
    legend.set_in_layout(False)


def _tight_plot_layout(fig) -> None:
    bottom = getattr(fig, "_copylot_legend_bottom", None)
    if bottom is None:
        fig.tight_layout(pad=0.8)
    else:
        fig.tight_layout(pad=0.8, rect=(0, bottom / fig.get_size_inches()[1], 1, 1))


def fit_figure_to_labels(
    fig,
    plt,
    *,
    dpi: int = 150,
    max_width_cm: float = MAX_AUTO_FIT_CM,
    max_height_cm: float = MAX_AUTO_FIT_CM,
    margin_inches: float = AUTO_FIT_MARGIN_INCHES,
    max_passes: int = 3,
    max_pixels: int = MAX_AUTO_FIT_PIXELS,
) -> tuple[float, float]:
    """Grow a crowded figure within the physical/pixel limits, preserving DPI.

    Requested dimensions are a minimum. Return the final size in inches.
    """
    resolved_dpi = max(72, int(dpi))
    tolerance = 1.0 / resolved_dpi
    max_width_in = cm_to_inch(max_width_cm)
    max_height_in = cm_to_inch(max_height_cm)
    pixel_budget = max(1, int(max_pixels))
    width, height = (float(value) for value in fig.get_size_inches())

    with _SAVE_LOCK, plt.rc_context({"savefig.bbox": None}):
        for _ in range(max(1, int(max_passes))):
            _tight_plot_layout(fig)
            fig.canvas.draw()
            bounds = fig.get_tightbbox(fig.canvas.get_renderer())
            if bounds is None:
                break
            width, height = (float(value) for value in fig.get_size_inches())
            overflow_x = max(0.0, -bounds.x0) + max(0.0, bounds.x1 - width)
            overflow_y = max(0.0, -bounds.y0) + max(0.0, bounds.y1 - height)
            tick_width, tick_height = _tick_label_space(fig, fig.canvas.get_renderer())
            overflow_x += tick_width
            overflow_y += tick_height
            if overflow_x <= tolerance and overflow_y <= tolerance:
                break
            target_width = width + (overflow_x + margin_inches if overflow_x > tolerance else 0.0)
            target_height = height + (overflow_y + margin_inches if overflow_y > tolerance else 0.0)
            target_width = min(target_width, max_width_in)
            target_height = min(target_height, max_height_in)
            # A tall canvas at 600 dpi can blow past the pixel budget long
            # before it reaches the physical ceiling; scale the growth back
            # rather than handing the browser a figure it cannot decode.
            pixels = (target_width * resolved_dpi) * (target_height * resolved_dpi)
            if pixels > pixel_budget:
                shrink = math.sqrt(pixel_budget / pixels)
                target_width = max(width, target_width * shrink)
                target_height = max(height, target_height * shrink)
                # Preserving an existing dimension can undo the proportional
                # reduction above. Cap the other dimension to the remaining budget.
                target_width = min(target_width, pixel_budget / resolved_dpi**2 / height)
                target_height = min(target_height, pixel_budget / resolved_dpi**2 / target_width)
            if target_width <= width + tolerance and target_height <= height + tolerance:
                break
            fig.set_size_inches(max(width, target_width), max(height, target_height))
            width, height = (float(value) for value in fig.get_size_inches())

        # Growing the canvas does not always close the gap: tight_layout gives
        # the extra room straight back to the axes, so a long y label can still
        # sit a hundredth of an inch outside after several passes, and the guard
        # in `to_png_bytes` refuses it. What is left is absorbed by pulling the
        # axes in, which removes the overflow by construction rather than by
        # hoping the next pass converges.
        for _ in range(2):
            fig.canvas.draw()
            bounds = fig.get_tightbbox(fig.canvas.get_renderer())
            if bounds is None:
                break
            width, height = (float(value) for value in fig.get_size_inches())
            left = max(0.0, -bounds.x0)
            right = max(0.0, bounds.x1 - width)
            bottom = max(0.0, -bounds.y0)
            top = max(0.0, bounds.y1 - height)
            if max(left, right, bottom, top) <= tolerance:
                break
            pad = margin_inches * 0.5
            box = fig.subplotpars
            new_left = box.left + (left + pad) / width
            new_right = box.right - (right + pad) / width
            new_bottom = box.bottom + (bottom + pad) / height
            new_top = box.top - (top + pad) / height
            # An impossible layout -- labels wider than any canvas -- would ask
            # for margins that leave no axes at all. Leave it alone rather than
            # raising here: the guard in `to_png_bytes` is the place that says
            # a figure does not fit, and it says so with a message the user can
            # act on.
            if new_right - new_left < 0.2 or new_top - new_bottom < 0.2:
                break
            fig.subplots_adjust(
                left=new_left, right=new_right, bottom=new_bottom, top=new_top
            )
        # The layout measured here is the one that must be rendered, so the
        # serializer must not run tight_layout over it again.
        fig._copylot_layout_final = True

    return width, height


def to_png_bytes(fig, plt, dpi: int = 150, tight: bool = True) -> bytes:
    """Serialize without clipping labels; grow crowded plots within safe limits."""
    buf = io.BytesIO()
    fmt = plot_output_format.get()
    save_kwargs = {"format": fmt, "dpi": max(72, int(dpi))}
    if getattr(fig, "_copylot_transparent_background", False):
        save_kwargs["transparent"] = True
    try:
        if fmt in _VECTOR_FORMATS:
            # A rasterized artist would drop a bitmap into a vector file, which
            # defeats the reason for asking for one.
            for artist in fig.findobj():
                if artist.get_rasterized():
                    artist.set_rasterized(False)
        with _SAVE_LOCK, plt.rc_context({"pdf.fonttype": 42, "savefig.bbox": None}):
            fig.canvas.draw()
            _fit_crowded_legend(fig, save_kwargs["dpi"])
            if getattr(fig, "_copylot_legend_bottom", None) is not None or (
                tight and not getattr(fig, "_copylot_layout_final", False)
            ):
                _tight_plot_layout(fig)
            fig.canvas.draw()
            bounds = fig.get_tightbbox(fig.canvas.get_renderer())
            width, height = fig.get_size_inches()
            tolerance = 1 / save_kwargs["dpi"]
            crowded_ticks = max(_tick_label_space(fig, fig.canvas.get_renderer())) > tolerance
            if crowded_ticks or (bounds is not None and (
                bounds.x0 < -tolerance or bounds.y0 < -tolerance
                or bounds.x1 > width + tolerance or bounds.y1 > height + tolerance
            )):
                fit_figure_to_labels(fig, plt, dpi=save_kwargs["dpi"])
                fig.canvas.draw()
                bounds = fig.get_tightbbox(fig.canvas.get_renderer())
                width, height = fig.get_size_inches()
                if bounds is not None and (
                    bounds.x0 < -tolerance or bounds.y0 < -tolerance
                    or bounds.x1 > width + tolerance or bounds.y1 > height + tolerance
                ):
                    raise ValueError(
                        "Plot labels do not fit even after automatic layout adjustment. "
                        "Increase width/height or reduce text size."
                    )
            if getattr(fig, "_copylot_legend_bottom", None) is not None:
                renderer = fig.canvas.get_renderer()
                for ax in fig.axes:
                    legend = ax.get_legend()
                    if legend is None:
                        continue
                    extent = legend.get_window_extent(renderer).transformed(fig.dpi_scale_trans.inverted())
                    if (extent.x0 < -tolerance or extent.y0 < -tolerance
                            or extent.x1 > width + tolerance or extent.y1 > height + tolerance):
                        raise ValueError("Plot legend does not fit. Increase width/height or reduce text size.")
            fig.savefig(buf, **save_kwargs)
    finally:
        plt.close(fig)
    notice = getattr(fig, "_copylot_notice", "")
    return PlotBytes(buf.getvalue(), notice) if notice else buf.getvalue()


def cm_to_inch(value_cm: float) -> float:
    value = float(value_cm)
    if not math.isfinite(value) or value <= 0:
        raise ValueError("Figure dimensions must be finite and greater than zero.")
    return value / 2.54


def cm_to_px(value_cm: float) -> int:
    return max(320, int(round((cm_to_inch(value_cm) * 96.0))))
