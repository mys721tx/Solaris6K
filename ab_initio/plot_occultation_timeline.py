from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

DEFAULT_INPUT_TSV = Path("moon_occultation_events.tsv")
DEFAULT_OUTPUT_PNG = Path("moon_occultation_timeline.png")


def plot_occultation_timeline(events_df, output_plot_path):
    if events_df.empty:
        return False

    timeline_df = events_df.copy()
    timeline_df["external_ingress_dt"] = pd.to_datetime(
        timeline_df["external_ingress"],
        errors="coerce",
    )
    timeline_df["external_egress_dt"] = pd.to_datetime(
        timeline_df["external_egress"],
        errors="coerce",
    )
    timeline_df["internal_ingress_dt"] = pd.to_datetime(
        timeline_df["internal_ingress"],
        errors="coerce",
    )
    timeline_df["internal_egress_dt"] = pd.to_datetime(
        timeline_df["internal_egress"],
        errors="coerce",
    )

    timeline_df = timeline_df.dropna(
        subset=["external_ingress_dt", "external_egress_dt"]
    )
    if timeline_df.empty:
        return False

    bodies = sorted(timeline_df["body"].unique())
    y_positions = {body: idx for idx, body in enumerate(bodies)}

    fig, ax = plt.subplots(figsize=(16, max(4, len(bodies) * 1.2)))

    for row in timeline_df.itertuples(index=False):
        y = y_positions[row.body]
        external_start = mdates.date2num(row.external_ingress_dt)
        external_end = mdates.date2num(row.external_egress_dt)

        ax.hlines(
            y,
            external_start,
            external_end,
            colors="#4C72B0",
            linewidth=8,
            alpha=0.55,
            zorder=1,
        )
        ax.scatter(
            [external_start, external_end],
            [y, y],
            color="#4C72B0",
            s=8,
            alpha=0.75,
            zorder=3,
        )

        if pd.notna(row.internal_ingress_dt) and pd.notna(row.internal_egress_dt):
            internal_start = mdates.date2num(row.internal_ingress_dt)
            internal_end = mdates.date2num(row.internal_egress_dt)

            ax.hlines(
                y,
                internal_start,
                internal_end,
                colors="#D64F4F",
                linewidth=4,
                alpha=0.95,
                zorder=2,
            )
            ax.scatter(
                [internal_start, internal_end],
                [y, y],
                color="#D64F4F",
                s=8,
                alpha=0.95,
                zorder=4,
            )

    partial_handle = Line2D(
        [0],
        [0],
        color="#4C72B0",
        lw=8,
        alpha=0.35,
        label="External",
    )
    complete_handle = Line2D(
        [0],
        [0],
        color="#D64F4F",
        lw=4,
        alpha=0.55,
        label="Internal",
    )

    ax.legend(
        handles=[partial_handle, complete_handle],
        frameon=True,
    )
    ax.set_title("Moon Occultation Timeline: External and Internal Intervals")
    ax.set_xlabel("Date")
    ax.set_yticks([y_positions[body] for body in bodies])
    ax.set_yticklabels([body.upper() for body in bodies])
    ax.grid(axis="x", linestyle="--", alpha=0.3)
    ax.xaxis.set_major_locator(mdates.AutoDateLocator())
    ax.xaxis.set_major_formatter(
        mdates.ConciseDateFormatter(ax.xaxis.get_major_locator())
    )
    ax.xaxis_date()

    fig.tight_layout()
    fig.savefig(output_plot_path, dpi=300)
    plt.close(fig)
    return True


def main():
    if not DEFAULT_INPUT_TSV.exists():
        raise FileNotFoundError(f"Missing occultation events TSV: {DEFAULT_INPUT_TSV}")

    events_df = pd.read_csv(DEFAULT_INPUT_TSV, sep="\t")
    wrote_plot = plot_occultation_timeline(events_df, DEFAULT_OUTPUT_PNG)
    if wrote_plot:
        print(f"Wrote occultation timeline plot to: {DEFAULT_OUTPUT_PNG}")
    else:
        print("No plottable occultation intervals were found.")


if __name__ == "__main__":
    main()
