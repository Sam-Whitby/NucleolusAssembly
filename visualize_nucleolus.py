#!/usr/bin/env python3
"""
visualize_nucleolus.py

Reads a nucleolus trajectory file produced by run_nucleolus and generates a
3-panel animated figure:
  Left:         particle positions coloured by polymer type with gradient overlay
  Top-right:    energy vs simulation step
  Bottom-right: cumulative exited complexes vs simulation step

Usage
-----
    # Watch a running simulation in real time (no file saved):
    python3 visualize_nucleolus.py <traj_file> --live

    # Display completed trajectory interactively:
    python3 visualize_nucleolus.py <traj_file>

    # Save to MP4 (recommended — much smaller/faster than GIF):
    python3 visualize_nucleolus.py <traj_file> --output run.mp4 --fps 10

Options
    --live               Real-time mode: poll file for new frames as simulation runs.
    --output FILE        Save animation to FILE (.mp4 strongly recommended over .gif).
    --fps N              Frames per second (default 5).
    --gradient-length L  Draw gradient overlay from x=0 to x=L (default: auto).
    --width W            Column width in y (default: auto).
    --skip N             Render only every N-th frame (default 1, batch mode only).
    --title TEXT         Figure title.

Trajectory format (written by run_nucleolus)
--------------------------------------------
    <N_particles>
    step=S energy=E exited=X L=LL W=WW nCopies=C
    <id> <poly_type> <x> <y> <copy>
    ...   (repeated N_particles times per frame)
"""

import argparse
import re
import sys

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.animation as animation
import numpy as np


_POLY_COLORS   = ["tab:blue", "tab:orange", "tab:green", "tab:red"]
_BACKBONE_COLOR = "#333333"
_BOND_LW        = 1.5


# --------------------------------------------------------------------------- #
# Parsing
# --------------------------------------------------------------------------- #

def _kv(header, key, default=None):
    m = re.search(rf'{key}=([^\s]+)', header)
    return m.group(1) if m else default


def parse_traj(path):
    """
    Return a list of complete frame dicts, each with:
        step, energy, exited, L, W, nCopies,
        particles: list of (id, poly_type, x, y, copy)
    Partial frames at the end of a file (e.g. while simulation is still running)
    are silently discarded so live mode never displays incomplete snapshots.
    """
    frames = []
    try:
        with open(path) as fh:
            lines = fh.readlines()
    except OSError:
        return frames

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        try:
            n = int(line)
        except ValueError:
            i += 1
            continue

        i += 1
        if i >= len(lines):
            break
        hdr = lines[i].strip()
        i += 1

        step    = int(_kv(hdr, 'step',    0))
        energy  = float(_kv(hdr, 'energy', 0.0))
        exited  = int(_kv(hdr, 'exited',  0))
        L       = float(_kv(hdr, 'L',     60))
        W       = float(_kv(hdr, 'W',     10))
        nCopies = int(_kv(hdr, 'nCopies', 4))

        particles = []
        for _ in range(n):
            if i >= len(lines):
                break
            parts = lines[i].strip().split()
            i += 1
            if len(parts) >= 5:
                particles.append((
                    int(parts[0]), int(parts[1]),
                    float(parts[2]), float(parts[3]),
                    int(parts[4]),
                ))

        if len(particles) == n:  # discard partial frames
            frames.append(dict(
                step=step, energy=energy, exited=exited,
                L=L, W=W, nCopies=nCopies, particles=particles,
            ))

    return frames


# --------------------------------------------------------------------------- #
# Figure construction (shared between batch and live modes)
# --------------------------------------------------------------------------- #

def _build_figure(grad_L, col_W, title, x_max=None):
    """
    Build the 3-panel figure layout.  Returns (fig, ax_anim, ax_ener, ax_exits).
    Static decorations (gradient background, vline at L, legend) are drawn here.
    Dynamic artists (scatter, bonds, ts lines) are created separately by
    _init_artists().
    """
    if x_max is None:
        x_max = grad_L + 5

    x_min = -1.0
    y_min = -0.5
    y_max = col_W + 0.5

    SCALE   = 0.38
    MAX_AW  = 15.0
    RIGHT_W = 2.8
    ML, MR  = 0.75, 0.30
    MT, MB  = 0.45, 0.55
    WS, HS  = 0.60, 0.45

    anim_h = (y_max - y_min) * SCALE
    anim_w = min((x_max - x_min) * SCALE, MAX_AW)
    if anim_w == MAX_AW:
        anim_h = (y_max - y_min) * anim_w / (x_max - x_min)

    rp_h  = max((anim_h - HS) / 2.0, 0.9)
    fig_w = ML + anim_w + WS + RIGHT_W + MR
    fig_h = MT + anim_h + MB

    def fx(x): return x / fig_w
    def fy(y): return y / fig_h

    fig    = plt.figure(figsize=(fig_w, fig_h))
    fig.suptitle(title, fontsize=10, y=1.0 - 0.1/fig_h)

    ax_anim  = fig.add_axes([fx(ML), fy(MB), fx(anim_w), fy(anim_h)])
    rx       = fx(ML + anim_w + WS)
    rw       = fx(RIGHT_W)
    rh       = fy(rp_h)
    ax_ener  = fig.add_axes([rx, fy(MB + rp_h + HS), rw, rh])
    ax_exits = fig.add_axes([rx, fy(MB),              rw, rh])

    # --- Animation panel static decorations ---
    ax_anim.set_xlim(x_min, x_max)
    ax_anim.set_ylim(y_min, y_max)
    ax_anim.set_xlabel("x  (column axis)", fontsize=9)
    ax_anim.set_ylabel("y", fontsize=9)

    grad_img = np.tile(
        np.clip(np.linspace(0, 1, 256), 0, 1)[np.newaxis, :], (10, 1)
    )
    ax_anim.imshow(grad_img,
                   extent=[0, grad_L, y_min, y_max],
                   origin="lower", aspect="auto",
                   cmap="Blues", alpha=0.20, zorder=0)
    ax_anim.axvline(0,      color="black",     lw=1.5, zorder=1)
    ax_anim.axvline(grad_L, color="steelblue", lw=1.5, linestyle="--", zorder=1)
    ax_anim.text(grad_L + 0.3, y_min + 0.3, f"L={grad_L:.0f}",
                 color="steelblue", fontsize=7)
    ax_anim.legend(
        handles=[mpatches.Patch(color=_POLY_COLORS[t], label=f"Polymer type {t}")
                 for t in range(4)],
        loc="upper right", fontsize=7, framealpha=0.7,
    )
    # Apply equal aspect last so that imshow's aspect="auto" and legend/vlines
    # don't override it.  adjustable="box" shrinks the axes box (not the data
    # limits) to satisfy the constraint.
    ax_anim.set_aspect("equal", adjustable="box")

    # --- Scalar panel labels ---
    ax_ener.set_xlabel("Step", fontsize=8)
    ax_ener.set_ylabel("Energy", fontsize=8)
    ax_ener.set_title("System energy", fontsize=8)
    ax_ener.tick_params(labelsize=7)

    ax_exits.set_xlabel("Step", fontsize=8)
    ax_exits.set_ylabel("Count", fontsize=8)
    ax_exits.set_title("Cumulative exited complexes", fontsize=8)
    ax_exits.tick_params(labelsize=7)

    return fig, ax_anim, ax_ener, ax_exits


def _init_artists(ax_anim, ax_ener, ax_exits):
    """
    Create all mutable artists.  Returns a dict with keys:
        scat, bond_lines, frame_text,
        ener_bg, ener_marker, exits_bg, exits_marker
    """
    scat        = ax_anim.scatter([], [], s=180, zorder=3,
                                   edgecolors="k", linewidths=0.4)
    frame_text  = ax_anim.text(0.02, 0.97, "", transform=ax_anim.transAxes,
                                fontsize=8, va="top", family="monospace")
    ener_bg,    = ax_ener.plot([],  [], color="#555555", lw=0.8, alpha=0.4)
    ener_marker,= ax_ener.plot([],  [], 'o', color="#e6194b", ms=5, zorder=5)
    exits_bg,   = ax_exits.plot([], [], color="#555555", lw=0.8, alpha=0.4)
    exits_marker,=ax_exits.plot([], [], 'o', color="#3cb44b", ms=5, zorder=5)
    return dict(
        scat=scat, bond_lines=[], frame_text=frame_text,
        ener_bg=ener_bg, ener_marker=ener_marker,
        exits_bg=exits_bg, exits_marker=exits_marker,
    )


def _update_particles(artists, ax_anim, fr):
    """Update scatter positions, backbone bond lines, and frame text for one frame."""
    pts    = fr['particles']
    xs     = [p[2] for p in pts]
    ys     = [p[3] for p in pts]
    colors = [_POLY_COLORS[p[1] % 4] for p in pts]

    scat = artists['scat']
    scat.set_offsets(np.column_stack([xs, ys]) if xs else np.empty((0, 2)))
    scat.set_color(colors)

    W = fr['W']
    for ln in artists['bond_lines']:
        ln.remove()
    artists['bond_lines'] = []
    by_id = {p[0]: p for p in pts}
    for p in pts:
        pid, _, px, py, copy_ = p
        if pid % 4 == 3:   # last segment of this polymer; pid+1 belongs to next polymer
            continue
        nxt = by_id.get(pid + 1)
        if nxt and nxt[4] == copy_:
            dx = nxt[2] - px
            dy = nxt[3] - py
            dy -= W * round(dy / W)   # minimum-image y
            if dx*dx + dy*dy < 2.5:
                ln, = ax_anim.plot([px, px + dx], [py, py + dy],
                                   color=_BACKBONE_COLOR, lw=_BOND_LW, zorder=2)
                artists['bond_lines'].append(ln)

    artists['frame_text'].set_text(
        f"step {fr['step']}  E={fr['energy']:.1f}  exited={fr['exited']}"
    )


# --------------------------------------------------------------------------- #
# Live mode
# --------------------------------------------------------------------------- #

def live_mode(path, args, grad_L, col_W):
    """
    Real-time viewer using FuncAnimation with an infinite frame generator.

    FuncAnimation owns the event loop (via plt.show), firing every `interval`
    ms.  The generator yields the current frame list when new frames are
    available, or None when the file hasn't changed.  This is more reliable
    than a manual plt.pause() loop, which doesn't guarantee a canvas repaint
    on macOS backends.
    """
    print(f"Live view — watching {path!r}.  Close the window to quit.", flush=True)

    fig, ax_anim, ax_ener, ax_exits = _build_figure(grad_L, col_W, args.title)
    artists = _init_artists(ax_anim, ax_ener, ax_exits)

    seen = [0]  # list so the closure can mutate it

    def frame_source():
        """Infinite generator: yields all_frames on new data, None otherwise."""
        while True:
            all_frames = parse_traj(path)
            if len(all_frames) > seen[0]:
                seen[0] = len(all_frames)
                yield all_frames
            else:
                yield None

    def update(all_frames):
        if all_frames is None:
            return  # nothing new; FuncAnimation redraws unchanged artists

        fr = all_frames[-1]
        ts_steps  = [f['step']   for f in all_frames]
        ts_energy = [f['energy'] for f in all_frames]
        ts_exited = [f['exited'] for f in all_frames]

        artists['ener_bg'].set_data(ts_steps, ts_energy)
        artists['exits_bg'].set_data(ts_steps, ts_exited)
        artists['ener_marker'].set_data([fr['step']], [fr['energy']])
        artists['exits_marker'].set_data([fr['step']], [fr['exited']])
        for ax in (ax_ener, ax_exits):
            ax.relim()
            ax.autoscale_view()
            if ts_steps:
                ax.set_xlim(min(ts_steps), max(ts_steps))

        _update_particles(artists, ax_anim, fr)

    # Keep a reference to _ani: if it were garbage-collected the animation stops.
    _ani = animation.FuncAnimation(  # noqa: F841
        fig, update,
        frames=frame_source(),
        interval=500,          # poll every 500 ms
        blit=False,
        cache_frame_data=False,
    )
    plt.show()


# --------------------------------------------------------------------------- #
# Batch mode
# --------------------------------------------------------------------------- #

def main():
    ap = argparse.ArgumentParser(description="Animate nucleolus trajectory.")
    ap.add_argument("traj", help="Trajectory file (output of run_nucleolus)")
    ap.add_argument("--live",   action="store_true",
                    help="Real-time mode: poll file for new frames as the simulation runs")
    ap.add_argument("--output", default=None,
                    help="Save animation (.mp4 recommended; .gif is large and slow)")
    ap.add_argument("--fps",    type=int,   default=5)
    ap.add_argument("--gradient-length", type=float, default=None,
                    help="Column length L for gradient overlay (default: read from file)")
    ap.add_argument("--width",  type=float, default=None,
                    help="Column width W in y (default: read from file)")
    ap.add_argument("--skip",   type=int,   default=1,
                    help="Render every N-th frame (batch mode only, default 1)")
    ap.add_argument("--title",  default="Nucleolus assembly simulation")
    args = ap.parse_args()

    if args.output and args.output.endswith(".gif"):
        print("Warning: GIF output can be very large and slow to render.\n"
              "         Consider --output file.mp4 instead (requires ffmpeg).",
              file=sys.stderr)

    # ---- Live mode ----
    if args.live:
        grad_L = args.gradient_length
        col_W  = args.width
        if grad_L is None or col_W is None:
            frames = parse_traj(args.traj)
            if frames:
                grad_L = grad_L or frames[0]['L']
                col_W  = col_W  or frames[0]['W']
        grad_L = grad_L or 60.0
        col_W  = col_W  or 10.0
        live_mode(args.traj, args, grad_L, col_W)
        return

    # ---- Batch mode ----
    print(f"Reading {args.traj} ...", flush=True)
    all_frames = parse_traj(args.traj)
    if not all_frames:
        print("No frames found in trajectory.", file=sys.stderr)
        sys.exit(1)
    frames = all_frames[::args.skip]
    print(f"  {len(all_frames)} frames total, rendering {len(frames)}.", flush=True)

    ts_steps  = [f['step']   for f in all_frames]
    ts_energy = [f['energy'] for f in all_frames]
    ts_exited = [f['exited'] for f in all_frames]

    grad_L = args.gradient_length or frames[0]['L']
    col_W  = args.width           or frames[0]['W']

    all_x = [p[2] for fr in frames for p in fr['particles']]
    x_max = (max(all_x) + 2.0) if all_x else (grad_L + 5)

    fig, ax_anim, ax_ener, ax_exits = _build_figure(grad_L, col_W, args.title, x_max=x_max)
    artists = _init_artists(ax_anim, ax_ener, ax_exits)

    # Draw static background time-series and fix axis ranges
    artists['ener_bg'].set_data(ts_steps, ts_energy)
    artists['exits_bg'].set_data(ts_steps, ts_exited)
    for ax in (ax_ener, ax_exits):
        ax.relim()
        ax.autoscale_view()
        ax.set_xlim(min(ts_steps), max(ts_steps))

    # Vertical cursor lines that track the current frame
    ener_vline  = ax_ener.axvline(0,  color="#e6194b", lw=0.8, ls="--", alpha=0.7)
    exits_vline = ax_exits.axvline(0, color="#3cb44b", lw=0.8, ls="--", alpha=0.7)

    def update(frame_idx):
        fr = frames[frame_idx]
        _update_particles(artists, ax_anim, fr)
        s = fr['step']
        artists['ener_marker'].set_data([s],  [fr['energy']])
        artists['exits_marker'].set_data([s], [fr['exited']])
        ener_vline.set_xdata([s, s])
        exits_vline.set_xdata([s, s])
        return ([artists['scat'], artists['frame_text'],
                 artists['ener_marker'], artists['exits_marker'],
                 ener_vline, exits_vline] + artists['bond_lines'])

    ani = animation.FuncAnimation(
        fig, update,
        frames=len(frames),
        interval=1000 // args.fps,
        blit=False,
    )

    if args.output:
        print(f"Saving to {args.output} ...", flush=True)
        writer = (animation.FFMpegWriter(fps=args.fps)
                  if args.output.endswith(".mp4")
                  else animation.PillowWriter(fps=args.fps))
        ani.save(args.output, writer=writer)
        print("Saved.")
    else:
        plt.show()


if __name__ == "__main__":
    main()
