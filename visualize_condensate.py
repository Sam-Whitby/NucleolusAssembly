#!/usr/bin/env python3
"""
visualize_condensate.py

Reads a condensate trajectory file produced by run_condensate and generates
a 3-panel animated figure:
  Left:         particle positions inside the circular condensate, coloured
                by polymer type, with radial gradient overlay and labels.
  Top-right:    energy vs simulation step.
  Bottom-right: cumulative recycled components vs simulation step.

Usage
-----
    python3 visualize_condensate.py <traj_file> [OPTIONS]

Options
    --output FILE       Save animation to FILE (.mp4 or .gif). Default: display.
    --fps N             Frames per second (default 5).
    --radius R          Override condensate radius from file.
    --skip N            Render only every N-th frame (default 1).
    --title TEXT        Figure title.

Trajectory format (written by run_condensate)
---------------------------------------------
    <N_particles>
    step=S energy=E exited=X R_c=RC cx=CX cy=CY nCopies=C coupling=MODE
    <id> <poly_type> <x> <y> <copy>
    ...   (repeated N_particles times per frame)
"""

import argparse
import re
import sys

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.animation as animation
import numpy as np


_POLY_COLORS  = ["tab:blue", "tab:orange", "tab:green", "tab:red"]
_BACKBONE_COLOR = "#333333"
_BOND_LW = 1.5


def _kv(header, key, default=None):
    m = re.search(rf'{key}=([^\s]+)', header)
    return m.group(1) if m else default


def parse_traj(path):
    frames = []
    with open(path) as fh:
        lines = fh.readlines()

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

        step    = int(_kv(hdr, 'step', 0))
        energy  = float(_kv(hdr, 'energy', 0.0))
        exited  = int(_kv(hdr, 'exited', 0))
        R_c     = float(_kv(hdr, 'R_c', 60))
        cx      = float(_kv(hdr, 'cx', R_c))
        cy      = float(_kv(hdr, 'cy', R_c))
        nCopies = int(_kv(hdr, 'nCopies', 4))
        coupling = _kv(hdr, 'coupling', 'product')

        particles = []
        for _ in range(n):
            if i >= len(lines):
                break
            parts = lines[i].strip().split()
            i += 1
            if len(parts) >= 5:
                particles.append((
                    int(parts[0]),   # id
                    int(parts[1]),   # poly_type
                    float(parts[2]), # x  (absolute lattice coords)
                    float(parts[3]), # y
                    int(parts[4]),   # copy
                ))

        frames.append(dict(
            step=step, energy=energy, exited=exited,
            R_c=R_c, cx=cx, cy=cy, nCopies=nCopies,
            coupling=coupling, particles=particles,
        ))

    return frames


def main():
    ap = argparse.ArgumentParser(description="Animate circular condensate trajectory.")
    ap.add_argument("traj", help="Trajectory file (output of run_condensate)")
    ap.add_argument("--output", default=None)
    ap.add_argument("--fps",    type=int,   default=5)
    ap.add_argument("--radius", type=float, default=None,
                    help="Override condensate radius R_c")
    ap.add_argument("--skip",   type=int,   default=1)
    ap.add_argument("--title",  default="Circular condensate assembly")
    args = ap.parse_args()

    print(f"Reading {args.traj} ...", flush=True)
    all_frames = parse_traj(args.traj)
    if not all_frames:
        print("No frames found.", file=sys.stderr)
        sys.exit(1)
    frames = all_frames[::args.skip]
    print(f"  {len(all_frames)} total frames, rendering {len(frames)}.", flush=True)

    R_c = args.radius if args.radius else frames[0]['R_c']
    cx  = frames[0]['cx']
    cy  = frames[0]['cy']

    ts_steps  = [f['step']   for f in all_frames]
    ts_energy = [f['energy'] for f in all_frames]
    ts_exited = [f['exited'] for f in all_frames]

    # Convert absolute lattice positions to positions relative to centre.
    def rel_xy(fr):
        return [(p[2] - cx, p[3] - cy) for p in fr['particles']]

    # ---------------------------------------------------------------------- #
    # Figure layout
    # ---------------------------------------------------------------------- #
    SCALE    = 0.045   # inches per lattice unit
    MAX_AW   = 8.0
    RIGHT_W  = 2.8
    ML, MR   = 0.60, 0.30
    MT, MB   = 0.50, 0.55
    WS, HS   = 0.55, 0.45

    diam   = 2.0 * R_c
    anim_h = min(diam * SCALE, MAX_AW)
    anim_w = anim_h   # square for the circle
    if anim_w > MAX_AW:
        anim_w = MAX_AW
        anim_h = MAX_AW

    rp_h   = max((anim_h - HS) / 2.0, 0.9)
    fig_w  = ML + anim_w + WS + RIGHT_W + MR
    fig_h  = MT + anim_h + MB

    def fx(x): return x / fig_w
    def fy(y): return y / fig_h

    fig = plt.figure(figsize=(fig_w, fig_h))
    coupling_label = frames[0].get('coupling', 'product')
    fig.suptitle(f"{args.title}  [coupling={coupling_label}]",
                 fontsize=10, y=1.0 - 0.1/fig_h)

    ax_anim  = fig.add_axes([fx(ML),               fy(MB),              fx(anim_w), fy(anim_h)])
    rx       = fx(ML + anim_w + WS)
    rw       = fx(RIGHT_W)
    rh       = fy(rp_h)
    ax_ener  = fig.add_axes([rx, fy(MB + rp_h + HS), rw, rh])
    ax_exits = fig.add_axes([rx, fy(MB),               rw, rh])

    # --- Animation panel ---
    pad = R_c * 0.15
    ax_anim.set_xlim(-R_c - pad, R_c + pad)
    ax_anim.set_ylim(-R_c - pad, R_c + pad)
    ax_anim.set_aspect('equal')
    ax_anim.set_xlabel("x − cx  (lattice units)", fontsize=9)
    ax_anim.set_ylabel("y − cy  (lattice units)", fontsize=9)

    # Radial gradient overlay (blue = denaturing at centre, white at edge)
    theta = np.linspace(0, 2 * np.pi, 360)
    r_vals = np.linspace(0, R_c, 200)
    T, R   = np.meshgrid(theta, r_vals)
    X_grid = R * np.cos(T)
    Y_grid = R * np.sin(T)
    C_grid = R / R_c   # γ value
    ax_anim.pcolormesh(X_grid, Y_grid, C_grid,
                       cmap='Blues_r', alpha=0.20, shading='auto', zorder=0)

    # Condensate boundary circle
    circle_patch = mpatches.Circle((0, 0), R_c,
                                    fill=False, edgecolor='steelblue',
                                    linewidth=1.5, linestyle='--', zorder=2)
    ax_anim.add_patch(circle_patch)
    ax_anim.text(R_c + 0.5, 0, f"R_c={R_c:.0f}",
                 color='steelblue', fontsize=7, va='center')

    # Injection zone marker
    inj_circle = mpatches.Circle((0, 0), 1.5,
                                  fill=False, edgecolor='red',
                                  linewidth=1.0, linestyle=':', zorder=3)
    ax_anim.add_patch(inj_circle)
    ax_anim.plot(0, 0, 'x', color='red', ms=5, zorder=4)  # centre hard wall

    # Legend
    legend_patches = [
        mpatches.Patch(color=_POLY_COLORS[t], label=f"Polymer type {t}")
        for t in range(4)
    ]
    ax_anim.legend(handles=legend_patches, loc='upper right', fontsize=7,
                   framealpha=0.7)

    scat = ax_anim.scatter([], [], s=140, zorder=5,
                            edgecolors='k', linewidths=0.4)
    bond_lines = []
    frame_text = ax_anim.text(0.02, 0.97, "", transform=ax_anim.transAxes,
                               fontsize=8, va='top', family='monospace')

    # --- Scalar panels ---
    ax_ener.plot(ts_steps, ts_energy, color='#555555', lw=0.8, alpha=0.4)
    ener_marker, = ax_ener.plot([], [], 'o', color='#e6194b', ms=5, zorder=5)
    ener_vline   = ax_ener.axvline(0, color='#e6194b', lw=0.8, ls='--', alpha=0.7)
    ax_ener.set_xlabel("Step", fontsize=8)
    ax_ener.set_ylabel("Energy", fontsize=8)
    ax_ener.set_title("System energy (excl. injection zone)", fontsize=8)
    ax_ener.tick_params(labelsize=7)
    ax_ener.set_xlim(min(ts_steps), max(ts_steps))

    ax_exits.plot(ts_steps, ts_exited, color='#555555', lw=0.8, alpha=0.4)
    exits_marker, = ax_exits.plot([], [], 'o', color='#3cb44b', ms=5, zorder=5)
    exits_vline   = ax_exits.axvline(0, color='#3cb44b', lw=0.8, ls='--', alpha=0.7)
    ax_exits.set_xlabel("Step", fontsize=8)
    ax_exits.set_ylabel("Count", fontsize=8)
    ax_exits.set_title("Cumulative recycled components", fontsize=8)
    ax_exits.tick_params(labelsize=7)
    ax_exits.set_xlim(min(ts_steps), max(ts_steps))

    # ---------------------------------------------------------------------- #
    # Update function
    # ---------------------------------------------------------------------- #
    def update(frame_idx):
        nonlocal bond_lines
        fr    = frames[frame_idx]
        pts   = fr['particles']
        rcx   = fr['cx']
        rcy   = fr['cy']

        xs     = [p[2] - rcx for p in pts]
        ys     = [p[3] - rcy for p in pts]
        colors = [_POLY_COLORS[p[1] % 4] for p in pts]

        scat.set_offsets(np.column_stack([xs, ys]))
        scat.set_color(colors)

        # Backbone bonds: consecutive ids within same copy and copy, dist² ≤ 2.5
        for ln in bond_lines:
            ln.remove()
        bond_lines = []
        by_id = {p[0]: p for p in pts}
        for p in pts:
            pid, _, px, py, copy_ = p
            nxt = by_id.get(pid + 1)
            if nxt and nxt[4] == copy_:
                dist2 = (nxt[2] - px)**2 + (nxt[3] - py)**2
                if dist2 < 2.5:
                    rx1, ry1 = px - rcx, py - rcy
                    rx2, ry2 = nxt[2] - rcx, nxt[3] - rcy
                    ln, = ax_anim.plot([rx1, rx2], [ry1, ry2],
                                       color=_BACKBONE_COLOR, lw=_BOND_LW, zorder=4)
                    bond_lines.append(ln)

        frame_text.set_text(
            f"step {fr['step']}  E={fr['energy']:.1f}  recycled={fr['exited']}"
        )

        s = fr['step']
        ener_marker.set_data([s], [fr['energy']])
        ener_vline.set_xdata([s, s])
        exits_marker.set_data([s], [fr['exited']])
        exits_vline.set_xdata([s, s])

        return ([scat, frame_text, ener_marker, ener_vline,
                 exits_marker, exits_vline] + bond_lines)

    ani = animation.FuncAnimation(
        fig, update,
        frames=len(frames),
        interval=1000 // args.fps,
        blit=False,
    )

    if args.output:
        print(f"Saving to {args.output} ...", flush=True)
        if args.output.endswith(".mp4"):
            writer = animation.FFMpegWriter(fps=args.fps)
        else:
            writer = animation.PillowWriter(fps=args.fps)
        ani.save(args.output, writer=writer)
        print("Saved.")
    else:
        plt.show()


if __name__ == "__main__":
    main()
