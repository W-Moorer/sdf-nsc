#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from skimage import measure
from skimage.draw import polygon as raster_polygon


EPS = 1e-9
PLANE_AXES = {
    "xy": (0, 1),
    "xz": (0, 2),
    "yz": (1, 2),
}


def normalize(v):
    arr = np.asarray(v, dtype=float)
    n = np.linalg.norm(arr)
    if n <= 1e-12:
        return None
    return arr / n


def pick_tangent(contact):
    n = normalize(contact["normal_world"])
    u_tau = normalize(contact.get("u_tau_pred_world", [0.0, 0.0, 0.0]))
    if u_tau is not None:
        return u_tau

    v_rel = normalize(contact.get("v_rel_world", [0.0, 0.0, 0.0]))
    if v_rel is not None:
        v_t = v_rel - np.dot(v_rel, n) * n
        v_t = normalize(v_t)
        if v_t is not None:
            return v_t

    seed = np.array([1.0, 0.0, 0.0], dtype=float)
    if abs(np.dot(seed, n)) > 0.95:
        seed = np.array([0.0, 1.0, 0.0], dtype=float)
    t = seed - np.dot(seed, n) * n
    return normalize(t)


def build_local_basis(contact):
    n = normalize(contact["normal_world"])
    t = pick_tangent(contact)
    if n is None or t is None:
        raise ValueError("Failed to build local basis from contact data.")
    b = normalize(np.cross(n, t))
    if b is None:
        raise ValueError("Degenerate local basis.")
    return t, n, b


def local_hessian(contact, t, n):
    h = np.asarray(contact["hessian_world"], dtype=float)
    basis = np.column_stack([t, n])
    return basis.T @ h @ basis


def load_obj(path: Path):
    vertices = []
    faces = []
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("v "):
                _, x, y, z = line.split()[:4]
                vertices.append((float(x), float(y), float(z)))
            elif line.startswith("f "):
                refs = line.split()[1:]
                face = []
                for ref in refs:
                    token = ref.split("/")[0]
                    if token:
                        face.append(int(token) - 1)
                if len(face) >= 3:
                    faces.append(face)
    return np.asarray(vertices, dtype=float), faces


def quat_to_matrix(quat):
    w, x, y, z = quat
    return np.array(
        [
            [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
            [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
            [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
        ],
        dtype=float,
    )


def transform_vertices(vertices, ref_pos, ref_quat):
    rotation = quat_to_matrix(ref_quat)
    return vertices @ rotation.T + np.asarray(ref_pos, dtype=float)


def project(points, plane):
    ax0, ax1 = PLANE_AXES[plane]
    return points[:, [ax0, ax1]]


def triangulate_face(face):
    for i in range(1, len(face) - 1):
        yield (face[0], face[i], face[i + 1])


def plane_intersections_with_triangle(tri_points, origin, plane_normal):
    d = [float(np.dot(p - origin, plane_normal)) for p in tri_points]
    if all(abs(di) < EPS for di in d):
        return []

    hits = []
    for i, j in ((0, 1), (1, 2), (2, 0)):
        p0 = tri_points[i]
        p1 = tri_points[j]
        d0 = d[i]
        d1 = d[j]

        if abs(d0) < EPS and abs(d1) < EPS:
            continue
        if abs(d0) < EPS:
            hits.append(p0)
            continue
        if abs(d1) < EPS:
            hits.append(p1)
            continue
        if d0 * d1 > 0.0:
            continue

        alpha = d0 / (d0 - d1)
        hit = p0 + alpha * (p1 - p0)
        hits.append(hit)

    unique_hits = []
    for hit in hits:
        if not any(np.linalg.norm(hit - other) < 1e-8 for other in unique_hits):
            unique_hits.append(hit)
    if len(unique_hits) == 2:
        return [unique_hits]
    return []


def mesh_plane_segments(vertices_world, faces, origin, plane_normal):
    segments = []
    for face in faces:
        for tri in triangulate_face(face):
            tri_points = [vertices_world[idx] for idx in tri]
            segments.extend(plane_intersections_with_triangle(tri_points, origin, plane_normal))
    return segments


def project_to_local(point_w, origin_w, t, n):
    d = np.asarray(point_w, dtype=float) - np.asarray(origin_w, dtype=float)
    return np.array([np.dot(d, t), np.dot(d, n)], dtype=float)


def segments_world_to_local(segments, origin, t, n):
    out = []
    for p0, p1 in segments:
        out.append(np.vstack([project_to_local(p0, origin, t, n), project_to_local(p1, origin, t, n)]))
    return out


def draw_local_segments(ax, segments_local, color, linewidth, alpha=1.0):
    for seg in segments_local:
        ax.plot(seg[:, 0], seg[:, 1], color=color, lw=linewidth, alpha=alpha)


def rasterize_projected_mesh(points_2d, faces, xlim, ylim, raster_size):
    mask = np.zeros((raster_size, raster_size), dtype=np.uint8)
    x0, x1 = xlim
    y0, y1 = ylim
    sx = (raster_size - 1) / max(1e-12, x1 - x0)
    sy = (raster_size - 1) / max(1e-12, y1 - y0)

    def to_px(pt):
        x, y = pt
        col = (x - x0) * sx
        row = (y1 - y) * sy
        return row, col

    for face in faces:
        if len(face) < 3:
            continue
        for tri in triangulate_face(face):
            pts = np.asarray([to_px(points_2d[idx]) for idx in tri], dtype=float)
            rr, cc = raster_polygon(pts[:, 0], pts[:, 1], shape=mask.shape)
            mask[rr, cc] = 1
    return mask


def mask_to_world_contours(mask, xlim, ylim):
    contours = []
    x0, x1 = xlim
    y0, y1 = ylim
    sx = max(1e-12, x1 - x0) / (mask.shape[1] - 1)
    sy = max(1e-12, y1 - y0) / (mask.shape[0] - 1)
    for contour in measure.find_contours(mask.astype(float), 0.5):
        rows = contour[:, 0]
        cols = contour[:, 1]
        xs = x0 + cols * sx
        ys = y1 - rows * sy
        contours.append(np.column_stack([xs, ys]))
    return contours


def draw_projected_outline(ax, points_2d, faces, xlim, ylim, color, linewidth, alpha, raster_size=1200):
    mask = rasterize_projected_mesh(points_2d, faces, xlim, ylim, raster_size=raster_size)
    contours = mask_to_world_contours(mask, xlim, ylim)
    for contour in contours:
        ax.plot(contour[:, 0], contour[:, 1], color=color, lw=linewidth, alpha=alpha)


def solve_surface_eta(h2, xi_values):
    h_tt = h2[0, 0]
    h_tn = h2[0, 1]
    h_nn = h2[1, 1]

    eta = np.zeros_like(xi_values)
    for i, xi in enumerate(xi_values):
        a = 0.5 * h_nn
        b = 1.0 + h_tn * xi
        c = 0.5 * h_tt * xi * xi

        if abs(a) < 1e-12:
            eta[i] = 0.0 if abs(b) < 1e-12 else -c / b
            continue

        disc = b * b - 4.0 * a * c
        disc = max(disc, 0.0)
        r1 = (-b + np.sqrt(disc)) / (2.0 * a)
        r2 = (-b - np.sqrt(disc)) / (2.0 * a)
        eta[i] = r1 if abs(r1) < abs(r2) else r2
    return eta


def draw_gap_arrow(ax, x, y0, y1, color, label):
    ax.annotate(
        "",
        xy=(x, y1),
        xytext=(x, y0),
        arrowprops=dict(arrowstyle="<->", lw=1.1, color=color),
    )
    ax.text(x + 0.0003, 0.5 * (y0 + y1), label, fontsize=8, color=color, va="center")


def draw_global_overview(ax, frames, cam_mesh, follower_mesh, plane="xy"):
    representative_index = 0
    for i, frame in enumerate(frames):
        if frame.get("contacts"):
            representative_index = i
            break

    rep_frame = frames[representative_index]
    cam_world = transform_vertices(cam_mesh[0], rep_frame["cam"]["x_ref_W"], rep_frame["cam"]["q_WRef"])
    follower_world = transform_vertices(
        follower_mesh[0], rep_frame["follower"]["x_ref_W"], rep_frame["follower"]["q_WRef"]
    )
    cam_2d = project(cam_world, plane)
    follower_2d = project(follower_world, plane)

    marker_points = []
    marker_labels = []
    support_points = [cam_2d, follower_2d]
    for idx, frame in enumerate(frames):
        if frame.get("contacts"):
            point = np.asarray(frame["contacts"][0]["point_slave_world"], dtype=float)
        else:
            point = np.asarray(frame["follower"]["x_ref_W"], dtype=float)
        point_2d = point[list(PLANE_AXES[plane])]
        marker_points.append(point_2d)
        marker_labels.append(chr(ord("A") + idx))
        support_points.append(point_2d[None, :])

    merged = np.vstack(support_points)
    mins = merged.min(axis=0)
    maxs = merged.max(axis=0)
    size = max(maxs[0] - mins[0], maxs[1] - mins[1])
    pad = 0.10 * size
    xlim = (mins[0] - pad, maxs[0] + pad)
    ylim = (mins[1] - pad, maxs[1] + pad)

    draw_projected_outline(ax, cam_2d, cam_mesh[1], xlim, ylim, color="#475569", linewidth=1.25, alpha=0.95)
    draw_projected_outline(ax, follower_2d, follower_mesh[1], xlim, ylim, color="#d97706", linewidth=1.3, alpha=0.95)

    if marker_points:
        marker_points = np.asarray(marker_points, dtype=float)
        ax.plot(marker_points[:, 0], marker_points[:, 1], color="#c2410c", lw=1.0, ls="--", alpha=0.7)
        for point, label in zip(marker_points, marker_labels):
            circ = plt.Circle(point, radius=0.012 * size, color="#ffffff", ec="#0f172a", lw=1.0, zorder=6)
            ax.add_patch(circ)
            ax.text(point[0], point[1], label, ha="center", va="center", fontsize=8.5, zorder=7)

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_frame_on(False)
    ax.set_title(f"Global Mechanism\nRepresentative pose: {chr(ord('A') + representative_index)}", fontsize=12)


def draw_contact_frame(ax, frame, step_size, half_window, cam_mesh, follower_mesh):
    contacts = frame.get("contacts", [])
    if not contacts:
        ax.text(0.5, 0.5, "No active contact", ha="center", va="center", transform=ax.transAxes, fontsize=11)
        return

    contact = contacts[0]
    t, n, b = build_local_basis(contact)
    h2 = local_hessian(contact, t, n)

    p0 = np.asarray(contact["point_master_surface_world"], dtype=float)
    p_slave = np.asarray(contact["point_slave_world"], dtype=float)
    p_pred = p_slave + step_size * np.asarray(contact["v_rel_world"], dtype=float)

    cam_vertices_world = transform_vertices(cam_mesh[0], frame["cam"]["x_ref_W"], frame["cam"]["q_WRef"])
    follower_vertices_world = transform_vertices(follower_mesh[0], frame["follower"]["x_ref_W"], frame["follower"]["q_WRef"])

    cam_segments_local = segments_world_to_local(
        mesh_plane_segments(cam_vertices_world, cam_mesh[1], p0, b), p0, t, n
    )
    follower_segments_local = segments_world_to_local(
        mesh_plane_segments(follower_vertices_world, follower_mesh[1], p0, b), p0, t, n
    )

    slave_local = project_to_local(p_slave, p0, t, n)
    pred_local = project_to_local(p_pred, p0, t, n)

    xi = np.linspace(-half_window, half_window, 401)
    eta_curve = solve_surface_eta(h2, xi)

    draw_local_segments(ax, cam_segments_local, color="#334155", linewidth=1.2, alpha=0.95)
    draw_local_segments(ax, follower_segments_local, color="#d97706", linewidth=1.4, alpha=0.95)
    ax.plot(xi, np.zeros_like(xi), color="#94a3b8", lw=1.2, ls="--", label="1st-order tangent")
    ax.plot(xi, eta_curve, color="#1d4ed8", lw=1.8, label="2nd-order Hessian manifold")

    ax.scatter([0.0], [0.0], s=20, color="#0f172a", zorder=5)
    ax.scatter([slave_local[0]], [slave_local[1]], s=24, color="#b91c1c", zorder=6)
    ax.scatter([pred_local[0]], [pred_local[1]], s=22, color="#ea580c", marker="s", zorder=6)

    n_len = half_window * 0.55
    ax.annotate(
        "",
        xy=(n_len * 0.12, n_len),
        xytext=(0.0, 0.0),
        arrowprops=dict(arrowstyle="->", lw=1.5, color="#0f766e"),
    )
    ax.text(n_len * 0.15, n_len * 0.95, "$n$", color="#0f766e", fontsize=9)

    v_rel_local = pred_local - slave_local
    ax.annotate(
        "",
        xy=(slave_local[0] + v_rel_local[0], slave_local[1] + v_rel_local[1]),
        xytext=(slave_local[0], slave_local[1]),
        arrowprops=dict(arrowstyle="->", lw=1.5, color="#c2410c"),
    )
    ax.text(
        slave_local[0] + v_rel_local[0] * 0.55,
        slave_local[1] + v_rel_local[1] * 0.55 + 0.0004,
        "$v_{rel}\\Delta t$",
        color="#c2410c",
        fontsize=8,
    )

    draw_gap_arrow(ax, slave_local[0], 0.0, slave_local[1], "#7c3aed", "$\\phi$")

    pred_curve_eta = solve_surface_eta(h2, np.array([pred_local[0]]))[0]
    draw_gap_arrow(ax, pred_local[0], pred_curve_eta, pred_local[1], "#0f766e", "$\\phi_{2nd}$")

    text = (
        f"$\\phi$ = {contact['phi']:.3e}\n"
        f"$C_{{curve}}$ = {contact['curvature_term']:.3e}\n"
        f"$\\phi_{{eff}}$ = {contact['phi_eff']:.3e}\n"
        f"$t^T H t$ = {h2[0,0]:.3e}"
    )
    ax.text(
        0.03,
        0.04,
        text,
        transform=ax.transAxes,
        fontsize=8,
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#cbd5e1", alpha=0.96),
    )

    ax.set_xlim(-half_window, half_window)
    ax.set_ylim(-half_window * 0.8, half_window * 0.8)


def main():
    parser = argparse.ArgumentParser(
        description="Render local contact diagrams using the paper's Hessian-SDF contact quantities."
    )
    parser.add_argument("--input", required=True, help="Snapshot JSON exported by baseline_cam_nsc")
    parser.add_argument("--output", required=True, help="Output figure path (.pdf or .svg recommended)")
    parser.add_argument("--window", type=float, default=0.012, help="Half window size in local tangent-normal units")
    parser.add_argument("--global-plane", default="xy", choices=sorted(PLANE_AXES.keys()), help="Projection plane for the global mechanism panel")
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with input_path.open("r", encoding="utf-8") as f:
        snapshot = json.load(f)

    frames = snapshot.get("frames", [])
    if not frames:
        raise SystemExit("No frames found in snapshot JSON.")

    step_size = float(snapshot.get("step_size", 0.0))
    cam_mesh = load_obj(Path(snapshot["cam_mesh_path"]))
    follower_mesh = load_obj(Path(snapshot["follower_mesh_path"]))

    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = ["Times New Roman", "Times", "DejaVu Serif"]

    width_ratios = [1.45] + [1.0] * len(frames)
    fig = plt.figure(figsize=(3.4 + 3.15 * len(frames), 3.5), constrained_layout=True)
    gs = fig.add_gridspec(1, len(frames) + 1, width_ratios=width_ratios)
    overview_ax = fig.add_subplot(gs[0, 0])
    axes = [fig.add_subplot(gs[0, i + 1]) for i in range(len(frames))]

    draw_global_overview(overview_ax, frames, cam_mesh, follower_mesh, plane=args.global_plane)

    for idx, (ax, frame) in enumerate(zip(axes, frames)):
        draw_contact_frame(ax, frame, step_size, args.window, cam_mesh, follower_mesh)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_frame_on(False)
        ax.set_title(f"{chr(ord('A') + idx)}   t = {frame['captured_time']:.3f} s", fontsize=11)

    legend_handles = [
        Line2D([0], [0], color="#475569", lw=1.1, label="Cam mesh outline"),
        Line2D([0], [0], color="#d97706", lw=1.2, label="Follower mesh outline"),
        Line2D([0], [0], color="#c2410c", lw=1.0, ls="--", label="Snapshot path"),
        Line2D([0], [0], color="#334155", lw=1.2, label="Cam mesh slice"),
        Line2D([0], [0], color="#d97706", lw=1.4, label="Follower mesh slice"),
        Line2D([0], [0], color="#94a3b8", lw=1.2, ls="--", label="1st-order tangent model"),
        Line2D([0], [0], color="#1d4ed8", lw=1.8, label="2nd-order Hessian model"),
        Line2D([0], [0], marker="o", color="#0f172a", lw=0, label="Master surface point"),
        Line2D([0], [0], marker="o", color="#b91c1c", lw=0, label="Slave sample"),
        Line2D([0], [0], marker="s", color="#ea580c", lw=0, label="Predicted point"),
    ]
    fig.legend(handles=legend_handles, loc="upper center", ncol=5, frameon=False, bbox_to_anchor=(0.53, 1.04))

    fig.savefig(output_path)
    print(f"Saved figure to: {output_path}")


if __name__ == "__main__":
    main()
