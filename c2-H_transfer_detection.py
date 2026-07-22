from pathlib import Path
import numpy as np
import time
import matplotlib.pyplot as plt
from ase.io import iread

# this code is used to detect the transfer of atoms between beads in a PIMD simulation. It reads the trajectory files, separates the coordinates based on bead number, and analyzes the transfer events.
# Additionally, it plots the distribution of distances between specified atoms and the delta distance (dist_HY - dist_XH) across all beads to visualize potential transfer events.
# Additionally, for reference distance including the delta distance will using the CLMD results to compare the distance distribution between the PIMD simulation and the classical MD simulation. The CLMD results are stored in the "1-CLMD" folder, and the distance distribution is calculated from the CLMD trajectory files.

base_path = Path("./")
output_root = base_path / "analysis-output"
show_plots = False

plt.rcParams.update({
    "axes.titlesize": 16,
    "axes.labelsize": 16,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "legend.fontsize": 16,
})

# ======= Manual input for describe the system =======
# Number of beads in the system
n_beads = 60
limit_steps = 100000 # The real number of steps in the simulation that want to analyze 
eq_step = 10000 # The real step number to start analyzing the trajectory after equilibration (must be less than limit_steps)
skip_steps = 20  # Real number of logged steps in the simulation (must the same with actual simulation) 
last_step = limit_steps // skip_steps
start_step = eq_step // skip_steps

# Mannual input the atom index to detect the transfer of atoms between beads
# the X-H-Y, where the X and Y are atom to consider the transfer of atoms between beads, and H is the atom to detect the transfer. The atom index is 0-based, e.g., 0 means the first atom in the system which is inline with the input file coor.xyz. The atom index is used to detect the transfer of atoms between beads in the PIMD simulation.
data_cases = {"1-m1n0":[0,8,6] , "2-m1n1":[0,7,5]} # N-H-O
verbose_events = False

def distance_r(a, b):
    """Calculate the distance between two points a and b in 3D space."""
    return np.sqrt(np.sum((a - b) ** 2))

def detect_transfer(path_to_beads_coor: Path, atom_indices: list, start_step: int, end_step: int) -> list:
    """Detect the transfer of atoms between beads in a PIMD simulation through distance differences."""
    delta = []
    NO_distances = []
    NH_distances = []
    HO_distances = []
    NO_distances_full = []
    for step, coordinates in enumerate(iread(path_to_beads_coor, index=f"{start_step}:{end_step}"), start=start_step):
        # Stream only the required frame range from disk once.
        # Get the positions of the atoms of interest
        pos_N = coordinates[atom_indices[0]].position
        pos_H = coordinates[atom_indices[1]].position
        pos_O = coordinates[atom_indices[2]].position

        # Calculate distances
        dist_NO = distance_r(pos_N, pos_O)
        dist_NH = distance_r(pos_N, pos_H)
        dist_HO = distance_r(pos_H, pos_O)
        

        # Measure the distance difference for transfer detection
        delta.append(dist_HO - dist_NH)
        NO_distances.append(dist_NO)
        NH_distances.append(dist_NH)
        HO_distances.append(dist_HO)
        NO_distances_full.append(dist_NO)
    return delta, NO_distances, NH_distances, HO_distances, NO_distances_full

def gaussian_smooth(data, sigma=1):
    """Smooth the data using a Gaussian filter."""
    radius = max(1, int(3 * sigma))
    x = np.arange(-radius, radius + 1)
    kernel = np.exp(-(x ** 2) / (2 * sigma ** 2))
    kernel /= kernel.sum()
    return np.convolve(data, kernel, mode="same")

for key in data_cases.keys():
    delta_all_beads = []
    NO_distances_all_beads = []
    HO_distances_all_beads = []
    case_start = time.perf_counter()
    case_output_dir = output_root / key
    case_output_dir.mkdir(parents=True, exist_ok=True)
    for i in range(1, n_beads + 1):
        path_to_beads_coor = base_path / f"{key}/2-PIMD/Result/bead-trajectories/bead_{i}.xyz"
        delta, NO_distances, NH_distances, HO_distances, NO_distances_full = detect_transfer(path_to_beads_coor, data_cases[key], start_step, last_step)
        delta_all_beads.append(np.asarray(delta))
        NO_distances_all_beads.append(np.asarray(NO_distances))
        HO_distances_all_beads.append(np.asarray(HO_distances))
    
    # refence distance distribution from CLMD simulation
    path_to_clmd_coor = base_path / f"{key}/1-CLMD/Result/coor.xyz"
    delta_clmd, NO_distances_clmd, NH_distances_clmd, HO_distances_clmd, NO_distances_full_clmd = detect_transfer(path_to_clmd_coor, data_cases[key], start_step, last_step)
    delta_clmd = np.asarray(delta_clmd)
    NO_distances_clmd = np.asarray(NO_distances_clmd)
    HO_distances_clmd = np.asarray(HO_distances_clmd)

    # Plot distribution of every bead on the same axes.
    all_values = np.concatenate(delta_all_beads)
    bins = 700 #np.linspace(all_values.min(), all_values.max(), 160)

    plt.figure(figsize=(10, 6))
    for bead_index, NO_distances in enumerate(NO_distances_all_beads, start=1):
        hist, edges = np.histogram(NO_distances, bins=bins, density=True)
        centers = 0.5 * (edges[:-1] + edges[1:])
        smoothed_hist = gaussian_smooth(hist, sigma=1.5)
        plt.plot(centers, smoothed_hist, alpha=0.18, linewidth=1.0)
    # Add average NO distance distribution across all beads and compare with CLMD results
    all_NO_distances = np.concatenate(NO_distances_all_beads)
    hist, edges = np.histogram(all_NO_distances, bins=bins, density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    smoothed_hist = gaussian_smooth(hist, sigma=1.5)
    plt.plot(centers, smoothed_hist, color="black", alpha=0.5, linewidth=1.5, label="Average-PIMD")
    # Compare with CLMD results
    hist, edges = np.histogram(NO_distances_clmd, bins=bins, density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    smoothed_hist = gaussian_smooth(hist, sigma=1.5)
    plt.plot(centers, smoothed_hist, color="red", alpha=0.5, linewidth=1.5, label="CLMD")

    plt.title(f"N-O Distance Distribution Across Beads - Case {key}")
    plt.xlabel("Distance between N and O (Å)")
    plt.ylabel("Probability Density")
    plt.ylim(0, )
    plt.legend()
    plt.tight_layout()
    plt.savefig(case_output_dir / "NO_distance_distribution.png", dpi=300)
    if show_plots:
        plt.show()
    plt.close()

# Let's analyze the average of XY distance vs step for each bead to see if there are any significant changes that might indicate a transfer event. We can plot the XY distance over time for each bead to visualize any potential transfer events.
    plt.figure(figsize=(10, 6))
    steps = np.arange(start_step, last_step)
    pimd_no_values = np.concatenate(NO_distances_all_beads)
    pimd_xy_steps = np.tile(steps, n_beads)

    # Using all PIMD bead points as scatter cloud and CLMD points for comparison.
    plt.axhline(np.mean(NO_distances_clmd), color="red", linestyle="--", linewidth=1.0, label="Average-CLMD")
    plt.axhline(np.mean(pimd_no_values), color="black", linestyle="--", linewidth=1.0, label="Average-PIMD")
    plt.xlim(start_step, last_step)
    #plt.ylim(0, )

    plt.title(f"N-O Distance Over Time (All PIMD Beads) - Case {key}")
    plt.xlabel("Step")
    plt.ylabel("Distance between N and O (Å)")
    plt.scatter(pimd_xy_steps, pimd_no_values, alpha=0.5, s=10, label="PIMD (all beads)")
    plt.scatter(steps, NO_distances_clmd, alpha=0.5, s=10, label="CLMD")
    plt.legend()
    plt.tight_layout() 
    plt.savefig(case_output_dir / "NO_distance_time_series.png", dpi=300)
    if show_plots:
        plt.show()
    plt.close()

# Delta R time-series comparison (PIMD vs CLMD), similar to XY time-series plot.
    plt.figure(figsize=(10, 6))
    pimd_delta_values = np.concatenate(delta_all_beads)
    pimd_delta_steps = np.tile(steps, n_beads)

    plt.axhline(np.mean(delta_clmd), color="red", linestyle="--", linewidth=1.0, label="CLMD")
    plt.axhline(np.mean(pimd_delta_values), color="black", linestyle="--", linewidth=1.0, label="Average-PIMD")
    plt.xlim(start_step, last_step)

    plt.title(f"δR Over Time (All PIMD Beads) - Case {key}")
    plt.xlabel("Step")
    plt.ylabel("δR = R_HO - R_NH (Å)")
    plt.scatter(pimd_delta_steps, pimd_delta_values, alpha=0.5, s=10, label="PIMD (all beads)")
    plt.scatter(steps, delta_clmd, alpha=0.5, s=10, label="CLMD")
    plt.legend()
    plt.tight_layout()
    plt.savefig(case_output_dir / "delta_time_series.png", dpi=300)
    if show_plots:
        plt.show()
    plt.close()

# ====== Delta R distribution across beads ======
    plt.figure(figsize=(10, 6))
    for bead_index, delta in enumerate(delta_all_beads, start=1):
        hist, edges = np.histogram(delta, bins=bins, density=True)
        centers = 0.5 * (edges[:-1] + edges[1:])
        smoothed_hist = gaussian_smooth(hist, sigma=1.5)
        plt.plot(centers, smoothed_hist, alpha=0.18, linewidth=1.0)
    
    # Add average delta distribution across all beads and compare with CLMD results
    all_delta = np.concatenate(delta_all_beads)
    hist, edges = np.histogram(all_delta, bins=bins, density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    smoothed_hist = gaussian_smooth(hist, sigma=1.5)
    plt.plot(centers, smoothed_hist, color="black", alpha=0.5, linewidth=1.5, label="Average Delta Distribution")
    
    # Compare with CLMD results
    hist, edges = np.histogram(delta_clmd, bins=bins, density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    smoothed_hist = gaussian_smooth(hist, sigma=1.5)
    plt.plot(centers, smoothed_hist, color="red", alpha=0.5, linewidth=1.5, label="CLMD Delta Distribution")

#    plt.axvline(0.117, color="black", linestyle="--", linewidth=1.0, label="Delta = 0")
    plt.title(f"δR Distribution Across Beads - Case {key}")
    plt.xlabel("δR = R_HO - R_NH (Å)")
    plt.ylabel("Probability Density")
    plt.ylim(0, )
    plt.legend()
    plt.tight_layout()
    plt.savefig(case_output_dir / "delta_distribution.png", dpi=300)
    if show_plots:
        plt.show()
    plt.close()

    # To measure the H moving away from the bonded Y atom, plot the dist_HY distribution across all beads and compare with CLMD results.
    plt.figure(figsize=(10, 6))
    all_dist_HO = np.concatenate(HO_distances_all_beads)
    for bead_index, dist_HO in enumerate(HO_distances_all_beads, start=1):
        hist, edges = np.histogram(dist_HO, bins=bins, density=True)
        centers = 0.5 * (edges[:-1] + edges[1:])
        smoothed_hist = gaussian_smooth(hist, sigma=1.5)
        plt.plot(centers, smoothed_hist, alpha=0.18, linewidth=1.0)
    hist, edges = np.histogram(all_dist_HO, bins=bins, density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    smoothed_hist = gaussian_smooth(hist, sigma=1.5)
    plt.plot(centers, smoothed_hist, color="black", alpha=0.6, linewidth=1.5, label="Average-PIMD")
    # Compare with CLMD results
    hist, edges = np.histogram(HO_distances_clmd, bins=bins, density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    smoothed_hist = gaussian_smooth(hist, sigma=1.5)
    plt.plot(centers, smoothed_hist, color="red", alpha=0.5, linewidth=1.5, label="Average-CLMD")
    plt.title(f"R_HO Distribution Across Beads - Case {key}")
    plt.xlabel("R_HO (Å)")
    plt.ylabel("Probability Density")
    plt.ylim(0, )
    plt.legend()
    plt.tight_layout()
    plt.savefig(case_output_dir / "dist_HO_distribution.png", dpi=300)
    if show_plots:
        plt.show()
    plt.close()

    pimd_xy_mean = float(np.mean(all_NO_distances))
    pimd_xy_std = float(np.std(all_NO_distances))
    clmd_xy_mean = float(np.mean(NO_distances_clmd))
    clmd_xy_std = float(np.std(NO_distances_clmd))
    pimd_delta_mean = float(np.mean(all_delta))
    pimd_delta_std = float(np.std(all_delta))
    clmd_delta_mean = float(np.mean(delta_clmd))
    clmd_delta_std = float(np.std(delta_clmd))

    summary_lines = [
        f"Case: {key}",
        f"Beads analyzed: {n_beads}",
        f"Step window: {start_step} to {last_step - 1}",
        "",
        "XY distance statistics (A)",
        f"PIMD mean: {pimd_xy_mean:.6f}",
        f"PIMD std: {pimd_xy_std:.6f}",
        f"CLMD mean: {clmd_xy_mean:.6f}",
        f"CLMD std: {clmd_xy_std:.6f}",
        f"Mean shift (PIMD - CLMD): {pimd_xy_mean - clmd_xy_mean:.6f}",
        "",
        "Delta = dist_HY - dist_XH statistics (A)",
        f"PIMD mean: {pimd_delta_mean:.6f}",
        f"PIMD std: {pimd_delta_std:.6f}",
        f"CLMD mean: {clmd_delta_mean:.6f}",
        f"CLMD std: {clmd_delta_std:.6f}",
        f"Mean shift (PIMD - CLMD): {pimd_delta_mean - clmd_delta_mean:.6f}",
    ]
    (case_output_dir / "interpretation_summary.txt").write_text("\n".join(summary_lines) + "\n")

    case_time = time.perf_counter() - case_start
    print(f"{key}: processed {n_beads} beads in {case_time:.2f} s")
    print(f"saved results to {case_output_dir}")




