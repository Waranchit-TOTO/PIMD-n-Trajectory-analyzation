from pathlib import Path

#print("This code separate the beads trajectory into different files based on the bead number. It reads the input trajectory file names coor.xyz, identifies the bead number, and writes each bead number's trajectory to a separate output file.")

base_path = Path("./")
data_cases = ["folder1", "folder2", "folderN"] # set your folder name that containing your trajectory here

# ======= Manual input for describe the system =======
# Number of beads in the system
n_beads = 60
limit_steps = 100000 # Real number of steps in the simulation that want to analyze 
skip_steps = 20  # Real number of logged steps in the simulation (must the same with actual simulation) 
n_steps = limit_steps // skip_steps


# ======== Function to get the number of atoms per system =======
def get_atoms_per_system(path_to_coor: Path, n_beads: int) -> int:
	with path_to_coor.open() as file_handle:
		total_atoms_first_step = int(file_handle.readline().strip())
	if total_atoms_first_step % n_beads != 0:
		raise ValueError(
			f"First step in {path_to_coor} has {total_atoms_first_step} atoms, "
			f"which is not divisible by {n_beads} beads"
		)

	return total_atoms_first_step // n_beads

info = {}
for case_name in data_cases:
	path_to_coor = base_path / f"{case_name}/2-PIMD/Result/coor.xyz" # Set your path pattern to your trajectory file 
	n_atoms = get_atoms_per_system(path_to_coor, n_beads)

	print(f"{case_name}: {n_atoms} atoms per system")
	info[case_name] = {
		"n_atoms": n_atoms,
	}
#print(info)

# ======== Function to separate the trajectory into different files based on bead number =======
def separate_trajectory_by_bead(path_to_coor: Path, n_beads: int, n_atoms: int, n_steps: int) -> None:
    with path_to_coor.open() as file_handle:
        for step in range(n_steps):
            # Read the number of atoms for this step
            total_atoms = int(file_handle.readline().strip())
            if total_atoms != n_beads * n_atoms:
                raise ValueError(
                    f"Step {step} in {path_to_coor} has {total_atoms} atoms, "
                    f"expected {n_beads * n_atoms} atoms"
                )

            # Read the comment line (usually contains step information)
            comment_line = file_handle.readline().strip()

            # Read the coordinates for all atoms
            coordinates = [file_handle.readline().strip() for _ in range(total_atoms)]

            # Write coordinates to separate files based on bead number
            for bead_index in range(n_beads):
                start_index = bead_index * n_atoms
                end_index = start_index + n_atoms
                bead_coordinates = coordinates[start_index:end_index]

                output_file_path = bead_folder_path / f"bead_{bead_index + 1}.xyz"
                with output_file_path.open("a") as output_file:
                    output_file.write(f"{n_atoms}\n")
                    output_file.write(f"{comment_line}\n")
                    output_file.write("\n".join(bead_coordinates) + "\n")


# Separate the trajectory for each case
for case_name in data_cases:
    path_to_coor = base_path / f"{case_name}/2-PIMD/Result/coor.xyz"
    n_atoms = info[case_name]["n_atoms"]
    bead_folder_path = path_to_coor.parent / "bead-trajectories"
    bead_folder_path.mkdir(exist_ok=True)
    separate_trajectory_by_bead(path_to_coor, n_beads, n_atoms, n_steps)
