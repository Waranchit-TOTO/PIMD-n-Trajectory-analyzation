# PIMD-n-Trajectory-analyzation
Mainly focusing on H transfer and PIMD-CLMD trajectory analyzing

The concept idea of "digest-beads-tranjectory.py"
    -The trajectory that work with this code is coordinate file trajectory from PIMD calculation. XYZ format with repeating coordinate with respect to number of beads in each trajectory step.
    -Code will digest each bead by finding the total atomic number of the cluster or system by "the total atomic number in first line of XYZ devided to number of beads (need to input by User)"
    -Then the trajectory will be divided and saving in the new folder names "bead-trajectories" in XYZ format.