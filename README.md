# PIMD-n-Trajectory-analyzation  

Mainly focusing on H transfer and PIMD-CLMD trajectory analyzing  

The concept idea of "c1-digest-beads-tranjectory.py"  
    -The trajectory that work with this code is coordinate file trajectory from PIMD calculation. XYZ format with repeating coordinate with respect to number of beads in each trajectory step.  
    -Code will digest each bead by finding the total atomic number of the cluster or system by "the total atomic number in first line of XYZ devided to number of beads (need to input by User)"  
    -Then the trajectory will be divided and saving in the new folder names "bead-trajectories" in XYZ format.  
    
The concept idea of "c2-H_transfer_detection.py"  
    -The trjectory of each beads obtained from "c1-digest-beads-tranjectory.py" will be read and calculate the "Delta_R" of specific three atoms where the center atom was consider the transfer between those two atoms. In this code I noted N, H, and O three atom to consider H moving between N and O atom. Where you can adopt to your case by think that are X, Y, Z atoms index. Thus the Delta_R will calculated by Delta_R = R_YZ - R_XY.  
  -This code will give the results of time series as Delta_R plot with respect to the step of the trajectory and also give the distribution of those Delta_R with respect to distance as x-axis which giving the detail of probability distribution over the beads and time. Where every analysis will also plot CLMD as the reference or information data to be analyzed. (You can copy the code and edit it to the direction you need, Happy to shaere! :D)  
  -This code also plot the distribution of specific bond length like R_OH to compare between CLMD and PIMD. (Let's change it the bond you need)
