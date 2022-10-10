import mdtraj as md
import numpy as np
import sys
import os

##Settings for input trajectories and output trajectory

#Define number of subunits in oligomer 1 and 2
oligomersize1 = int(sys.argv[1])
oligomersize2 = int(sys.argv[2])

#Length of output oligomer
target_length = int(sys.argv[3]) #Should be <oligomersize1+oligomersize2-4 and >oligomersize1

#Number of frames in new trajectory
num_frames = 15001

#Time-step between frames in new trajectory
time_step = 4000


##Some information use for choosing subunits and domains for superposing and slicing

#Number of atoms per SPOP subunit
atoms_per_subunit = 5241

#Start of BTB/BACK in monomer
first_atom_BTBBACK = 2580


##Settings for checking clashes and/or contact number between newly joined subunits

#Set clashes_test = True to test if any atoms are closer than a distance cut-off when joining new subunits
clashes_test = True

#Interatomic distance cut-off for steric clash - if atoms are closer than this distance, it counts as a clash
atomic_clash_distance=0.04  #0.4



#Set contacts_test = True to test if there are more contacts when joining new subunits than there normally is 
contacts_test = True

#Atom-atom distance cutoff for atomic contact
contact_distance_cutoff = 5

#Number of frames (evenly spaced) in input trajectory to use for calculating the standard number of contacts 
#Set low enough that frames are uncorrelated, so normal standard deviation can be used instead of block-error analysis
nr_frames_contacts_ref = 30

print("Combining %i-mer and %i-mer trajectories into %i-mer trajectory" % (oligomersize1, oligomersize2, target_length))

#Load trajectories of SPOP1 and SPOP2
print("Loading trajectories")
SPOP_traj_1 = md.load('SPOP_%imer/prodrun_AAbackmapped.xtc' % oligomersize1, top='SPOP_%imer/prodrun_AAbackmapped.gro' % oligomersize1)
SPOP_traj_2 = md.load('SPOP_%imer/prodrun_AAbackmapped.xtc' % oligomersize2, top='SPOP_%imer/prodrun_AAbackmapped.gro' % oligomersize2)

#Remove last two subunits of oligomer1, so you don't get "outer" subunits inside the output oligomer
oligomersize1 -= 2
target_length_slice_SPOP1 = np.arange(atoms_per_subunit*oligomersize1)
SPOP_traj_1 = md.Trajectory.atom_slice(SPOP_traj_1, target_length_slice_SPOP1)

#Determine max output oligomer size and do sanity check of target length
new_oligomersize = oligomersize1 + oligomersize2 - 2
assert target_length<=new_oligomersize and target_length>oligomersize1, "Target length is not possible with given oligomer sizes."

#Take an atom slice out of SPOP2 trajectory to reach target oligomer length
#Keep final "outer" subunit, so slice the oligomer from the beginning
if target_length < new_oligomersize:
    #Define new oligomer2 size
    size_diff = new_oligomersize - target_length
    oligomersize2_new = oligomersize2 - size_diff
    
    #Slice out the first subunits of oligomer2, to get the right length
    target_length_slice_SPOP2 = np.arange(atoms_per_subunit*size_diff,atoms_per_subunit*oligomersize2)
    SPOP_traj_2 = md.Trajectory.atom_slice(SPOP_traj_2, target_length_slice_SPOP2)
    
    #Redefine oligomersize2 and new_oligomersize to fit with the slicing and target length
    oligomersize2 = oligomersize2_new
    new_oligomersize = oligomersize1 + oligomersize2 - 2

#Number of BTB/BACK atoms
atoms_per_BTBBACK = atoms_per_subunit-first_atom_BTBBACK

#Define atoms for superposing in SPOP1:
first_BTB_BACK_start_1 = atoms_per_subunit*(oligomersize1-1)-atoms_per_BTBBACK
first_BTB_BACK_end_1 = atoms_per_subunit*(oligomersize1-1)

second_BTB_BACK_start_1 = atoms_per_subunit*oligomersize1-atoms_per_BTBBACK
second_BTB_BACK_end_1 = atoms_per_subunit*oligomersize1

BTB_BACK_1 = np.concatenate([np.arange(first_BTB_BACK_start_1, first_BTB_BACK_end_1), np.arange(second_BTB_BACK_start_1, second_BTB_BACK_end_1)], axis=None)

#Define atoms for superposing in SPOP2
first_BTB_BACK_start_2 = atoms_per_subunit-atoms_per_BTBBACK
first_BTB_BACK_end_2 = atoms_per_subunit

second_BTB_BACK_start_2 = atoms_per_subunit*2-atoms_per_BTBBACK
second_BTB_BACK_end_2 = atoms_per_subunit*2

BTB_BACK_2 = np.concatenate([np.arange(first_BTB_BACK_start_2, first_BTB_BACK_end_2), np.arange(second_BTB_BACK_start_2, second_BTB_BACK_end_2)], axis=None)

#Define atoms for slicing out dimer in SPOP2
slice_start_2 = atoms_per_subunit*2
slice_end_2 = atoms_per_subunit*oligomersize2

SPOP_2_atoms_slice = np.arange(slice_start_2, slice_end_2)

#Find average number of contacts and std dev between subunit 2 and 3 in SPOP1 (or SPOP2 if SPOP1 is not big enough)
#This is used to set the cutoff for number of contacts allowed between newly-joined subunits when stitching

if contacts_test == True:
    
    if oligomersize1 >= 4:
        contacts_test_traj = SPOP_traj_1
        print("Finding tolerable number of contacts between joined subunits based on first BACK/BACK interface in SPOP oligomer 1")
    elif oligomersize2 >= 4:
        contacts_test_traj = SPOP_traj_2  
        print("Finding tolerable number of contacts between joined subunits based on first BACK/BACK interface in SPOP oligomer 2")
    else:
        contacts_test = False
        print("Checking contact number is not possible, as there is no BACK/BACK interface for reference in either input oligomer. Only possible with >=4-mer")

if contacts_test == True:
    
    #Make 2D array with all combinations of atoms in second and third subunit for calculating distances
    SPOP_second_subunit = np.arange(atoms_per_subunit, atoms_per_subunit*2)
    SPOP_third_subunit = np.arange(atoms_per_subunit*2, atoms_per_subunit*3)

    atoms_list_1 = np.concatenate([SPOP_second_subunit]*len(SPOP_second_subunit),axis=None)
    atoms_list_2 = np.repeat(SPOP_third_subunit, len(SPOP_second_subunit), axis=None)

    atom_pairs = np.append(atoms_list_1.reshape(-1,1), atoms_list_2.reshape(-1,1), axis=1)
    
    #For each frame, choose all distances shorter than contact cutoff as contacts
    #Count number of contacts in the frame and append to list
    nr_contacts_vs_frame = []
    stepsize_contacts_ref = int(len(contacts_test_traj)/nr_frames_contacts_ref)
    for i in np.arange(0,len(contacts_test_traj),step=stepsize_contacts_ref):
        distances = md.compute_distances(contacts_test_traj[i], atom_pairs=atom_pairs)
        contacts_distances = distances[distances<contact_distance_cutoff]
        nr_contacts = len(contacts_distances)
        nr_contacts_vs_frame.append(nr_contacts)
        
    #Find average contact number between subunits and stddev
    nr_contacts_vs_frame = np.array(nr_contacts_vs_frame)
    contacts_avg = np.average(nr_contacts_vs_frame)
    contacts_stddev = np.std(nr_contacts_vs_frame)
    
    #Set tolerable number of contacts when stitching
    contacts_clash_cutoff_upper = int(contacts_avg + 2*contacts_stddev)
    contacts_clash_cutoff_lower = int(contacts_avg - 2*contacts_stddev)
    print("Tolerable number of contacts between joined subunits is between " + str(contacts_clash_cutoff_lower) + " and " + str(contacts_clash_cutoff_upper))

#Make 2D array with all combinations of atoms in SPOP1 last subunit and SPOP2 first subunit for calculating distances in newly joined subunits
if clashes_test == True or contacts_test == True:

    SPOP1_last_subunit = np.arange((oligomersize1-1)*atoms_per_subunit, oligomersize1*atoms_per_subunit)
    SPOP2_first_subunit = np.arange(oligomersize1*atoms_per_subunit, (oligomersize1+1)*atoms_per_subunit)

    atoms_list_1 = np.concatenate([SPOP1_last_subunit]*len(SPOP1_last_subunit),axis=None)
    atoms_list_2 = np.repeat(SPOP2_first_subunit, len(SPOP1_last_subunit), axis=None)

    atom_pairs = np.append(atoms_list_1.reshape(-1,1), atoms_list_2.reshape(-1,1), axis=1)

#Loop over trajectory frames
for i in range(num_frames):
    print("Aligning frame " + str(i))
    
    clash = True

    while clash == True:
    
        #Get random frame from SPOP traj1
        random_frame_1 = np.random.randint(0,len(SPOP_traj_1))
        SPOP_frame_1 = SPOP_traj_1[random_frame_1]
    
        #Get random frame from SPOP traj2
        random_frame_2 = np.random.randint(0,len(SPOP_traj_2))
        SPOP_frame_2 = SPOP_traj_2[random_frame_2]
    
        #Superpose last dimer in SPOP1 to first dimer in SPOP2
        SPOP_frame_1 = md.Trajectory.superpose(SPOP_frame_1, SPOP_traj_2, frame=random_frame_2, atom_indices=BTB_BACK_1, ref_atom_indices=BTB_BACK_2)
    
        #Slice out first dimer of SPOP2
        SPOP_frame_2 = md.Trajectory.atom_slice(SPOP_frame_2, SPOP_2_atoms_slice)
    
        #Combine SPOP1 and SPOP2 into one frame
        combined = md.Trajectory.stack(SPOP_frame_1, SPOP_frame_2)
    
        if clashes_test == True or contacts_test == True:
    
            #Calculate interatomic distances
            distances = md.compute_distances(combined, atom_pairs=atom_pairs)  

            #Get shortest distance
            shortest_distance = np.amin(distances)

            #Calculate number of contacts (atom pairs with distance < contact_distance_cutoff)
            contacts_distances = distances[distances<contact_distance_cutoff]
            nr_contacts = len(contacts_distances)

        #Check for clashes and number of contacts
        if contacts_test == True and clashes_test == True and shortest_distance <= atomic_clash_distance and (nr_contacts < contacts_clash_cutoff_lower or nr_contacts > contacts_clash_cutoff_upper):
            print("There was a steric clash and there was an intolerable number of contacts! Selecting new structures.")
            print("Shortest interatomic distance was: " + str(shortest_distance*10) + 'Å')
            print("Number of contacts was " + str(nr_contacts))
        elif clashes_test == True and shortest_distance <= atomic_clash_distance:
            print("There was a steric clash! Selecting new structures.")
            print("Shortest interatomic distance was: " + str(shortest_distance*10) + 'Å')
        elif contacts_test == True and (nr_contacts < contacts_clash_cutoff_lower or nr_contacts > contacts_clash_cutoff_upper):
            print("There was an intolerable number of contacts! Selecting new structures.")
            print("Number of contacts was " + str(nr_contacts))
        else:
            clash = False
    
    #Change timestamp
    combined.time = i*time_step
    
    #Add to trajectory (or make new trajectory if first frame)
    if i==0:
        traj_combined = combined
    else:
        traj_combined = md.Trajectory.join(traj_combined, combined)

#Center stitched oligomer in joined trajectory
print("Centering new trajectory")
traj_combined = md.Trajectory.center_coordinates(traj_combined)

#Print some info on the new traj
print("New trajectory is: " + str(traj_combined))

print("Writing xtc and gro files")

#Make new dir for writing xtc and gro files
os.system('mkdir SPOP_%imer' % new_oligomersize)

#Write xtc and temporary gro files
outfile_xtc = 'SPOP_%imer/prodrun_AAbackmapped.xtc' % new_oligomersize
outfile_gro = 'SPOP_%imer/prodrun_AAbackmapped.gro' % new_oligomersize
md.Trajectory.save_xtc(traj_combined, outfile_xtc)
md.Trajectory.save_gro(traj_combined[0], outfile_gro)

print("Finished!")
