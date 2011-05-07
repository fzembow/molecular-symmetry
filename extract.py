#!/usr/bin/env python

##
## for testing purposes: extract N records from the beginning of an SDF file that contains many
## and prints it to stdout
##

import re, math, os
from numpy import array, sort, empty, empty_like, float32, int32, around, matrix, asarray

#whether to use GPU-accelerated atomic distance calculations
USE_GPU = False
if USE_GPU:
    import pycuda.driver as cuda
    import pycuda.autoinit
    from pycuda import compiler, gpuarray
    
    block_dim=(8,8,1)
    THREAD_X = 8.0
    THREAD_Y = 8.0
    
    mod = compiler.SourceModule("""
      #include <math.h>
    
      __global__ void distances(float *a, float *result, int N)
      {
        int idx = threadIdx.x + __umul24(blockIdx.x,8);
        int idy = threadIdx.y + __umul24(blockIdx.y,8);            
        
        if(idx < N && idy < N)
        {
        result[idx + __umul24(idy,N)] = floorf(sqrt(pow(a[__umul24(idx,3)] - a[__umul24(idy,3)],2) + pow(a[__umul24(idx,3) + 1] - a[__umul24(idy,3) + 1],2) + pow(a[__umul24(idx,3) + 2] - a[__umul24(idy,3) + 2],2)) * 1000 +  0.5) / 1000;
        }
      }
      """)
      
    #get the function on the GPU
    func = mod.get_function("distances")

#whether to use multiple processes, sharing molecules with a Queue
#NOTE: you can't do multiprocessing with gpu acceleration unless you have 
#multiple GPUs
USE_MULTIPROCESSING = True

if USE_MULTIPROCESSING:
    from multiprocessing import Queue, Process, cpu_count
    NUM_PROCESSES = cpu_count()
    WORK_SIZE = 2
    QUEUE_TIMEOUT = 0.25

#whether to show timing information
TIMING = False

#whether to exclude single atoms from the result (which have infinite sym.)
EXCLUDE_SINGLE_ATOM = True

#how many decimal places to keep when computing distance
PRECISION = 3

#how much error is tolerated in computing sameness
TOLERANCE = 0.075

#for processing SDF files
#the specification can be found at http://www.symyx.com/downloads/public/ctfile/ctfile.pdf
DELIMITER_RE = re.compile("$$$$")
#MOLFILE_DELIMITER = "V2000"
ATOM_HEADER_RE = re.compile("^(\d+)\s+(?:\d+\s+){7,8}V2000")
#ATOM_RE = "^\t?(\D?\d+\.\d+) +(\D?\d+\.\d+) +(\D?\d+\.\d+) (\w+)"
ID_FIELD_RE = re.compile("^> <PUBCHEM_COMPOUND_CID>")

def main():
    
    import sys
    if len(sys.argv) != 2 and len(sys.argv) != 3:
        print "usage: python extract_sdf.py (input)"
        sys.exit(1)
    
    filename = sys.argv[1]
    if len(sys.argv) == 3:
        filename2 = sys.argv[2]
    else:
        filename2 = None

    if TIMING:
        #start keeping track of time
        import time
        start = time.time()
    
    #run calculation on all molecules with multiple processes
    if USE_MULTIPROCESSING:
        process_multiple(filename, filename2)
        
    #run calculation on all molecules    
    else:
        process_single(filename)
    
    if TIMING:    
        print (time.time() - start), "seconds elapsed (total)"
        
def process_single(filename):
    '''
    single process detection of symmetry
    '''
    for molfile in extract_molfiles(filename):
        sym = calculate_fold_symmetry(molfile['atoms'])
        if sym > 1:
            type = calculate_rotation(molfile)
            print molfile['id'], sym, type
            
def process_multiple(filename, secondfile=None):
    '''
    multiprocess detection of symmetry
    '''
    if TIMING:
        import time
        start = time.time()
            
    #if we are syncing with a 3D file
    if secondfile is not None:
        two_molfiles = True
        second_generator = extract_molfiles(secondfile)
        molfile2 = second_generator.next()
    else:
        two_molfiles = False
    
    #add all molecules in a file to the queue, in lists of size WORK_SIZE
    q = Queue()
    work = []
    for molfile in extract_molfiles(filename):
        
        #if we have a molfile from the second file, then append that to the work unit
        if two_molfiles:
            if molfile2['id'] == molfile['id']:
                molfile['atoms3d'] = molfile2['atoms']
                try:
                    molfile2 = second_generator.next()
                except StopIteration:
                    two_molfiles = False
      
        work.append(molfile)
        if len(work) < WORK_SIZE:
            continue
        q.put(work)
        work = []
    
    #for imperfect rounding
    if len(work) > 0:
        q.put(work)
    
    if TIMING:  
        print (time.time() - start), "seconds loading molecule file"
    
    #spawn out the processes
    processes = []
    for i in range(NUM_PROCESSES):
        p = Process(target=worker, args=(q,))
        processes.append(p)
        p.start()

    #wait for all of the processes to finish before returning
    for p in processes:
        p.join()


def worker(q, ):
    '''
    multiprocess worker, gets its molfiles from queue q, prints to stdout
    '''
    
    #keep fetching items from the queue to multiprocess
    from Queue import Empty
    while True:
        try:
            work = q.get(True, QUEUE_TIMEOUT)    
            #try all of the molfiles in this work unit, stop when we encounter symmetry.
            # work_unit[0] is the 2D molfile
                
            for molfile in work:
                sym = calculate_fold_symmetry(molfile['atoms'])
                if sym > 1:
                    print molfile['id'], sym
                    break
                    
        #if there are no more items on the Queue then we are done
        except Empty:
            return

def extract_molfiles(filename):
    '''
    given an SDF file (filename), this yields dictionaries of {'id', 'atoms'}
    
    atoms is a list of 4-tuples (x,y,z,element)
    '''
    f = open(filename, "r")
    
    molecule = {}
    seen_atom = False
    
    #iterate through the compound and dump out molfile information one at a time
    while True:
        try:
            line = f.next().lstrip()
        except StopIteration:
            break
        
        #parse all of a molecule's atoms once we encounter one, to save on regexp time
        if not seen_atom:
            
            header = ATOM_HEADER_RE.match(line)
            if header:
                
                if len(header.group(1)) == 5:
                    num_atoms = int(line[0:2])
                else:
                    num_atoms = int(line[0:3])#header.group(1))
                atom_list = []
                
                for i in range(num_atoms):
                                        
                    line = f.next().lstrip()
                    line = re.split("\s+",line,4)
                    
                    #for j in range(3):
                    #    np_atom_list[i][j] = float(line[0])
                    
                    atoms = (float(line[0]), float(line[1]), float(line[2]))#, line[3])
                    #atoms = (float(line[0]), float(line[1]), float(line[2]), line[3])
                    atom_list.append(atoms)
                
                #molecule['atoms'] = np_atom_list
                molecule['atoms'] = atom_list
                seen_atom = True

        elif ID_FIELD_RE.match(line):
            pubchem_id = int(f.next().lstrip())
            molecule["id"] = pubchem_id
        elif DELIMITER_RE.match(line):
            
            yield molecule
            molecule = {}
            seen_atom = False
    
    f.close()
    
def arrays_equal(a,b):
    '''
    determines whether two arrays are equal, elementwise. 
    NOTE: assumes the arrays are the same length, and that the first element is always the same
    '''

	#Method that's currently working.
	#for i in range(1,len(a)):
    #    if math.fabs(a[i] - b[i]) > TOLERANCE:
    #        return False
    
    false_count = 0.0

    for i in range(1, len(a)):
        if math.fabs(a[i] - b[i]) > TOLERANCE:
            false_count += 1.0

    full = len(a) - 1.0

    if  (false_count / full) < 0.1:
        return True
    else:
        return False
    
def check_max_symmetry_support(a, N):
    '''
    finds the largest-fold symmetry that a given atom is capable of supporting, when in the center of a molecule
    [ 0.     0.62   0.62   0.62   0.62 ] -> 4
    '''
    
    index = 1
    best_count = 99999
    while index < N-1:
        count = 1
        while index < N-1 and math.fabs(a[index] - a[index+1]) < TOLERANCE: 
            count += 1
            index += 1
        
        best_count = min(count, best_count)
        
        index += 1
        
    return best_count
    
def min(a, b):
    '''
    returns the lesser of two values
    '''
    if a < b:
        return a
    return b

##TODO: incorporate the element as another element of distance
def dist(coord1, coord2):
    '''
    calculates the distance between two (x,y,z,...) tuples (coord1, coord2)
    '''
    
    #TODO: don't need to do the square root -- but need to adjust TOLERANCE to compensate. however, this only saves like 2 seconds on the entire calculation... forget it
    
    return round(math.sqrt( (coord1[0] - coord2[0]) ** 2 + (coord1[1] - coord2[1]) ** 2 + (coord1[2] - coord2[2]) ** 2 ),PRECISION)
    
def calculate_rotation(molecule):
    atoms = molecule['atoms']

    centroid = find_centroid(atoms)

    new_coordinates = []

    for atom in atoms:
        mat = matrix([[atom[0] - centroid[0]], [atom[1] - centroid[1]], [atom[2] - centroid[2]]])
        # find new coordinates by checking z_
        new_coord = find_new_coord(mat)
        new_coordinates.append((new_coord[0,0] + centroid[0], new_coord[1,0]+centroid[1], new_coord[2,0] + centroid[2]))
        
    atoms_s =  sort(atoms, axis=0)
    new_coordinates_s = sort(new_coordinates, axis=0)
    print atoms
    print atoms_s
    #if arrays_equal(atoms_s, new_coordinates_s):
    #    return "C2 along z-axis"
    #else:
    #    return "No C2"

    #print centroid

def find_new_coord(mat):
   C2_z_mat = matrix('-1 0 0; 0 -1 0; 0 0 1')
   new_coord = C2_z_mat * mat
   #return asarray(new_coord)
   return asarray(new_coord)

def find_centroid(atoms):
    x_sum = 0.0
    y_sum = 0.0
    z_sum = 0.0
    count = 0
    for atom in atoms:
        count += 1
        x_sum += atom[0]
        y_sum += atom[1]
        z_sum += atom[2]

    return (x_sum/count, y_sum/count, z_sum/count)

def calculate_fold_symmetry(atoms):
    
    #number of atoms
    N = len(atoms)
        
    #single-atom molecules are excluded
    if EXCLUDE_SINGLE_ATOM and N == 1:
        return 1
    
    #distance matrix on CPU
    distances = array([[0.0] * N] * N)
    if USE_GPU:
        
        #transfer data to GPU
        import time
        start = time.time()
        
        a = array(atoms, float32)
        a_gpu = cuda.mem_alloc(a.nbytes)
        cuda.memcpy_htod(a_gpu, a)
        
        grid_dim=(int(math.ceil(N/THREAD_X)),int(math.ceil(N/THREAD_Y)) )
        
        #create results array
        distances_gpu = gpuarray.empty((N, N), float32)
        
        func(a_gpu, distances_gpu, int32(N), block=block_dim, grid=grid_dim)
        
        distances = distances_gpu.get()
        
    else:
        for i,atom in enumerate(atoms):
            for j,target in enumerate(atoms):        
                if i != j and j >= i:   #lower left corner of matrix
                    distances[i][j] = distances[j][i] = dist(atom, target) 
        
    #NOTE: at this point, each row of the distances array
    # contains the distance from that atom to each of the atoms:
    #
    # distance = distances[src][target]
    
    # sort distance matrix
    distances = sort(distances)
    distances = sort(distances, axis=0)    
    # now distances are sorted by row and also by column
    
    #TODO:do we need to do this for all atoms or can we stop at a certain point?
    
    #check for symmetry
    fold_symmetry = 99999 #the smallest global symmetry found
    index = 0
    while index < (N-1):
        current_symmetry = 1
        while index < N-1 and arrays_equal(distances[index],distances[index+1]):
            index += 1
            current_symmetry += 1
                
        #check if it's an atom in the middle of symmetry
        if current_symmetry == 1:
            current_symmetry = check_max_symmetry_support(distances[index], N)
        
        #find the lowest feasible symmetry
        fold_symmetry = min(fold_symmetry, current_symmetry)
        if fold_symmetry == 1:
            break
        
        index += 1
    
    return fold_symmetry
    
    #TODO: calculate planes of symmetry?
	# figure out vizualization
	# figure out threshold
	
    
if __name__ == "__main__":
    main()