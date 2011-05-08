
DELIMITER = "$$$$"

def main():
    
    import sys
    if len(sys.argv) < 2:
        print "usage: python utils.py (input) (number to extract)"
        print "or"
        print "usage: python utils.py (input_file) //for downloading images"
        sys.exit(1)
    
    filename = sys.argv[1]
    
    if len(sys.argv) == 2:
        download_images(filename)
        sys.exit(0)
    
    #below, TRUNCATE is being run
    try:
        n = int(sys.argv[2])
    except ValueError:
        print "second parameter must be an integer"
        sys.exit(1)
        
    if n < 1:
        print "number to extract must be > 0"
        sys.exit(1)
    
    truncate_sdf(filename, n)

def download_images(filename):

    '''
    for a given results file, downloads all of the images of the associated compounds
    '''

    import subprocess
    folder_name = filename.split(".")[0] + '-img'
    target_image = folder_name + "/%i.png"
    url = "http://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?t=l&cid=%i"
    
    subprocess.call(["rm","-rf", folder_name])
    subprocess.call(["mkdir",folder_name])
    
    f = open(filename, "r")
    
    for line in f:
        l = line.strip().split("\t")
        if len(l) == 4:
            sym = int(l[2])
            cid = int(l[0])
            
            subprocess.call(["wget", "-O", target_image % cid, url % cid])
    
    f.close()

def truncate_sdf(filename, n):
    '''
    given an SDF file (filename) and an integer (n), this spits out n molecules' entries
    to standard out
    
    useful for generating small test-sets
    '''
    f = open(filename, "r")
    
    for line in f:
        line = line.strip()
        print line
        
        if line.startswith(DELIMITER):
            
            n -= 1
            if n == 0:
                f.close()
                return
    
    f.close()
    
if __name__ == "__main__":
    main()