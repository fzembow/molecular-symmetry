
DELIMITER = "$$$$"

def main():
    
    import sys
    if len(sys.argv) != 3:
        print "usage: python extract_sdf.py (input) (number to extract)"
        sys.exit(1)
    
    filename = sys.argv[1]
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