

#Project for Python for Genomic Data Science

#First step - Create a dictionary of the fasta file and count the number of records.

def records(filename):
    #open file
    try:
        f = open(filename)
    except IOError:
        print("File %s does not exist!!" %filename)
#create a dictionary of the ids and sequences 
    seqs = {}
    for line in f:
    #let's discard the newline at the end(if any)
        line = line.rstrip()
    # distinguish header from sequence
        if line[0] == '>':
            words = line.split()
            name = words[0][1:]
            seqs[name] = ''
        else: # sequence, not header
            seqs[name] = seqs[name] + line
    f.close()
    #no. of records in the dictionary or file
    print("No. of records is " + str(len(seqs)))
         

# Second step - To sort and count the length of the sequences(values). Get the longest and shortest ones.

def longshortseq(filename):
     #open file
    try:
        f = open(filename)
    except IOError:
        print("File %s does not exist!!" %filename)
#create a dictionary of the ids and sequences 
    seqs = {}
    for line in f:
    #let's discard the newline at the end(if any)
        line = line.rstrip()
    # distinguish header from sequence
        if line[0] == '>':
            words = line.split()
            name = words[0][1:]
            seqs[name] = ''
        else: # sequence, not header
            seqs[name] = seqs[name] + line
    f.close()
    # count the longest and shortest sequences in the dictionary.
    import math
    longest = 0
    shortest = math.inf
    longseqs = {}
    shortseqs = {}
    for name, seq in seqs.items():
        if len(seq) > longest:
            longest = len(seq)
        if len(seq) < shortest:
            shortest = len(seq)
    for name, seq in seqs.items():
        if len(seq) == longest:
            longseqs[name] = seq
        if len(seq) == shortest:
            shortseqs[name] = seq
    for name, seq in longseqs.items():
        print(name, seq)
        print("The lenth of the longest sequence is: %s" %len(seq))
    for name, seq in shortseqs.items():
        print(name, seq)
        print("The length of the shortest sequence is: %s" %len(seq))


# Third step - i) To identify all ORFs in each sequence of a fasta file
# ii) what is the length of the longest ORF in the file? 
# iii) What is the identifier of the sequence containing the longest ORF? 
# iv) For a given sequence identifier, what is the longest ORF contained in the sequence represented by that identifier? 
# v) What is the starting position of the longest ORF in the sequence that contains it? The position should indicate the character number in the sequence.
    
def find_orf(frame, filename):
    if frame == 1:
        fr = 0
    elif frame == 2:
        fr = 1
    elif frame == 3:
        fr = 2
    else:
        print("Type a valid frame number.")   
    orf = [] 
    
    from Bio import SeqIO
    for record in SeqIO.parse(open(filename), "fasta"):
        name, sequence = record.id, str(record.seq)
        start_indexes = []
        stop_indexes = [] 
        stops =["TAA", "TAG", "TGA"]          
    #Find all ATG indexes
        for i in range(fr, len(sequence), 3):
            if sequence[i:i+3] == "ATG":
                start_indexes.append(i)
              
    # Find all stop codon indexes
        for i in range(fr, len(sequence), 3):
            if sequence[i:i+3] in stops:
                stop_indexes.append(i)
               
             
        mark = 0
        for i in range(0,len(start_indexes)):
            for j in range(0, len(stop_indexes)):
                if start_indexes[i] < stop_indexes[j] and start_indexes[i] > mark:
                    orf.append(["Id:"+name,"Start Position:"+str(start_indexes[i]),"Length of Sequence:"+str(len(sequence[start_indexes[i]:stop_indexes[j]+3]))])
                    mark = stop_indexes[j]+3
        
    return orf

def findlongorf(frame, filename):
    
    if frame == 1:
        fr = 0
    elif frame == 2:
        fr = 1
    elif frame == 3:
        fr = 2
    else:
        print("Type a valid frame number.")   
    orf = [] 
    
    from Bio import SeqIO
    for record in SeqIO.parse(open(filename), "fasta"):
        name, sequence = record.id, str(record.seq)
        start_indexes = []
        stop_indexes = [] 
        stops =["TAA", "TAG", "TGA"]          
    #Find all ATG indexes
        for i in range(fr, len(sequence), 3):
            if sequence[i:i+3] == "ATG":
                start_indexes.append(i)
              
    # Find all stop codon indexes
        for i in range(fr, len(sequence), 3):
            if sequence[i:i+3] in stops:
                stop_indexes.append(i)
        longorf = 0       
        lenorf = 0     
        mark = 0
        for i in range(0,len(start_indexes)):
            for j in range(0, len(stop_indexes)):
                if start_indexes[i] < stop_indexes[j] and start_indexes[i] > mark:
                    lenorf = len(sequence[start_indexes[i]:stop_indexes[j]+3])
                    orf.append(["Id:"+name,"Start Position:"+str(start_indexes[i]),"Length of Sequence:"+str(len(sequence[start_indexes[i]:stop_indexes[j]+3]))])
                    mark = stop_indexes[j]+3
                    if lenorf > longorf:
                        longorf = lenorf
                        longid = name
                        longstart = start_indexes[i]
                    else:
                        pass
                        
        print("Id: " +longid + " Longest ORF is " +str(longorf) + " starting at position "+str(longstart))
        
    return orf

def find_orfseq(seqId):
     
    orf = [] 
     
    from Bio import SeqIO
    for record in SeqIO.parse(open("dna2.fasta"), "fasta"):
        name, sequence = record.id, str(record.seq)
        if name == seqId:
            for f in range(0,3):
                start_indexes = []
                stop_indexes = [] 
                stops =["TAA", "TAG", "TGA"]          
    #Find all  ATG indexes
                for i in range(f, len(sequence), 3):
                    if sequence[i:i+3] == "ATG":
                        start_indexes.append(i)
              
    # Find all stop codon indexes
                for i in range(f, len(sequence), 3):
                    if sequence[i:i+3] in stops:
                        stop_indexes.append(i)
                longorf = 0
                lenorf = 0
                mark = 0
                for i in range(0,len(start_indexes)):
                    for j in range(0, len(stop_indexes)):
                        if start_indexes[i] < stop_indexes[j] and start_indexes[i] > mark:
                            lenorf = len(sequence[start_indexes[i]:stop_indexes[j]+3])
                            orf.append(["RF:"+str(f+1),"Start Position:"+str(start_indexes[i]),"Length of Sequence:"+str(len(sequence[start_indexes[i]:stop_indexes[j]+3]))])
                            mark = stop_indexes[j]+3
                            if lenorf >= longorf:
                                longorf = lenorf
                                longstart = start_indexes[i]
                                rf = f+1        
                print("Reading frame "+str(rf)+" Longest ORF is " +str(longorf) + " starting at position "+str(longstart))
        
            return orf

# Fourth Step - i) Given a length n, identify all repeats of length n in all sequences in the FASTA file. 
# ii) How many times each repeat occurs in the file?
# iii) Which is the most frequent repeat of a given length?


def Arunfind_repeats(n, filename):
   
    
    sequence = ""
    position_tracker = {}
    
    from Bio import SeqIO 
    with open(filename, "rU") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequence += str(record.seq)
    iseq={}              
    for i in range(len(sequence)-n+1):
        repeat = sequence[i:i+n]
        if (repeat not in iseq):
            iseq[repeat]=1
            for j in range(len(sequence)+1):
                 # The comparison sequence to compare to the target.
                compseq=""
                if j!=i:
                    compseq=sequence[j:j+n]
                    
                if compseq == repeat:
                    if (compseq in position_tracker) == False:
                        position_tracker[compseq] = [i,j] 
                      
                    else:  
                       k=position_tracker[compseq]
                       if i not in k:
                           (position_tracker[compseq]).append(i)
                       elif j not in k:
                            (position_tracker[compseq]).append(j)
                             
    for k in position_tracker:
        sorted(position_tracker, key=lambda k: len(position_tracker[k]), reverse=True)
        if (len(position_tracker[k])) > 3:
            print (k, len(position_tracker[k]))
        
    
   

def find_repeats(n, filename):
    import operator
    exp = []
    sequence = ""
    repeats = {}
    from Bio import SeqIO
    with open(filename, "rU") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequence += str(record.seq)
    
    for i in range(len(sequence)):
        repeat = sequence[i:i+n]
        if repeat not in exp:
            for j in range(len(sequence)):
                compseq = sequence[j:j+n] # The comparison sequence to compare to the target.
                if j > i and repeat == compseq:
                    if repeat in repeats:
                        repeats[repeat] = repeats[repeat]+1
                    elif repeat not in repeats:
                        repeats[repeat] = 2
            exp.append(repeat)
    sorted_repeats = sorted(repeats.items(), key=operator.itemgetter(1))
    print(sorted_repeats)

