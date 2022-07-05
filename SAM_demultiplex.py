# Used to demultiplex an Illumina Hi-seq run using either one or two indexing primers
# Written for Python 3.5+
# Author: Michael A. Schon

import os
import sys
import argparse

####################
### USER OPTIONS ###
####################

def get_arguments():
    '''
    Processes input arguments and returns an argparse object.
    Allowed options:
        use_indices     - 'i5','i7','both'
        experiment_type - 'SE','PE'
    For paired-end (PE) experiments, two FASTQ files will be output for each sample, containing mate pairs 1 and 2.
    '''
    desc = (
        "Takes a SAM file and a table of index sequences,"
        "and demultiplexes the SAM into FASTQ file(s) for"
        "each sample in the table."
    )
    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # add arguments to the ArgumentParser
    parser.add_argument(
        '-S','--sam', dest='sam', type=str, 
        help='filepath to multiplexed SAM file',
        required=True
    )
    parser.add_argument(
        '-T', '--table', dest='table', type=str, 
        help='reference table of indices',
        required=True
    )
    parser.add_argument(
        '-I','--indices', dest='indices', type=str, 
        help='which index feature(s) to demultiplex on',
        choices=['i5','i7','both'], default='both'
    )
    parser.add_argument(
        '-M','--mismatch', dest='mismatch', type=int, 
        help='number of mismatches to allow in index sequence',
        default=0, choices=[0,1,2,3]
    )
    parser.add_argument(
        '-O', '--output', dest='output', type=str, 
        help='output directory',
        default='.'
    )
    parser.add_argument(
        '--paired', dest='paired', action='store_true',
        help='is the sample paired-end?'
    )
    args = parser.parse_args()
    return args


#################
### FUNCTIONS ###
#################

def get_indices(index_table, use_i5, use_i7):
    '''
    Generates a dictionary of index sequences
    given a filepath to an index table and a
    specification of which indices to use (i5, i7, or both).
    '''
    file = open(index_table,'r')
    indices = {}
    for line in file:
        if len(line) == 0 or line[0] == '#': continue
        try:
            label,ind_i7,ind_i5 = line.rstrip().split('\t')
            if len(ind_i5) > 8:ind_i5 = ind_i5[:8]
            if len(ind_i7) > 8:ind_i7 = ind_i7[:8]
            if use_i5 and use_i7:
                indices[ind_i5+'_'+ind_i7] = label
            elif use_i5:
                indices[ind_i5] = label
            elif use_i7:
                indices[ind_i7] = label
            else:
                print(use_i5, use_i7)
                print("ERROR: cannot use neither index to demultiplex")
                sys.exit(1)
        except:
            label,ind = line.rstrip().split('\t')
            indices[ind]=label
    
    return indices


def open_output_files(output_folder, index_dict, paired):
    if not output_folder.endswith('/'):
        output_folder += '/'
    
    output_files = {}
    output_files['unassigned'] = open(output_folder+'unassigned.fastq','w')
    for k,v in index_dict.items():
        if paired:
            output_files[k+'/1'] = open(output_folder+v+'.1.fastq','w')
            output_files[k+'/2'] = open(output_folder+v+'.2.fastq','w')
        else:
            output_files[k] = open(output_folder+v+'.fastq','w')
    
    return output_files


def below_hamming_distance(a, b, distance):
    '''
    Compares strings a and b and returns a bool:
    True if Hamming distance is <= specified distance,
    else False.
    '''
    assert len(a) == len(b)
    if a == b: return True
    d = 0
    for i in range(len(a)):
        d += int(a[i] != b[i])
        if d > distance:
            return False
    
    return True


def best_index_match(index, list_of_indices, distance):
    '''
    Given an index sequence list of sequences to compare to,
    returns an index from the list if and only if:
    the Hamming distance to this index is <= distance,
    and this is true of no other indices in the list.
    '''
    match = 'unassigned'
    if distance == 0:
        return match
    
    for listed_index in list_of_indices:
        if below_hamming_distance(index, listed_index, distance):
            if match == 'unassigned':
                match = listed_index
            else:
                return 'unassigned'
    
    return match
    

#########################
### EVALUATE SAM FILE ###
#########################

def main():
    args = get_arguments()
    use_i5 = args.indices in ['i5','both']
    use_i7 = args.indices in ['i7','both']
    indices = get_indices(args.table, use_i5, use_i7)
    index_sequences = sorted(indices.keys())
    if not os.path.exists(args.output):
        os.mkdir(args.output)
    
    output_files = open_output_files(args.output, indices, args.paired)
    sam_file = open(args.sam)
    flags = {'4':'', '77':'/1', '141':'/2'}
    for line in sam_file:
        if line[0] == '@':continue
        qname,flag,rname,pos,mapq,cigar,mrnm,mpos,isize,seq,qual,ind_i5,ind_i5q,ind_i7,rg,qt = line.rstrip().split('\t')
        ind_i5 = ind_i5.replace('B2:Z:','')[:8]
        ind_i7 = ind_i7.replace('BC:Z:','')[:8]
        line_index = (['',ind_i5][use_i5]+'_'+['',ind_i7][use_i7]).strip('_')
        if line_index in indices.keys():
            output_files[line_index+flags[flag]].write('\n'.join(['@'+qname+flags[flag],seq,'+',qual])+'\n')
        else:
            match = best_index_match(line_index, index_sequences, args.mismatch)
            if match == 'unassigned':
                output_files['unassigned'].write('\n'.join(['@'+qname+flags[flag],seq,'+',qual])+'\n')
            else:
                output_files[match+flags[flag]].write('\n'.join(['@'+qname+flags[flag],seq,'+',qual])+'\n')

    for i in indices.keys():
        if args.paired:
            output_files[i+'/1'].close()
            output_files[i+'/2'].close()
        else:
            output_files[i].close()


if __name__ == '__main__':
    main()
