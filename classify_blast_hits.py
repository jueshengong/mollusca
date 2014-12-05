#!/usr/bin/env python


import os
import logging
import argparse
import gzip
import array
import os.path
import subprocess
import numpy

#############################
## GLOBAL CONSTANTS #########
#############################

PARENT_NAME = 'Mollusca'
PARENT_TAXA_ID = '6447'
PARENT_LIST = '/home/juesheng.ong/planktons/taxid/mollusca/members.txt'

PARENT_CODE    = 1
NO_PARENT_CODE = 2

GI_PROT_DMP = '/home/juesheng.ong/NCBI_taxa/gi_taxid_prot/gi_taxid_prot.dmp'

MIN_BIT_SCORE = 100
MAX_E_VALUE = 1e-8
BIT_DELTA = 100
OTHERS_SPEC_LIST = '/home/juesheng.ong/NCBI_taxa/other_spec_list/'


###Simplify blast output#####


#################################################
######## UTILITY FUNCTION #######################
#################################################

def iterFasta(path):
    """Iterates over the sequences of a fasta file"""
    logging.info("Loading fasta files from %s" % path)
    name = None
    seq = []
    if path.endswith('.gz') or path.endswith('.gz"'):
        handle = gzip.open(path)
    else:
        handle = open(path)
    for line in handle:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name:
                yield (name, ''.join(seq))
            name = line[1:]
            seq = []
        else:
            seq.append(line)
    if name:
        yield (name, ''.join(seq))
    handle.close()


def giReference(gi_path):
    
    def getMaxGI(gi_path):
        handle = open(gi_path,'r')
        maxval = 0
        for line in handle:
            line = line.strip()
            if not line:
                continue
            fields = line.split('\t')
            pos = int(fields[0])
            if pos > maxval:
                maxval = pos
        handle.close()
        return maxval
    
    maxval = getMaxGI(gi_path)            
    d = maxval + 1  ### FIXME calculate yourself
    testgi = numpy.zeros(d)
    handle = open(gi_path,'r')
    logging.info('Start reading gi_prot dmp file')
    for line in handle:
        line = line.strip()
        if not line:
            continue
        fields = line.split('\t')
        tax_id = int(fields[1])
        pos = int(fields[0])
        testgi[pos] = tax_id
    handle.close()
    logging.info('GI conversion list has been created in memory')
    return testgi


def loadFamily(path):
    '''This function loads the member of a group into a list'''
    family_list = []
    path = str(path)
    print "Openning handle.."
    handle = open(path, 'r')
    for line in handle:
        line = line.strip()
        if line == '':
            continue
        tax_id = line
        family_list.append(tax_id)
    print "completely loaded member list in %s size %d" %(path, len(family_list))
    handle.close()
    return family_list


def prepareOtherTaxaDict(species_list_dir):
    '''This function loads member lists of a selection of groups into a dictionary object'''
    logging.info('Parsing classified group of members into dictionary')
    others_dict = {}
    for root, dirs, filenames in os.walk(species_list_dir):
        for f in filenames:
            species = f.split('_')[0]
            members = loadFamily(os.path.join(root,f))
            #handle = open(os.path.join(root,f),'r')
            #for line in handle:
            #    q = line.replace('\n','')
            #    members.append(q)
            others_dict[species] = members
            #handle.close()
        return others_dict

#####################################################################################
######## MAIN CLASSIFICATION PROCESS: PARSING, MERGING, SORTING, CLASSIFYING ########
#####################################################################################
                
def filterBestBit(args, querries):
    '''This function filters a blast dictionary based on 
    bitscore of each querry hits. Only hits with the highest bitscore will be kept'''
    
    def sortHits(h1, h2):
        if h1[1] > h2[1]:
            return 1
        if h1[1] < h2[1]:
            return -1
        return 0
    
    bestBitClassification = {}
    for qName in querries:
        hits = querries[qName]
        hits.sort(sortHits)
        minBitScore     = args.min_bit
        maxEvalue       = args.max_e_val
        bestBit         = 0
        bestName        = 'noHit'
        setGroup        = 1
        bestEvalue      = 1
        for hName, isGroup, bitscore, evalue in hits:
            if bitscore < minBitScore or evalue < maxEvalue:
                continue
            if bitscore > bestBit:
                bestBit     = bitscore                
                bestName    = hName
                bestEvalue  = evalue
                setGroup     = 1
        hit = (bestName, setGroup, bestBit, bestEvalue)
        bestBitClassification[qName] = hit
    logging.info('Filtered merged blast records for highest bitscore hits')
    logging.info('Found %d maximal bitscore for queries hits' % len(querries))
    return bestBitClassification

def importBlastDict(filtered_querries, querries):    
    '''This function merges the two blast query dictionaries into one'''
    logging.info("checking number of querries in both sets")
    if not len(filtered_querries) == len(querries):
        logging.debug("Size of querry sets are not the same. This might be caused by missing data. Please check the raw blast output files.")
    else:
        logging.info("Querry sets compatible. Perform merging..")
    ## p and n are used to keep track of successful merge and missing queries ##
    p = 0
    n = 0
    for qName in filtered_querries:
        slot = filtered_querries[qName]
        if qName in querries:
            querries[qName].append(slot)
            p += 1    
        else:
            querries[qName] = [slot]   #make sure query with no refseq hit, is still in list form [(hits)]
            n += 1        
    logging.info('Summary: %d successful merges and %d missing query entries' %(p,n))
    return querries
        
                
def parseBlast(blast_out):
    """Parse the blast results to be used later to prepare training and testing
    set with unambiguously classified sequences"""

    logging.info('Parsing blast results')
    if blast_out.endswith('.gz'):
        handle = gzip.open(blast_out)
    else:
        handle = open(blast_out,'r')
    querries = {}
    n        = 0
    for line in handle:
        line = line.strip()
        if not line:
            continue
        if line[0] == '#':
            continue
        fields = line.split()
        qName  = fields[0]
        hName  = fields[1]
        evalue = float(fields[10])
        isGroup = 0   ### 0: Unsure and 1: Confirmed
        bitscore = float(fields[11])
        if not qName in querries:
            querries[qName] = []
            n += 1
        hit = (hName, isGroup, bitscore, evalue)
        querries[qName].append(hit)
    logging.info('Parsed %d blast records' % n)
    logging.info('Found %d queries hits' % len(querries))
    handle.close()
    return querries


def classifyFromBlast(querries, args, gi_list, parent_id, parent_name, parent_list):
    """Classify blast results into ambiguous and unambiguous sequences from
    host and symbiont"""

    logging.info('Classifying using Blast results')
    logging.info('Parent ID is %s' %parent_id)
    c = 0
    blastClassification = {}
    trainingClassification = {}
    othersTrainClassification = {}
    othersBlastClassification = {}
    parentTrained      = 0
    noParentTrained    = 0
    parentClassified   = 0
    noParentClassified = 0
    # refError is a counter to track hits that cannot be associated with a tax-id#
    refError           = 0
    taxa_result = []
    for qName in querries:
        c += 1
        hits = querries[qName]
        #FOR DEBUG
        #print len(hits)
        #hits.sort(sortHits)
        hasParent       = False
        noParent        = False
        bestParent      = ''
        bestNonParent   = ''
        othersClassified = '' ##This is for the non-Parent hits summary later
        childClassified = ''
        #parentBestEvalue    = -1
        #noParentBestEvalue  = -1
        parentBestBit       = 0
        noParentBestBit     = 0
        minBitScore         = args.min_bit  #pls revert to 100
        maxEvalue           = args.max_e_val
        for element in hits:
            if len(element) == 1:
                print element 
                print "Error with length of hit array"
            hName    = element[0]
            isGroup  = element[1]
            bitscore = element[2]
            evalue   = element[3]
            if bitscore < minBitScore or evalue > maxEvalue:
                continue
            ###Variables are defined below because if minBitScore not met, dont have to look through GI for tax_id###
            if isGroup == 0:
                gi = hName.split('|')[1]
                gi = int(gi)
                tax_id = str(gi_list[gi])
                tax_id = tax_id.replace('.0','') # FIXME Check that only 0 after '.'
                #print tax_id
                if tax_id == '0':
                    ###to keep track of untraceable gi
                    refError += 1
                    continue
                if tax_id in parent_list:
                    ## DEBUG ##FIXME
                    if not tax_id in taxa_result:
                        taxa_result.append(tax_id) 
                    hasParent        = True
                    if bitscore > parentBestBit:
                        parentBestBit = bitscore
                        bestParent    = hName
                        childClassified = tax_id
                elif not tax_id in parent_list:
                    noParent        = True
                    if bitscore > noParentBestBit:
                        noParentBestBit = bitscore
                        bestNonParent    = hName
                        #Note: othersClassified will store taxid that will be parsed to a dictionary for non-parent hits##
                        othersClassified = tax_id
            elif isGroup == 1:
                hasParent = True
                if bitscore > parentBestBit:
                    parentBestBit = bitscore
                    bestParent    = hName
                    childClassified = parent_id
        ###PLS REMOVE THIS AFTER TESTING
        #To handle cases where best hit is the database hits, entry on othersClassified will be the parent taxa id itself
        if parentBestBit == 0 and noParentBestBit == 0:
            continue
        ###FOR THE STRICTLY UNAMBIGUOUS CASE
        if hasParent and not noParent:
            trainingClassification[qName]    = PARENT_CODE
            blastClassification[qName]       = PARENT_CODE
            othersTrainClassification[qName] = childClassified
            othersBlastClassification[qName] = childClassified
            parentTrained                   += 1
            parentClassified                += 1             
        elif noParent and not hasParent:
            trainingClassification[qName]    = NO_PARENT_CODE
            blastClassification[qName]       = NO_PARENT_CODE
            othersTrainClassification[qName] = othersClassified
            othersBlastClassification[qName] = othersClassified
            noParentTrained                 += 1
            noParentClassified              += 1
        
        if hasParent and noParent:
            ###TODO might want to add BitRatio : bitRatio = float(parentBestBit)/float(noParentBestBit)####
            bitDelta = parentBestBit - noParentBestBit
            ###FOR THE NON-STRINGENT CLASSIFICATION CASE
            if bitDelta >= 0:
                trainingClassification[qName] = PARENT_CODE
                othersTrainClassification[qName] = childClassified
                parentTrained                += 1
            else:
                trainingClassification[qName]    = NO_PARENT_CODE
                othersTrainClassification[qName] = othersClassified
                noParentTrained                 += 1
            ###FOR THE STRINGENT CLASSIFICATION CASE    
            if bitDelta > args.bit_delta:
                blastClassification[qName] = PARENT_CODE
                othersBlastClassification[qName] = childClassified
                parentClassified          += 1
            elif bitDelta <= -args.bit_delta:
                blastClassification[qName] = NO_PARENT_CODE 
                othersBlastClassification[qName] = othersClassified
                noParentClassified        += 1                
        #print "completed %d queries" %c
    print '%d variety of parent taxa members in query' %len(taxa_result)
    print taxa_result
    print "length of others classified"
    print len(othersTrainClassification)
    logging.debug("Others Blast size %d while counter shows %d " %(len(othersBlastClassification),noParentClassified))
    logging.info('Found %d non-traceable gi entries' % refError)
    logging.info('Found %d unambiguous hits' % len(trainingClassification))
    logging.info('Found %d likely %s only hits' % (parentTrained,parent_name))
    logging.info('Found %d likely non-%s only hits' % (noParentTrained,parent_name))
    logging.info('Found %d blast-classified %s hits' % (parentClassified,parent_name))
    logging.info('Found %d blast-classified non-%s hits' % (noParentClassified,parent_name))
    return trainingClassification, blastClassification, othersTrainClassification, othersBlastClassification


def seqSplit(args, trainingClassification, blastClassification):
    """Write the unambiguously classified sequences into four fasta files:
    training.fasta for host sequences, testing.fasta for host sequences,
    training.fasta for symb sequences and testing.fasta for symb sequences."""

    logging.info('Preparing training sequences')
    m = 0
    j = 0
    input_name    = args.input
    short_name    = '_'.join([input_name.split('_')[0],input_name.split('_')[1]]) 
    handle        = open(args.input,'r')
    parentTrain   = open(os.path.join(args.out_dir, short_name + '_With_' + PARENT_NAME +'_training.fasta'), 'w')
    noParentTrain = open(os.path.join(args.out_dir, short_name + '_No_'+ PARENT_NAME + '_training.fasta'), 'w')
    parentBlast   = open(os.path.join(args.out_dir, short_name + '_With_' + PARENT_NAME +'_blast_train.fasta'), 'w')
    noParentBlast = open(os.path.join(args.out_dir, short_name + '_No_' + PARENT_NAME +'_blast_train.fasta'), 'w')
    short_seq_dict = {}
    for name, seq in iterFasta(args.input):
        identity = (name.split(' ')[0])
        if len(seq) < args.seq_size:
            # Only need to initialize once!
            short_seq_dict[identity] = 1
            continue
        seqClass = trainingClassification.get(identity, 0)
        if seqClass == PARENT_CODE:
            parentTrain.write('>%s\n%s\n' % (identity, seq))
            m += 1
        elif seqClass == NO_PARENT_CODE:
            noParentTrain.write('>%s\n%s\n' % (identity, seq))
            j += 1
    for name, seq in iterFasta(args.input):
        if len(seq) < args.seq_size:
            continue
        identity = (name.split(' ')[0])
        seqClass = blastClassification.get(identity, 0)
        if seqClass == PARENT_CODE:
            parentBlast.write('>%s\n%s\n' % (identity, seq))
            m += 1
        elif seqClass == NO_PARENT_CODE:
            noParentBlast.write('>%s\n%s\n' % (identity, seq))
            j += 1
    handle.close()
    parentTrain.close()
    noParentTrain.close()
    return short_seq_dict


####################################################################
####  PREPARATION OF OTHER CLASSIFICATION SUMMARIES   ##############
####################################################################



def prepareOtherSummary(args, other_taxid_dict, othersTrainClass, othersBlastClass, short_seq_dict):
    ###total non-Parent hits: for Training  #########
    input_name = args.input
    short_name = '_'.join([input_name.split('_')[0],input_name.split('_')[1]])
    c = 0
    d = 0
    short_train = 0
    short_blast = 0
    ##Function currently not in use##
    def countTable(other_taxid_dict):
        count_dict = {}
        for group in other_taxid_dict:
            count_dict[group] = 0
        count_dict['others'] = 0
        return count_dict
        
    count_train = []
    for qName in othersTrainClass:
        if qName in short_seq_dict:
            short_train += 1
            continue
        hitId = othersTrainClass[qName]
        hitId = hitId.replace('.0','')       
        found = False
        for entry in other_taxid_dict:
            group = str(entry)
            if hitId in other_taxid_dict[entry]:
                count_train.append('%s: %s'%(group, hitId))
                found = True
                c += 1
                break
        if not found:
            count_train.append('others: %s' %hitId)   
    count_blast = []           
    for qName in othersBlastClass:
        if qName in short_seq_dict:
            short_blast += 1
            continue
        hitId = othersBlastClass[qName]
        hitId = hitId.replace('.0','') 
        #print hitId        
        found = False
        for entry in other_taxid_dict:
            group = str(entry)
            if hitId in other_taxid_dict[entry]:
                #print "FOUND OTHER HITS' FAMILY"
                count_blast.append('%s: %s'%(group, hitId))
                found = True
                d += 1
                break
        if not found:
            count_blast.append('others: %s' %hitId)
    logging.debug("ignored %d short general training sequence and %d short blast train sequence" %(short_train, short_blast))
    print "other taxa classified: Train - %d and Blast - %d" %(c,d)
    other_train_summary = open(os.path.join(args.out_dir, short_name) + '_others_train_summary.txt','w')
    other_blast_summary = open(os.path.join(args.out_dir, short_name) + '_others_blast_summary.txt','w')
    for element in count_train:
        other_train_summary.write('%s\n'%element)
    for element in count_blast:
        other_blast_summary.write('%s\n'%element)        
    other_train_summary.close()
    other_blast_summary.close()    
    logging.info('Summary for Other hits have been completed.')
    
    


def mainArgs():
    """Process command-line arguments"""

    parser = argparse.ArgumentParser(description='Perform SVM Classification of Host and Symbiont (or Parasite) Sequences')
    parser.add_argument('input',
                        help='The input fasta sequences')
    #subparsers = parser.add_subparsers(help='Choose between option_1 or option_2 input format')
    #group = parser.add_mutually_exclusive_group()
    #parser_1 = subparsers.add_parser('option_1', help='Provide raw protein sequences, and perform blast before preparation for SVM')
    #parser_2 = subparsers.add_parser('option_2', help='Provide blast results as an input, directly start the preparation for SVM')
    parser.add_argument('-o',
                        '--out_dir',
                        type=str,
                        default='tmp',
                        help='The output directory path')
    parser.add_argument('-r',
                        '--refseq_dir',
                        type=str,
                        help='The blast refseq directory path')
    parser.add_argument('-m',
                        '--merged_dir',
                        type=str,
                        default='',
                        help='The merged blast output directory path. By default, it is assumed that the merged blast outputs sit in the same directory \
                        as the fasta files')   
    parser.add_argument('-b',
                        '--min_bit',
                        type=int,
                        default=MIN_BIT_SCORE,
                        help='The minimum bitscore to filter blast hits.') 
    parser.add_argument('-d',
                        '--bit_delta',
                        type=int,
                        default=BIT_DELTA,
                        help='The Delta threshold for bitscore difference between classified groups') 
    parser.add_argument('-e',
                        '--max_e_val',
                        type=float,
                        default=MAX_E_VALUE,
                        help='The maximum e-value threshold to filter blast hits.') 
    parser.add_argument('-s',
                        '--seq_size',
                        type=int,
                        default=0,
                        help='The minimum sequence length to be parsed into the training set.') 
    args = parser.parse_args()
    return args

def main():
    args = mainArgs()
    
    ## Configuration for .log ##########################
    logName = os.path.basename(args.input) + '_' + PARENT_NAME + '_'+ args.out_dir +'_classifyBlast.log'
    logging.basicConfig(level=logging.DEBUG, format=("%(asctime)s - %(funcName)s - %(message)s"), filename=logName, filemode='w')
    
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(funcName)s - %(message)s")
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)    
    ##################################################
    
    ## Defining Path suffix Variables #################
    refseq_prefix = '_blastx_refseq_prot.txt'
    refseq_path   = os.path.basename(args.input) + refseq_prefix
    refseq_path   = os.path.join(args.refseq_dir, refseq_path)
    merged_prefix = '_merged_blast.txt'
    merged_path   = os.path.basename(args.input) + merged_prefix
    dmp_path = GI_PROT_DMP
    parent_dir = PARENT_LIST
    
    if not os.path.exists(args.out_dir):
        logging.info("Output direcotry does not exist. Creating the output directory..")
        os.makedirs(args.out_dir)
        logging.info("Output directory successfully created")
    
    if args.merged_dir:
        merged_path   = os.path.join(args.merged_dir, merged_path)
    else:        
        merged_path   = os.path.join(os.getcwd(), merged_path)    

    ## Parsing of Blast Output, COnversion of gi to taxid for member lookup (Classification) and Sorting via bitscore ## 
    logging.info("loading parent list")
    parent_list = loadFamily(parent_dir)
    filtered_querries = parseBlast(merged_path)
    filtered_querries = filterBestBit(args, filtered_querries) #Reduce memory usage 
    querries = parseBlast(refseq_path)
    querries = importBlastDict(filtered_querries, querries)  
    #### NOTE: othersTrainClassification, othersBlastClassification are both taxa-ID dictionaries for non-parent querries #### 
    gi_list = giReference(dmp_path)
    trainingClass, blastClass, othersTrainClass, othersBlastClass = classifyFromBlast(querries, args, gi_list, PARENT_TAXA_ID, PARENT_NAME, parent_list)
    print "length of othersTrainClass dict printing.."
    print len(othersTrainClass)
    print "length of othersBlastClass dict printing.."
    print len(othersBlastClass)
    
    logging.debug("classification completed, writing to fasta. classified query size is %d" %len(trainingClass))
    logging.info("Setting: minimum sequence length for parsing is set to be %d" %args.seq_size)
    short_seq_dict = seqSplit(args, trainingClass, blastClass)
    logging.info('Main Classification process completed.')
    ## Preparation of classification summary for the non-Parent hits ###########
    other_dict = prepareOtherTaxaDict(OTHERS_SPEC_LIST)
    prepareOtherSummary(args, other_dict, othersTrainClass, othersBlastClass, short_seq_dict) 
    logging.info('Sub-classification summary for other hits have been successfully completed.')
    logging.info('Exiting...')
    
if __name__ == '__main__': 
    main()
# ENABLING THIS GIVES : /etc/bashrc: line 106: /opt/python/bin/virtualenvwrapper.sh: No such file or directory

