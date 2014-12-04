#!/usr/bin/python

import os
import sys
import argparse



##PLEASE CHANGE PATH ACCORDINGLY

###These files can be obtained form the NCBI taxanomy website
NAMES_PATH = "/home/juesheng.ong/NCBI_taxa/names.dmp"
NODES_PATH = "/home/juesheng.ong/NCBI_taxa/nodes.dmp"



class Nodes:
    """Nodes"""
    def __init__(self):
        self.tax_id = 0       # Number of the tax id.
        self.parent = 0       # Number of the parent of this node
        self.children = []    # List of the children of this node
        self.tip = 0          # Tip=1 if it's a terminal node, 0 if not.
        self.name = ""        # Name of the node: taxa if it's a terminal node, numero if not.       
    def genealogy(self):      # Trace genealogy from root to leaf
        ancestors = []        # Initialise the list of all nodes from root to leaf.
        tax_id = self.tax_id  # Define leaf
        while 1:
            if name_object.has_key(tax_id):
                ancestors.append(tax_id)
                tax_id = name_object[tax_id].parent
            else:
                break
            if tax_id == "1":
                # If it is the root, we reached the end.
                # Add it to the list and break the loop
                ancestors.append(tax_id)
                break
        return ancestors # Return the list



#############################
#                           #
#   Read taxonomy files     #
#                           #
#############################


def read_taxa(names_path, entry_type):
    name_dict = {}          # Initialise dictionary with TAX_ID:NAME
    name_dict_reverse = {}  # Initialise dictionary with NAME:TAX_ID

    # Load  NCBI names file ("names.dmp")
    name_file =  open(names_path,"r")
    while 1:
        line = name_file.readline()
        if line == "":
            break
        line = line.rstrip()
        line = line.replace("\t","")
        tab = line.split("|")
        if tab[3] == entry_type:
            tax_id, name = tab[0], tab[1]     # Assign tax_id and name ...
            name_dict[tax_id] = name          # ... and load them
            name_dict_reverse[name] = tax_id  # ... into dictionaries
    name_file.close()
    return name_dict, name_dict_reverse

######################
# Load taxonomy  ####
#####################

def load_taxa(nodes_path, name_dict):
    
    # Define taxonomy variable
    name_object = {}
    # get nodes.dmp handle
    taxonomy_file = open(nodes_path,"r")
    while 1:
        line = taxonomy_file.readline()
        if line == "":
            break
        #print line
        line = line.replace("\t","")
        tab = line.split("|")
        
        tax_id = str(tab[0])
        tax_id_parent = str(tab[1])
        division = str(tab[4])

        # Define name of the taxid
        name = "unknown" #initialise name
        if tax_id in name_dict:
            name = name_dict[tax_id]
        
        if not name_object.has_key(tax_id):
            name_object[tax_id] = Nodes()
        name_object[tax_id].tax_id   = tax_id        # Assign tax_id
        name_object[tax_id].parent   = tax_id_parent # Assign tax_id parent
        name_object[tax_id].name     = name          # Assign name
        
        if  tax_id_parent in name_object:
            children = name_object[tax_id_parent].children  # If parent is is already in the object
            children.append(tax_id)                  # ...we found its children.
            name_object[tax_id_parent].children = children  # ... so add them to the parent
    taxonomy_file.close()
    return name_object
    
###TODO Tidy up into neat functions: potentially use Pycogent's module instead ##

def writeToFile(clist, out_name):
    c = 0
    handle = open(out_name,'w')
    for i in clist:
        #sci_name = name_dict[i]
        handle.write('%s\n'%i)
        c += 1
        if c % 2000 == 0:
            print "parsed %d species entries to %s" % (c, out_name)
    handle.close()


def readLineage(args, tax_id, name_dict, nodes_path):
    """read the lineage of the dict_object, and return the list of
    organisms that have ancestor:[tax_id] in them."""
    global name_object #export dictionary to global to use the geneology() function
    name_object = load_taxa(nodes_path, name_dict)
    parent_id = str(tax_id)
    descendant_list = []
    print "name_object size is: %d" %len(name_object)
    if args.list_immed:
        print "testing: immediate parent"
        print name_object[args.tax_id].parent
        print "testing: immediate children"
        print name_object[args.tax_id].children
    for entry in name_object:
        query_id = str(name_object[entry].tax_id)     
        if parent_id in name_object[query_id].genealogy():
            descendant_list.append(query_id)
            if len(descendant_list) % 1000 == 0:
                print "parent found: %d"%len(descendant_list)
    return descendant_list


#######################################            
#### RUNNING FUNCTIONS ################
#######################################



def mainArgs():
    """Process command-line arguments"""

    parser = argparse.ArgumentParser(description='Perform SVM Classification of Host and Symbiont (or Parasite) Sequences')
    parser.add_argument('tax_id',
                        help='The NCBI taxa_id of the queried species')
    #subparsers = parser.add_subparsers(help='Choose between option_1 or option_2 input format')
    #group = parser.add_mutually_exclusive_group()
    #parser_1 = subparsers.add_parser('option_1', help='Provide raw protein sequences, and perform blast before preparation for SVM')
    #parser_2 = subparsers.add_parser('option_2', help='Provide blast results as an input, directly start the preparation for SVM')
    parser.add_argument('-o',
                        '--out_name',
                        type=str,
                        default='desc_list.txt',
                        help='The output filename of the list containing all descendant of the queried species')                 
    parser.add_argument('-a',
                        '--list_immed',
                        type=bool,
                        default=False,
                        help='print immediate parent and children of the species')
    args = parser.parse_args()
    return args



def main():
    args = mainArgs()
    #Node = Nodes()
    name_dict, name_dict_reverse = read_taxa(NAMES_PATH, "scientific name")
    desc_list = readLineage(args, args.tax_id, name_dict, NODES_PATH)
    writeToFile(desc_list, args.out_name)
    print "Parsing of descendant id entries successfully completed."

if __name__ == '__main__':
    main()
