PSYTRANS EXTENSION: Mollusca Classification Project
========

Extension of Psytrans to classify mollusca genes in contaminated biological sample

CLASSIFYING BLAST OUTPUTS INTO TRAINING SETS
SCRIPT NAME
       classify_blast_hits.py 

SYNOPSIS
       python classify_blast_hits.py [INPUT_FASTA] [-O OUTPUT_DIR][-R REFSEQ_DIR] [-M MERGED_DIR] [OPTIONS]

DESCRIPTION
       classify_blast_hits.py uses blast hits statistics to classify whether queried sequences belong to a specific family of species.
       The script takes as input the queried fasta file and the directory of its relevant blast results.
       The queries will be compared to these two files using BLASTX.
       Alternatively, the user can provide the output of pre-computed BLASTX searches (in tabular format: -outfmt 6 or 7).
       The classification process will then produce two version of training sets(.fasta) to be used for psytrans applications: 1) A general less-stringent training set 2) A more stringent training set based on bitscore differences .


OPTIONS
   Generic Program Information

       -h, --help
              Print a usage message briefly summarizing the command-line options.

EXAMPLE
       The user is required to have already computed and stored the blast results of the queried fasta sequence of interest. 
       The blast results are assumed to be given separately (as input) for the refseq_prot database and the taxanomically related reference database.
       The user will have to specify the relevant blast parameters in performing the classification and provide the path to the blast results.
       The script can be executed as follows:

              python classify_blast_hits.py [INPUT_FASTA] [-O OUTPUT_DIR][-R REFSEQ_DIR] [-M MERGED_DIR] [OPTIONS]

       The creation fo the GI to taxid conversion table will take around 5-10mins. Whereas the blast classification process take no more than 20 mins to complete (based on 70k queries). When the script finishes, the training fasta sequences will be created in the output directory along with summaries of the breakdown of identifiable groups of species. 

-----------------------------------------------------------------------------------------------------------------

OBTAINING DESCENDANT OF A SPECIFIC PARENT SPECIES
SCRIPT NAME
       get_descendant.py 

SYNOPSIS
       python get_descendant.py [species_taxid] [-O OUTPUT_FILENAME][OPTION]

DESCRIPTION
       get_descendant.py parses the NCBI taxanomy tree into dictionaries and obtain the relevant descendants of the species family in the form of a list.
       The script takes as input the taxa ID of the species of interest. path to the nodes.dmp and names.dmp files will have to correctly defined within the script.
       The script then produces a list of the taxa_id of all descendant of the queried species which will be stored in a text(.txt) file. This file is an essential input for the classification process in classify_blast_hits.py.  


OPTIONS
   Generic Program Information

       -h, --help
              Print a usage message briefly summarizing the command-line options.

EXAMPLE
       The user is required to have read-access to the NCBI_taxanomy names.dmp and nodes.dmp files. 
       The taxa_id for the group of species of interest will also have to be obtained.
       The script can be executed as follows(take for example the group of species of interest are 'mollusca'):

              python get_descendant.py '6447' -o mollusca_list.txt

       The script will take no more than 5mins to run. The list of descendant on the species family of interest will be written into the output_file with each descendant's taxa_id corresponding one line in the .txt file.

------------------------------------------------------------------------------------------------------------------



AUTHOR
       Written by Jue-Sheng Ong.

REPORTING BUGS
       Report bugs at sylvain.foret@anu.edu.au
       Psytrans repository <https://github.com/sylvainforet/psytrans>

COPYRIGHT
       Copyright © 2014 Sylvain Forêt & Jue-Sheng Ong.

       psytrans  is a free  software and comes with ABSOLUTELY NO WARRANTY.  You are welcome to redistribute it under the terms of the GNU General Public License
       versions 3 or later.  For more information about these matters see http://www.gnu.org/licenses/.
