#!/usr/bin/env python3


import os
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import translate
from Bio.Seq import reverse_complement
from datetime import date


###### parse arguments

def parseArgs():
    """Parse command line arguments"""

    import argparse

    try:
        parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter ,
            description="Convert gffs used as input to panaroo into pathologic file format that can be processed for pathway tools PGDB (using mpwt)")
        parser.add_argument('-g',
			    '--gff',
			    action='store',
                nargs='+',
                required=True,
			    help='Gff file(s) to search <required>')
        parser.add_argument('-e',
			    '--eggnog',
			    action='store',
                required=True,
			    help='EGGNOG csv file for each of the protein sequence groups from panaroo')
        parser.add_argument('-p',
			    '--panaroo',
			    action='store',
                required=True,
			    help='Panaroo gene_presence_absence.csv file')
        parser.add_argument('-s',
			    '--species',
			    action='store',
			    help='Species information as csv file - strain,genus_species name (with header)')
        parser.add_argument('-o',
			    '--output_dir',
			    action='store',
                default="formatted_gbks",
			    help='Output directory, default=[formatted_gbks]. Will not write over a previously existing output folder!')

    except:
        print("An exception occurred with argument parsing. Check your provided options.")
        traceback.print_exc()

    return parser.parse_args()


def __read_gene_presence_absence__(input_file):

    ## read in the gene_presence_absence.csv file


    data = open(input_file)
    lines = data.readlines()
    data.close()


    gene_groups = {}
    # making a dictionary with dictionary in it - key is group name - first field (delimited using comma)
    # dicts within are strains (value is list of strains), and description (value is string of the description - (first field (; deliminated) in the 3rd fiel (, deliminated)) )

    for line in lines:
        toks = line.strip().split(',')

        ##get the unique gene/group name
        group = toks[0]

        ## get the description:
        description = toks[2].split(';')[0]

        ## now get the strains:
        strains  = toks[3:]
        ### remove empty items in the list:
        strains = list(filter(None, strains))



        ## add to the dictionary:

        gene_groups[group] = {}
        gene_groups[group]["strains"] = strains
        gene_groups[group]["description"] = description


    return(gene_groups)


def __genes_to_gene_groups__(input_file):
    ## read in the gene_presence_absence.csv file


    data = open(input_file)
    lines = data.readlines()
    data.close()


    genes = {}
    # making a dictionary with dictionary in it - key is group name - first field (delimited using comma)
    # dicts within are strains (value is list of strains), and description (value is string of the description - (first field (; deliminated) in the 3rd fiel (, deliminated)) )

    for line in lines:
        toks = line.strip().split(',')

        ##get the unique gene/group name
        group = toks[0]

        ## get the description:
        description = toks[2].split(';')[0]

        ## now get the strains:
        gene_ids  = toks[3:]
        ### remove empty items in the list:
        gene_ids = list(filter(None, gene_ids))

        for gene_id in gene_ids:
            ## add to the dictionary:

            genes[gene_id] = {}
            genes[gene_id]["group"] = group
            genes[gene_id]["description"] = description


    return(genes)


def __read_eggnog_data__(input_file):

    data = open(input_file)
    lines = data.readlines()
    data.close()

    eggnog_dict = {}

    for line in lines:

        toks = line.strip().split(",")

        id = toks[0]
        pref_name = toks[5]
        GO_terms = toks[6].split(';')
        EC_terms = toks[7].split(';')
        function = toks[20]

        eggnog_dict[id] = {}
        eggnog_dict[id]["pref_name"] = pref_name
        eggnog_dict[id]["GO_terms"] = GO_terms
        eggnog_dict[id]["EC_terms"] = EC_terms
        eggnog_dict[id]["function"] = function

    return(eggnog_dict)


def __gff_to_patho__(gff,panaroo_dict,eggnog_dict,output_dir):

    ## 1. read in gff and get the annotation section

    data = open(gff)
    lines = data.readlines()
    data.close()

    strain = gff.split('/')[-1].split(".gff")[0].replace(".","_")

    os.makedirs("{output_dir}/{strain}".format(output_dir=output_dir,strain=strain))


    # = open("{output_dir}/{strain}/{strain}.pf".format(output_dir=output_dir,strain=strain),"a") # set the output file writing

    ## 2. run through contig by contig

    #set current contig :
    current_contig = 'contig'


    for line in lines:
        if line == "##FASTA\n":
            return # reached the end of annotation
        elif line.startswith("#"):
            continue
        else: # read and format the new file:
            toks = line.strip().split("\t")
            # make sure is a CDS:
            if toks[2] == "CDS":
                #get details from each line:
                contig = toks[0].replace("|","_").replace(".","_")
                # check to see if new contig and create new file:
                if contig != current_contig:
                    ### change pipe to underscore:
                    contig = contig.replace("|","_").replace(".","_")
                    ###change name back
                    current_contig = contig

                    with open("{output_dir}/{strain}/{contig}.pf".format(output_dir=output_dir,contig=contig,strain=strain),"w") as output:
                        output.write(";;;;;;;;;;;;;;;;;;;;;;;;;\n;; {contig}\n;;;;;;;;;;;;;;;;;;;;;;;;;\n".format(contig=contig))

                #also change the name here as well if required
                contig = contig.replace("|","_").replace(".","_")

                start = toks[3]
                end = toks[4]
                #print(line)
                gene_id = toks[8].split("ID=")[1].split(';')[0]


                ## now find the gene id in the panaroo dictionary to get the group name:
                if gene_id in panaroo_dict:
                    group = panaroo_dict[gene_id]["group"]
                    description = panaroo_dict[gene_id]["description"]
                    # print(gene_id)
                    # print(group)

                    ### now find the information for the group from the eggnog_dict:

                    if group in eggnog_dict:
                        pref_name = eggnog_dict[group]["pref_name"]
                        GO_terms = eggnog_dict[group]["GO_terms"]
                        EC_terms = eggnog_dict[group]["EC_terms"]
                        function = eggnog_dict[group]["function"]

                        #write out:
                        with open("{output_dir}/{strain}/{contig}.pf".format(output_dir=output_dir,strain=strain,contig=contig),"a") as output:
                            output.write("ID\t{gene_id}\n".format(gene_id=gene_id))
                            #pref name
                            if pref_name != "NA":
                                output.write("NAME\t{pref_name}\n".format(pref_name=pref_name))
                            #start and stop
                            output.write("STARTBASE\t{start}\n".format(start=start))
                            output.write("ENDBASE\t{end}\n".format(end=end))
                            #function
                            output.write("FUNCTION\t{function}\n".format(function=function))
                            #EGGNOG and EC terms:
                            for GO in GO_terms:
                                output.write("DBLINK\t{GO}\n".format(GO=GO))
                            for EC in EC_terms:
                                if EC != "NA":
                                    output.write("EC\t{EC}\n".format(EC=EC))
                            #product type
                            output.write("PRODUCT-TYPE\tP\n")
                            # write the // at the end:
                            output.write("//\n")


                        # print(pref_name)
                        # print(GO_terms)
                        # print(EC_terms)
                        # print(function)

                    # else just use the description that panaroo uses.....
                    else:
                        #print("group {group} not here".format(group = group))
                        #print(description)
                        with open("{output_dir}/{strain}/{contig}.pf".format(output_dir=output_dir,strain=strain,contig=contig),"a") as output:
                            output.write("ID\t{gene_id}\n".format(gene_id=gene_id))
                            output.write("STARTBASE\t{start}\n".format(start=start))
                            output.write("ENDBASE\t{end}\n".format(end=end))
                            output.write("FUNCTION\t{function}\n".format(function=description))
                            output.write("PRODUCT-TYPE\tP\n")
                            # write the // at the end:
                            output.write("//\n")


def __gff3_fasta_to_singlefastas__(gff_file,output_dir):

        strain = gff_file.split('/')[-1].split(".gff")[0].replace(".","_")


        fasta_file_sequence = []

        with open(gff_file) as f:
            fasta = False
            for line in f:
                if fasta:
                    fasta_file_sequence.append(line)
                    continue
                if line.startswith("##FASTA"):
                    fasta = True
                    continue
                if fasta == False:
                    continue


        #sequences = ''.join(fasta_file_sequence)

        entries_list = (''.join(fasta_file_sequence)).split('>')
        #print(entries_list)

        for entry in entries_list:
             contig = entry.split('\n')[0].replace("|","_").replace(".","_")
             with open("{output_dir}/{strain}/{contig}.fasta".format(output_dir=output_dir,strain=strain,contig=contig),"w") as output:
                  output.write('>'+entry)


def __gff_to_fasta__(gff_file, type):
    '''convert a gff file with the appended FASTA to
    gff_file = input gff file
    fasta_file = output file
    output: a protein coding FASTA file '''
    sequences = {}
    contigs = {}
    out = open("tmp_fasta.fa", "w")
    with open(gff_file) as f:
        fasta = False
        for line in f:
            if fasta:
                out.write(line)
                continue
            if line.startswith("##FASTA"):
                fasta = True
                continue
            if line.startswith("#"):
                continue
            toks = line.strip().split("\t")
            if toks[2] != "CDS":
                continue
            name = toks[-1].split("ID=")[1].split(';')[0]
            if toks[0] not in contigs:
                contigs[toks[0]] = []
            contigs[toks[0]].append({"name": name, "start": int(toks[3])-1, "stop": int(toks[4]), "strand": toks[6]})
    out.close()

    ## read the contigs and save the final fasta file
    with open("tmp_fasta.fa") as handle:
        for values in SimpleFastaParser(handle):
            curr_contig = values[0].split()[0]
            #print contigs.keys()
            if curr_contig not in contigs: # no CDSs in this contig
                continue
            #print 'still continuing'
            for cds in contigs[curr_contig]:
                #out.write(">" + cds["name"] + "\n")
                seq = values[1][cds["start"]:cds["stop"]]
                if cds["strand"] == "-":
                    seq = reverse_complement(seq)
                if type == 'prot':
                    sequences[cds["name"]] = translate(seq)
                    #out.write(translate(seq) + "\n")
                if type == 'nuc':
                    sequences[cds["name"]] = seq
                    #out.write(seq + "\n")
    os.remove("tmp_fasta.fa")
    return sequences



def __gff_to_gbk__(gff,panaroo_dict,eggnog_dict,output_dir,species_file,protein_sequences):

    data = open(gff)
    lines = data.readlines()
    data.close()

    strain = gff.split('/')[-1].split(".gff")[0].replace(".","_")


    # get the species information:
    species_data = open(species_file)
    species_lines = species_data.readlines()
    species_data.close()


    ### get the species now:
    for s in species_lines:
        toks = s.strip().split(',')
        s_strain = toks[0]
        species = toks[1]

        if strain == s_strain:
            break # now have the species - break from the loop

    print(species)


    os.makedirs("{output_dir}/{strain}".format(output_dir=output_dir,strain=strain))

    # get the gff information into a dictionary - each contig will be a new dictionary:


    contigs = {}

    fasta = False
    fasta_contig = ''
    contig = ''

    for line in lines:
        if fasta:
            #print(contigs.keys())
            ### add the sequence to the dictionary for that specific contig
            if line.startswith('>'):
                # is a new contig:
                fasta_contig = line.strip().split('>')[1]
                ## create a new entry in the dictionary (if it exists - otherwise make it):
                if fasta_contig in contigs:
                    contigs[fasta_contig]["sequence"] = []
                else: # not any features in the contig:
                    contigs[fasta_contig] = {}
                    contigs[fasta_contig]["sequence"] = []

            else:
                contigs[fasta_contig]["sequence"].append(line.strip())
                ### now have added the nuc sequence to the dictionary.
            continue
        if line.startswith("##FASTA"):
            fasta = True
            continue
        if line.startswith("#"):
            continue

        toks = line.strip().split("\t")
        if toks[2] != "CDS":
            continue

        contig = toks[0]
        if contig not in contigs:
            contigs[contig] = {} # set new dictionary for this contig if it doesn't already exist:

        cds_name = toks[8].split("ID=")[1].split(";")[0]
        start = toks[3]
        end=toks[4]
        direction=toks[6]
        #create new dictionary for this cds:
        contigs[contig][cds_name] = {}
        contigs[contig][cds_name]["start"] = start
        contigs[contig][cds_name]["end"] = end
        contigs[contig][cds_name]["direction"] = direction




        ### add the details from the panaroo dictionary and the eggnog dictionary:
        ## now find the gene id in the panaroo dictionary to get the group name:
        cds_name = cds_name
        if cds_name in panaroo_dict:
            ##set the name for group:
            group = panaroo_dict[cds_name]["group"]
            ##now add them to the contigs dict
            contigs[contig][cds_name]["group"] = panaroo_dict[cds_name]["group"]
            contigs[contig][cds_name]["description"] = panaroo_dict[cds_name]["description"]

            ### now find the information for the group from the eggnog_dict:

            if group in eggnog_dict:
                contigs[contig][cds_name]["pref_name"] = eggnog_dict[group]["pref_name"]
                contigs[contig][cds_name]["GO_terms"] = eggnog_dict[group]["GO_terms"]
                contigs[contig][cds_name]["EC_terms"] = eggnog_dict[group]["EC_terms"]
                contigs[contig][cds_name]["function"] = eggnog_dict[group]["function"]


    ## now write out this data into the desired genbank that mpwt/pathway tools wants:

    ## write out per contig:
    open("{output_dir}/{strain}/{strain}.gbff".format(output_dir=output_dir,strain=strain),"w")


    ## move contig by contig:

    f = open("{output_dir}/{strain}/{strain}.gbff".format(output_dir=output_dir,strain=strain),"a")
    for contig in contigs.keys():
        ## get the length of the contig:
        formatted_contig =  contig.replace("|","_").replace(".","_")
        contig_length = len(''.join(contigs[contig]["sequence"]))
        #get the date:
        today = date.today()
        d1 = today.strftime("%d-%b-%Y")
        ## write this line out:
        f.write("LOCUS\t{contig}\t{length} bp\tDNA\tlinear\tBCT\t{date}\n".format(contig=formatted_contig,length=str(contig_length),date=d1))
        ##DEFINITION:
        f.write("DEFINITION\t{species} {strain} DNA, complete genome\n".format(species=species.replace("_"," "),strain=strain))
        ##ACCESSION
        f.write("ACCESSION\t{contig}\n".format(contig=formatted_contig))
        ##VERSION
        f.write("VERSION\t{contig}\n".format(contig=formatted_contig))
        ##KEYWORDS
        f.write("KEYWORDS\t.\n")
        ##SOURCE
        f.write("SOURCE\t.\n")
        ##ORGANISM
        f.write("ORGANISM\t{species} {strain}".format(species=species.replace("_"," "),strain=strain))
        ## the rest:
        f.write("\tBacteria; Proteobacteria; Gammaproteobacteria; Enterobacterales;\t")
        f.write("\tEnterobacteriaceae; Serratia.\n")
        ## FEATURES:
        f.write("FEATURES\t\tLocation/Qualifiers\n")
        ##write the source
        f.write("\tsource\t1..{length}\n".format(length=str(contig_length)))
        f.write('\t\t/organism="{species}"\n\t\t/mol_type="genomic_DNA"\n\t\t/strain="{strain}"'.format(species=species.replace("_"," "),strain=strain))



        ### now move through CDS by CDS:
        for cds, cds_items in contigs[contig].items():
            if cds == "sequence":
                continue

            ### now write parts out for each cds feature:
            #first check if forward or reverse:
            if cds_items["direction"] == "+":
                f.write('\tgene\t{start}..{end}\n\t\t/locus_tag="{cds}"\n'.format(start=cds_items["start"],end=cds_items["end"],cds=cds))
                f.write('\tmRNA\t{start}..{end}\n\t\t/locus_tag="{cds}"\n'.format(start=cds_items["start"],end=cds_items["end"],cds=cds))
                f.write('\tCDS\t{start}..{end}\n\t\t/locus_tag="{cds}"\n'.format(start=cds_items["start"],end=cds_items["end"],cds=cds))
            elif cds_items["direction"] == "-":
                f.write('\tgene\tcomplement({start}..{end})\n\t\t/locus_tag="{cds}"\n'.format(start=cds_items["start"],end=cds_items["end"],cds=cds))
                f.write('\tmRNA\tcomplement({start}..{end})\n\t\t/locus_tag="{cds}"\n'.format(start=cds_items["start"],end=cds_items["end"],cds=cds))
                f.write('\tCDS\tcomplement({start}..{end})\n\t\t/locus_tag="{cds}"\n'.format(start=cds_items["start"],end=cds_items["end"],cds=cds))

            ## now write the GO terms and others:
            #write the product first:
            f.write('\tproduct')
















    f.close()




def main():


    args = parseArgs()

    os.makedirs(args.output_dir)



    #Order of things to do:

    #1. Read in the panaroo data frame and create a dictionary of the "gene group" against the names of the "genes" in the gff file - also need gene decriptions


    #dict = __read_gene_presence_absence__(args.panaroo)

    panaroo_dict = __genes_to_gene_groups__(args.panaroo)

    #2. Then read in the eggnog data and create a dictionary with the "gene_group" with the different GO terms that correspond with the genes

    eggnog_dict = __read_eggnog_data__(args.eggnog)



    for gff in args.gff:
        #protein_sequences = __gff_to_fasta__(gff,"prot")
        #__gff_to_gbk__(gff,panaroo_dict,eggnog_dict,args.output_dir,args.species,protein_sequences)
        __gff_to_patho__(gff,panaroo_dict,eggnog_dict,args.output_dir)
        __gff3_fasta_to_singlefastas__(gff,args.output_dir)

    #3. Use this data with each gff file to convert them to genbanks with the same information for each grouped gene








if __name__ == '__main__':
    main()
