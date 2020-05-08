#! /home/djwilliams/software/miniconda3/bin/python


#### read in the pathway inference data and return a new file in csv with the details wanted


def parseArgs():
    """Parse command line arguments"""

    import argparse

    try:
        parser = argparse.ArgumentParser(
            description='Parse name matching data from pathway tools into a csv file for loading into R')

        parser.add_argument('-i',
            '--input',
            action='store',
            help='input file')
        parser.add_argument('-o',
            '--output',
            action='store',
            help='Output file')

    except:
        print("An exception occurred with argument parsing. Check your provided options.")
        traceback.print_exc()

    return parser.parse_args()




def __parse_name_data__(input_file,output_file):

    data = open(input_file)
    lines = data.readlines()
    data.close()


    unambiguous = False
    ambiguous = False

    genes = {}

    for line in lines:
        if line.startswith("-----------"):
            continue

        if line.startswith("Unambiguous reaction assignments:"):
            unambiguous = True
            continue

        if line.startswith("Ambiguous function name matches (no EC# or GO term) that should be examined:"):
            ambiguous = True
            continue


        ### now go through:

        if unambiguous == False:
            continue
        elif unambiguous == True and ambiguous == False:
            ## start parse data of the genes:

            gene_id = line[23:40].strip() # get the geneid
            #print(gene_id)
            reactions = line[107:].strip().split(', ') # and all the reactions
            #print(reactions)

            ##add gene_id to the dictionary if not already in it:
            if gene_id not in genes:
                genes[gene_id] = []
                for reaction in reactions:
                    reaction = reaction.split(' ')[-1] # ignore crap left over if there
                    if reaction not in genes[gene_id]:
                        genes[gene_id].append(reaction.strip()) # add the reaction to the list of reactions that gene is involved in
            elif gene_id in genes:
                for reaction in reactions:
                    reaction = reaction.split(' ')[-1] # ignore crap left over if there
                    if reaction not in genes[gene_id]:
                        genes[gene_id].append(reaction.strip()) # add the reaction to the list of reactions that gene is involved in


        if ambiguous == False:
            continue

        if ambiguous == True:
            continue




    #print(genes)

    #write out these :

    with open(output_file,"w") as output:
        output.write("gene,reaction\n")

    #go through genes:
    for gene, reactions in genes.items():
        for reaction in reactions:
            with open(output_file,"a") as output:
                output.write("{gene},{reaction}\n".format(gene=gene,reaction=reaction))


def main():

    args = parseArgs()

    __parse_name_data__(args.input,args.output)


if __name__ == '__main__':
    main()
