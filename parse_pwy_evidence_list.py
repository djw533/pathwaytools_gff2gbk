
import os

def parseArgs():
    """Parse command line arguments"""

    import argparse

    try:
        parser = argparse.ArgumentParser(
            description='Parse name matching data from pathway tools into a csv file for loading into R')

        parser.add_argument('-i',
            '--input',
            action='store',
            required = True,
            help="input directory that's run the full program")
        parser.add_argument('-o',
            '--output',
            action='store',
            default = "total_pwys_and_reactions_per_genome_output.csv",
            help='Output file')

    except:
        print("An exception occurred with argument parsing. Check your provided options.")
        traceback.print_exc()

    return parser.parse_args()



def parse_pwy_evidence_list(input):

    data = open(input)
    lines = data.readlines()



    output_dict = {}


    for line in lines:
        if line.startswith(";"):
            continue
        else:
            toks = line.strip().strip("()").split() # remove the newlines, any parentheses, and any quotes

            pathway = toks[0].replace('"','') # set the pathway Name
            reactions = 'c({reactions})'.format(reactions=",".join(toks[1:]))

            #add to dict
            output_dict[pathway] = reactions


    return(output_dict)



def main():

    args = parseArgs()

    #create the output file

    with open(args.output, "w") as output:
        output.write("strain,pathway,reactions\n")

    #1 - parse through each file

    for root, dirs, files in os.walk(args.input):

        if len(dirs) > 0 & dirs[0] == "1.0":
            strain = root.split('/')[-1] # should be the name of the folder

        if root.endswith("/1.0/reports"): # i.e. is the base dir with the files in it:
            input_file = "{root}/pwy-evidence-list.dat".format(root = root)

            # now parse this

            pwy_evidence = parse_pwy_evidence_list(input_file)

            #write this out:
            with open(args.output,"a") as output:
                for pathway, reactions in pwy_evidence.items():
                    output.write("{strain},{pathway},{reactions}\n".format(
                    strain = strain, pathway = pathway, reactions = reactions
                    ))






if __name__ == '__main__':
    main()
