#! /home/djwilliams/software/miniconda3/bin/python


#### read in the pathway inference data and return a new file in csv with the details wanted


def parseArgs():
    """Parse command line arguments"""

    import argparse

    try:
        parser = argparse.ArgumentParser(
            description='Parse pathway inference data from pathway tools into a csv file for loading into R')

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


def __parse_inf_data__(input_file,output_file):

    data = open(input_file)
    lines = data.readlines()
    data.close()

    # create dictionary for the pathways in the file:
    pathways = {}

    # now create a list to put the file into, with the breaks between elements in the list where the line starts with a '(:PATHWAY'
    stripped_lines = []

    for line in lines:
        stripped_lines.append(line.strip('\n'))

    pathways_to_parse = ''.join(stripped_lines).split("(:PATHWAY ")

    ##
    # now go through each and put into dictionaries
    for pathway_details in pathways_to_parse:
        if pathway_details == "":
            continue #ignore the first one - as is nothing.
        toks = pathway_details.split(':')[1:] # take all parts but the name of the pathway
        pathway_name = pathway_details.split(':')[0].strip()

        ### now add the pathway to the dictionary as a new dictionary:
        pathways[pathway_name] = {}

        ## work through toks - with all the details:  - for each there should be a key and value in the list - add these to the dictionary:
        for tok in toks:
            #print(tok)
            key = tok.split(' ')[0]
            value = ' '.join(tok.split(' ')[1:]).strip()



            ## special action required if the key is "EXPLANATION-CODE":
            if key == "EXPLANATION-CODE" and value.startswith("(") and  value.endswith(")"):
                #remove the beginning and end parentheses
                value = value[1:-1]
                #print(value)



                ### get the value for the explanation code:
                if '(' in value:
                    explanation_code = value.split('(')[0].strip()
                    pathways[pathway_name]["EXPLANATION-CODE"] = explanation_code
                else:
                    explanation_code = value.strip()
                    pathways[pathway_name]["EXPLANATION-CODE"] = explanation_code

                #print(explanation_code)

                if explanation_code == "CHECK-FOR-KEY-NON-REACTIONS":
                    value_to_keep = ' '.join(value.split('(')[1:]).strip(')')
                    pathways[pathway_name]["CHECK-FOR-KEY-NON-REACTIONS"] = value_to_keep
                else:

                    # define the open and ends of the parentheses in value, character by character
                    isopen = 0
                    passed_parantheses = False
                    regions = []
                    for i in range(0,len(value),1): # go through numbers of the length of the string:
                        if value[i] == "(":
                            if isopen == 0:
                                start = i
                            passed_parantheses = True
                            isopen += 1
                        if value[i] == ")":
                            isopen -= 1
                            if passed_parantheses == True and isopen == 0:
                                stop = i
                                ##add this to the regions_to_split:
                                regions.append(value[start+1:stop])

                    ##now go through these
                    for region in regions:
                        if len(region.split()) > 1:
                            f = region.split()
                            key = f[0].strip()
                            value = ' '.join(f[1:]).strip().strip('(').strip(')')
                            pathways[pathway_name][key] = value

            else:
                pathways[pathway_name][key] = value.strip(')')



    ### now go through the file

    #print(pathways)

    ## now want to get all the unique names of each key that are in all the dictionaries (these will be headers in the csv file - we don't mind having columns with NA values):

    key_names = []



    for key in pathways.keys():
        items = pathways[key].keys()
        # add to key names if not in it
        for item in items:
            if item not in key_names:
                key_names.append(item)


    # go back through and create new keys if required:
    for key in pathways.keys():
        items = pathways[key].keys()
        # add to key names if not in it
        for name in key_names:
            if name not in items:
                pathways[key][name] = ''



    key_names.sort()

    ## start writing out:
    with open(output_file,"w") as output:
        output.write("PATHWAY,{colnames}\n".format(colnames=','.join(key_names)))

    ## through the dictionary:
    for key in pathways.keys():
        #print(key)
        elements_to_write = []
        for k in key_names:
            elements_to_write.append(pathways[key][k])
        ## replace empty values with NA:
        elements_to_write = ["NA" if x == '' else x for x in elements_to_write]
        #write out:
        #print(elements_to_write)
        with open(output_file,"a") as output:
            output.write("{pathway_name},{values}\n".format(pathway_name=key,values = ','.join(elements_to_write)))






def main():

    args = parseArgs()

    __parse_inf_data__(args.input,args.output)


if __name__ == '__main__':
    main()
