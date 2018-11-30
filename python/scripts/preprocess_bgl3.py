## Extract information from Bgl3 deep mutational scan
## We only extract single and double mutations. 
## For each mutation we also count how many times it appears in each file
import os
import functools
import pandas as pd

data_dir = "data"

input_dir = os.path.join(data_dir, "seq-fcn-data/Bgl3_random_mutagenesis")
output_csv = os.path.join(data_dir, "Bgl3_single_mutations.csv")

files = {'unsorted':'unlabeled_Bgl3_mutations.txt',
         'positive':'positive_Bgl3_mutations.txt',
         'positive_hitemp':'positive_Bgl3_mutations_hitemp.txt'}

def extract_mutations(filename, comma_count, directory=input_dir):
    "Return all lines without a comma in it"
    with open(os.path.join(directory, filename)) as f:
        return [line.strip() for line in f if line.count(",") == comma_count]

def extract_single_mutations(filename, *args, **kw):
    return extract_mutations(filename, 0, *args, **kw)

def extract_double_mutations(filename, *args, **kw):
    return extract_mutations(filename, 1, *args, **kw)

def convert_mutation_list_to_df(mutations_list, count_label):
    """ Take in a list of mutations and count how many of there in a dataframe
        The number of counts appears in the count_label column
    """
    df = pd.Series(mutations_list)
    unique_counts = df.groupby(df).nunique()
    retdf = pd.DataFrame( {"mutations": unique_counts.index, 
                          count_label:unique_counts} )
     # the mutations become the index of the dataframe. Lets reset that
    retdf.reset_index(drop=True, inplace=True)
    return(retdf)

# Create pandas dataframes for each file which counts the numbers of mutations
single_mutations = {k:convert_mutation_list_to_df(
    extract_single_mutations(files[k]), 
    count_label=k) for k in files.keys()}

# Now lets merge all the dataframes together
df = functools.reduce(lambda left, right: pd.merge(left, right, on="mutations"), 
                      single_mutations.values())

# Save it as csv
df.to_csv(output_csv)
