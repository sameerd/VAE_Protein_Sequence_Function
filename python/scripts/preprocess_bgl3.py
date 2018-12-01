## Extract information from Bgl3 deep mutational scan
## We only extract single and double mutations. 
## For each mutation we also count how many times it appears in each file
import os
import functools
import pandas as pd

data_dir = "data"
working_dir = os.path.join(data_dir, "working")

input_dir = os.path.join(data_dir, "seq-fcn-data/Bgl3_random_mutagenesis")
single_output_file = os.path.join(working_dir, "Bgl3_single_mutations.tsv")
double_output_file = os.path.join(working_dir, "Bgl3_double_mutations.tsv")

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
        This works for both single and double mutations
    """
    df = pd.Series(mutations_list)
    unique_counts = df.groupby(df).nunique()
    retdf = pd.DataFrame( {"mutations": unique_counts.index, 
                          count_label:unique_counts} )
     # the mutations become the index of the dataframe. Lets reset that
    retdf.reset_index(drop=True, inplace=True)
    return(retdf)


## Single Mutations
# Create dataframes for each file which counts the numbers of mutations
single_mutations = {k:convert_mutation_list_to_df(
    extract_single_mutations(files[k]), 
    count_label=k) for k in files.keys()}

# Now lets merge all the dataframes together
single_df = functools.reduce(lambda left, right: 
        pd.merge(left, right, on="mutations", how="outer"),
        single_mutations.values())

# Save it as csv
single_df.to_csv(single_output_file, sep="\t")

## Double Mutations
double_mutations = {k:convert_mutation_list_to_df(
    extract_double_mutations(files[k]), 
    count_label=k) for k in files.keys()}


# Check that all the double mutations are indexed correctly
# Looks like this always passes so there is no need to always check this
for key, df in double_mutations.items():
  # Extract the first and the second index where the mutations occur
  y = df["mutations"].str.extractall(
          r"\D(?P<index1>[\d]*)\D,\D(?P<index2>[\d]*)\D")
  y.index1 = pd.to_numeric(y.index1)
  y.index2 = pd.to_numeric(y.index2)
  # check that the first mutation is never larger than or equal 
  # to the second
  if y.loc[y.index1 >= y.index2].size != 0:
    print("The ordering of the double mutations in ", key, 
          " is not correct!")

# Now lets merge all the dataframes together
double_df = functools.reduce(lambda left, right: 
        pd.merge(left, right, on="mutations", how="outer"),
        double_mutations.values())

double_df.to_csv(double_output_file, sep="\t")


