import pandas as pd

def parse_hmmer_domtblout(filename):
    # Reads input line by line
    data = []
    with open(filename, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                data.append(line.split())

    # Converts list to DataFrame
    df2 = pd.DataFrame(data)

    # Define column names
    cols = ["target_name", "target_accession", "tlen", "query_name", "accession", "qlen", "E-value", 
            "score", "bias", "#", "of", "c-Evalue", "i-Evalue", "domain_score", "domain_bias",
            "hmm_start", "hmm_end", "ali_start", "ali_end", "env_from", "env_to", "acc"]
    
    df1 = df2.iloc[:,:len(cols)]
    df1.columns=cols

    # check if df1 empty
    if df1.empty:
        print("The dataframe is empty.")
        df1.columns = cols + ["description of target"]
        return df1 #return empty df

    # Subtracts already read columns to get 'description of target'
    df2[len(df1.columns)] = df2[df2.columns[len(df1.columns):]].apply(
        lambda x: ' '.join(x.dropna().astype(str)), axis=1)
    
    # Combine dataframes
    df = pd.concat([df1, df2[len(df1.columns)]], axis=1)
    df.columns = cols + ["description of target"]

    return df
