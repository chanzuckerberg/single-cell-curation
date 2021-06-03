import pandas as pd
import sys
import os

def main():
    
    sample_info_file = './data/original/sample_info.tsv'
    out_file = './data/original/counts.tsv.gz'
    
    sample_info = pd.read_csv(sample_info_file, sep='\t')
    
    for i in range(sample_info.shape[0]):
    #for i in range(2):
        print(i)
        current = pd.read_csv(os.path.join(os.path.dirname(sample_info_file), sample_info.loc[i, 'file']), sep='\t', index_col=[0,1])
        to_keep = ['Unnamed' not in k for k in current.columns]
        current = current.loc[:, to_keep]
        current.columns = current.columns + '_' + sample_info.loc[i, 'donor'] + '_' + sample_info.loc[i, 'tissue'] + '_' + sample_info.loc[i, 'state']
        
        if i == 0:
            table = current
        else:
            table = table.join(current)
    
    table.reset_index(inplace=True)
    table.drop('Accession', axis=1, inplace=True)
    table.to_csv(out_file, sep='\t', compression='gzip', index=False)
    
    
if __name__=="__main__":
    main()
        
