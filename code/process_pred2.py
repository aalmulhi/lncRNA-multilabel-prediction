import sys;
import pandas as pd
import numpy as np
import pickle
import matplotlib

#reading prediction probabilities 
def read_pred(fn):
    df= pd.read_csv(fn,sep=",",index_col='Unnamed: 0')
    return df

# Load anst-desct information for liver
def load_dict(fn):
    with open(fn, 'rb') as handle:
        d= pickle.load(handle)
    return d

#check cascading
def check_cascading(df,d):
    my_lst=[]
    for i, row in df.iterrows():
        for j, column in row.iteritems():
            for k,v in d.items():
                if column==1 and j==k:
                    flag= True
                    for iv in v:
                        if df.loc[i,iv] != 1:
                            flag=False
                            break
                    if(flag==True):
                        my_lst.append(i)
    l=set(my_lst)
    return list(l)

def main():
        if len(sys.argv) < 4:
                print "Enter files as following:"
                print "pro_pred.py prediction_prob_file Liver_anst_desc_dict_file threshold"
        else:
		fn1 = open(sys.argv[1]);
		fn2 = sys.argv[2]
		th= float(sys.argv[3])

                name_p = sys.argv[1].strip().rsplit(".",1)[0].strip().rsplit("/",1)[-1];
		pred_prob=read_pred(fn1)
		liver_dict=load_dict(fn2)
		
                liver_lncRna_lst=list(pred_prob.index.values)
		liver_tag_lst=list(pred_prob.columns.values)
                n= len(liver_lncRna_lst)
                d= len(liver_tag_lst)

                #setting predictions with thresholds
                pred_prob_mtx=pred_prob.values
                liver_preds_test= np.zeros(shape=(n, d))
                liver_preds_test[pred_prob_mtx>= th] = 1
                liver_preds_test[pred_prob_mtx< th] = 0
                liver_preds_test_df= pd.DataFrame(data=liver_preds_test, index=liver_lncRna_lst, columns=liver_tag_lst)
                
                #if label predicted as 1 make sure thie prediction is constistant with ansestors
                lncRnaLst_lvr=check_cascading(liver_preds_test_df,liver_dict)
                liver_preds_test_df=liver_preds_test_df.loc[lncRnaLst_lvr]

                #drop tag columns with no predictions
                liver_preds_test_df=liver_preds_test_df.loc[:, (liver_preds_test_df != 0).any(axis=0)]

                #save prediction with threshold=th
                liver_preds_test_df.to_csv("results/"+name_p+"_"+str(th)+".txt",sep=",")

                

                             
main();
