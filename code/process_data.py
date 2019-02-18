import sys;
import pandas as pd
import numpy as np
import networkx as nx
# Reading ENSG with the associated tags saparated by tab
# ENSG1;tag1;tag2;tag3
# ENSG2,tag2;tag5;tag7 ...
def read_ENSG_file(fn,fn_tag):
    df=pd.read_csv(fn,sep='\t', header=None)
    ID_=df.values
    ID_dict={}
    for i in ID_:
        if i[0] in ID_dict:
            ID_dict[i[0]].append(i[1].split(";"))
        else:
            ID_dict[i[0]] = [i[1].split(";")]
    ID_lst=[]
    ID_lst= ID_dict.keys()
    return ID_lst,ID_dict

# Reading tags with thier definitions
# tag_ID \t tag_name
def read_tag_file(fn,tn):
    df=pd.read_csv(fn,sep='\t', header=None)
    mtx=df.values
    n,d=mtx.shape
    mtx_lst=mtx[1:,0]
    lst= list(set(mtx_lst))
    return lst

# Reading lncRNA file
def read_lncRNA(fn,ENSG_lst):
        lncRNA_lst =fn.read().splitlines()
        lncRNA_matched_lst=list(set(lncRNA_lst) & set(ENSG_lst))
        ENSG_lst_Excld=[item for item in ENSG_lst if item not in lncRNA_lst]
        return lncRNA_lst,lncRNA_matched_lst,ENSG_lst_Excld

# Reading LIHC file
def disease_matrix_processing(fn,lncRNA_IDs_Ex,lncRNA_lst):
    #LIHC
    df= pd.read_csv(fn,sep='\t')
    df1=df.fillna(df.mean())
    
    #training LIHC excluding the lncRNAs
    df2=df1.loc[:,lncRNA_IDs_Ex]

    #testing LIHC including only lncRNAs
    df1_col=list(df1.columns.values)
    match_lncRNAs=list(set(df1_col) & set(lncRNA_lst))
    df3=df1.loc[:,match_lncRNAs]

    return df2,df3


#Create tag matrix using ENSG gene list and liver list
def create_tag_df(tag_lst,gene_lst):
    tag_df=pd.DataFrame(0.0,columns=tag_lst,index=gene_lst)
    return tag_df

#Fill tag matrix in which ENSG gene and tag has a 1 if they are assosiated, otherwise 0
def fill_tags(tag_df,tag_dict):
    for k,v in tag_dict.items():
        for v1 in v:
            for v2 in v1:
                try:
                    tag_df.loc[[k],[v2]]=1.0
                except:
                    pass
    return tag_df

#Filter_Process: Drop ENSG genes with no assosiatd tag amd vice versa
def filter_process(LIHC_train_df,tag_df):
    tag_df1=tag_df[tag_df.values.sum(axis=1) != 0]
    gene_index_lst=list(tag_df1.index.values)
    tag_df1=tag_df1.loc[:, (tag_df1!= 0).any(axis=0)]
    LIHC_train_df=LIHC_train_df[gene_index_lst]
    tag_df1=tag_df1.sort_index(axis=1)
    return LIHC_train_df,tag_df1

#Storing training and testing data
def store_df(df,fn):
    df.to_csv(fn,sep=",")

#store lists
def store_lst(fn,lst):
    np.savetxt(fn, lst, delimiter=",", fmt='%s')
    
#Reading depth information 
def read_depth(fn):
    df=pd.read_csv(fn,sep='\t', header=None)
    mtx=df.values
    mtx=mtx[1::]
    lst=[]
    for i in mtx:
        i=i.tolist()
        lst.append(i)
    for l in lst:
        l[1]=int(l[1])
    return lst

#Storing depth for the tags we will use in training
def create_depth_df(Hg_lst,tag_lst):
    Hg_df=pd.DataFrame(Hg_lst,columns=['Label','Height'])
    Hg_df['Label'] =Hg_df['Label'].str.replace('_',':')
    Hg_df=Hg_df.loc[Hg_df['Label'].isin(tag_lst)]
    Hg_df.set_index('Label',inplace=True)
    Hg_df=Hg_df.sort_index()
    return Hg_df


#reading parent-child information
def read_tree_ids(fn):
    df=pd.read_csv(fn,sep='\t', header=None)
    df[0]=df[0].str.replace('_',':')
    df[1]=df[1].str.replace('_',':')
    tree=df.values
    tree_dict={}
    for i in tree:
        if i[0] in tree_dict:
            tree_dict[i[0]].append(i[1].split("|"))
        else:
            try:
                tree_dict[i[0]] = [i[1].split("|")]
            except:
                pass 
    return tree_dict

#filtring tree-dict to have only tags show in my data
def filter_dict_keys(merged_tree_dict,tag_lst):
    f_dict_k=dict((k, merged_tree_dict[k]) for k in tag_lst if k in merged_tree_dict)

    f_dict_v={}
    for k,v in f_dict_k.items():
        for iv in v:
            if iv[0] in tag_lst:
                if k in f_dict_v:
                    f_dict_v[k].append(iv)
                else:
                    f_dict_v[k]=iv               
    return f_dict_v

#Using networkx for ansesstors and descendant
def calc_anst_desc_dict(f_dict):
    G= nx.DiGraph()
    for k, v in f_dict.items():
        ancestor= v[0]
        descendant = k
        G.add_edge(ancestor, descendant)
    anst_desc_dict={}
    for k, v in f_dict.items():
        anst=(nx.ancestors(G, k))
        anst_desc_dict[k]=list(anst)
    return anst_desc_dict

#Cascading tags to increase the relevance of tags and to get the TE tag matrix
def cascading_process(Liver_tag_df,anst_desc_dict):
    for i, row in Liver_tag_df.iterrows():
        for j, column in row.iteritems():
            for k,v in anst_desc_dict.items():
                if column==1 and j==k:
                    for iv in v:
                        try:
                            Liver_tag_df.loc[i,iv] = 1.0
                        except:
                            pass
    return Liver_tag_df

#store anst_desc dict
def store_dict(fn,d):
    '''
    with open(fn, "w") as file:
        for k, v in d.items():
            dictionary_content = k + "\t" + str(v) + "\n"
            file.write(dictionary_content)
    print "done"
    '''
    with open(fn, 'wb') as handle:
        pickle.dump(d, handle, protocol=pickle.HIGHEST_PROTOCOL)
def main():
        if len(sys.argv) < 11:
                print "Enter files as following:"
                print "row_data.py ENSG_HPID_file ENSG_DOID_file Liver_HPID_file Liver_DOID_file lncRNA_ENSG_file LIHC_matrix_file HPID_Height_file DOID_Height_file HPID_tree DOID_tree"
        else:
		fn1 = open(sys.argv[1]);
		fn2 = open(sys.argv[2]);
		fn3 = open(sys.argv[3]);
		fn4 = open(sys.argv[4]);
		fn5 = open(sys.argv[5]);
		fn6 = open(sys.argv[6]);
		fn7 = open(sys.argv[7]);
		fn8 = open(sys.argv[8]);
		fn9 = open(sys.argv[9]);
		fn10 = open(sys.argv[10]);
		
		print "=============================================================="
                print "Preprocessing Data...."
                print "=============================================================="
                ENSG_HPID_lst,ENSG_HPID_dict= read_ENSG_file(fn1,"ENSG_HPID");
                ENSG_DOID_lst,ENSG_DOID_dict= read_ENSG_file(fn2,"ENSG_DOID");

                ENSG_IDs= list(set(ENSG_DOID_lst + ENSG_HPID_lst))
        
                Liver_HPID=read_tag_file(fn3,"HPID")    
                Liver_DOID=read_tag_file(fn4,"DOID")

                Liver_tags= list(set(Liver_HPID + Liver_DOID))
                Liver_tags= [e.replace('_', ':') for e in Liver_tags]

                lncRNA_IDs,lncRNA_IDs_matched,lncRNA_IDs_Ex= read_lncRNA(fn5,ENSG_IDs)

                LIHC_train,LIHC_test=disease_matrix_processing(fn6,lncRNA_IDs_Ex,lncRNA_IDs)
    
                Liver_tag_df= create_tag_df(Liver_tags,lncRNA_IDs_Ex)
                Liver_tag_df= fill_tags(Liver_tag_df,ENSG_HPID_dict)
                Liver_tag_df= fill_tags(Liver_tag_df,ENSG_DOID_dict)
                
                LIHC_train,Liver_tag_df= filter_process(LIHC_train,Liver_tag_df)
                tag_lst= list(Liver_tag_df.columns.values)
                gene_lst= list(Liver_tag_df.index.values)
                lncRna_test_lst= (list(LIHC_test.columns.values))
                print "Step1: filtring process is done.."

                DOID_Hg_lst= read_depth(fn7)
                HPID_Hg_lst= read_depth(fn8)
                Hg_lst= DOID_Hg_lst + HPID_Hg_lst
                Hg_df= create_depth_df(Hg_lst,tag_lst)

                HPID_tree_dict=read_tree_ids(fn9)
                DOID_tree_dict=read_tree_ids(fn10)
                merged_tree_dict= DOID_tree_dict.copy()
                merged_tree_dict.update(HPID_tree_dict)
                f_tree_dict= filter_dict_keys(merged_tree_dict,tag_lst)
                anst_desc_dict= calc_anst_desc_dict(f_tree_dict)
                Liver_TE_tag_df= cascading_process(Liver_tag_df,anst_desc_dict)

                print "Step2: cascading process is done.."

                
                store_df(LIHC_train,'pro_data/LIHC_train.csv')
                store_df(LIHC_test,'pro_data/LIHC_test.csv')
                store_df(Liver_tag_df,'pro_data/Liver_tag.csv')
                store_df(Hg_df,'pro_data/Liver_Depth.csv')
                store_df(Liver_TE_tag_df,'pro_data/Liver_TE_tag.csv')
                store_lst("pro_data/tag_lst.csv", tag_lst)
                store_lst("pro_data/lncRNA_lst.csv", lncRna_test_lst)
                store_dict("pro_data/Liver_anst_desc_dict.pickle",anst_desc_dict)
		
                print "Saving files is done .."
                            
main();
