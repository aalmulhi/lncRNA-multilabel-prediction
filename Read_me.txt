
Genome-Scale Multilabel Prediction for Tissue-specific Annotation of Long
Non-coding RNAs
Submitted to Great Lakes Bioinformatics Conference 2019
Almulhim, Aljohara
------------------------------------------------------------------------------------------------------
PreProcessing: Filtring and Cascading processes
------------------------------------------------------------------------------------------------------
To Run: python pro_data.py ENSG_HPID_labels ENSG_DOID_labels Liver_HPIDs Liver_DOIDs lncRNA_ENSG LIHC_matrix HPID_height DOID_height HPID_child_parent DOID_child_parent

ENSG_HPID_labels.txt: The labeled ENSG_IDs with HPID_tags from literature
File format: ENSG_ID1 HPID1;HPID2;...
	     ENSG_ID2 HPID3;HPID2;...
ENSG_DOID_labels.txt: The labeled ENSG_IDs with DOID_tags from literature
File format: ENSG_ID1 DOID1;DOID2;...
	     ENSG_ID2 DOID3;DOID2;...
Liver_HPIDs.txt: List of HPID tags for liver tissue 
File format: HPID \t definition 
Liver_DOIDs.txt: List of DOID tags for liver tissue 
File format: DOID \t definition 
lncRNA_ENSG.txt: list of lncRNA_ENSG_IDs
LIHC_matrix.txt: matrix of liver patients samples with ENSG_IDs separated by \t
HPID_height.txt: list of HPIDs with depth information
File format: HPID1 \t D1
	     HPID2 \t D2
DOID_height.txt: list of DOIDs with depth information
File format: DOID1 \t D1
	     DOID2 \t D2
HPID_child_parent.txt: list of descendant-ancestor for HPID tag
File format: HPID_child1 \t HPID_parent1|HPID_parent2|HPID_parent3 ..
	     HPID_child1 \t HPID_parent5|HPID_parent7|HPID_parent4 ..    
DOID_child_parent.txt: list of descendant-ancestor for DOID tag
File format: DOID_child1 \t DOID_parent1|DOID_parent2|DOID_parent3 ..
	     DOID_child1 \t DOID_parent5|DOID_parent7|DOID_parent4 ..  

------------------------------------------------------------------------------------------------------
Model: Neural Network
------------------------------------------------------------------------------------------------------
To run: Python NN_model.py LIHC_train Tag_matrix Tag_depth LncRNA_test tag_list lncRNA_lst Learning_rate Number_of_epochs Batch_size Hidden_unit_number log_file
LIHC_train: Our training data in which rows are number of patients samples and the columns are ENSG_gene IDs
Tag_matrix: pass the tag_matrix or the TE_tag_matrix where the rows are ENSG_gene IDs and the columns are the tags
Tag_depth: list of tags and their depth(Height) based on the tag tree
LncRNA_test: Our testing data where we predict for each lncRNAs number of tags. It has number of patients samples as rows and lncRNA_IDs as columns
Liver_tag_lst= Liver tags list
lncRNA_lst = lncRNA list
Learning_rate: 0.0001-0.01
Hidden_unit_numbers: 100-400
Batch_size: 50-150
Number_of_epochs: 100-2500
log_file: Log of loss function values
------------------------------------------------------------------------------------------------------
Predictions: Thresholds and consistency 
------------------------------------------------------------------------------------------------------
python process_pred.py prediction_prob liver_dict p
prediction_prob: the predictions probabilities
liver_dict: ancestor-descendant liver dictionary
p: Threshold 
