
Genome-Scale Multilabel Prediction for Tissue-specific Annotation of Long
Non-coding RNAs
Submitted to Great Lakes Bioinformatics Conference 2019
Almulhim, Aljohara
------------------------------------------------------------------------------------------------------
PreProcessing: Filtring and Cascading processes
------------------------------------------------------------------------------------------------------
To Run: python pro_data.py data/ENSG_HPID_labels.txt data/ENSG_DOID_labels.txt data/Liver_HPIDs.txt data/Liver_DOIDs.txt data/lncRNA_ENSG.txt data/LIHC_matrix.txt data/HPID_height.txt data/DOID_height.txt data/HPID_child_parent.txt data/DOID_child_parent.txt

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
To run: python NN_model.py LIHC_train_data Tag_matrix Height_file LIHC_test_data Liver_tag_lst lncRNA_lst learning_rate epoch_num batch_size hidden_units_num results/log.txt
LIHC_test_data: Our training data in which rows are number of patients samples and the columns are ENSG_gene IDs
Tag_matrix: pass the tag_matrix or the TE_tag_matrix where the rows are ENSG_gene IDs and the columns are the tags
Height file: list of tags and their depth(Height) based on the tag tree
LIHC_test_data: Our testing data where we predict for each lncRNAs number of tags. It has number of patients samples as rows and lncRNA_IDs as columns
Liver_tag_lst= Liver tags list
lncRNA_lst = lncRNA list
Learning rate: 0.0001-0.01
Hidden unit numbers: 100-400
Batch size: 50-150
Number of epochs: 100-200

------------------------------------------------------------------------------------------------------
Predictions: Thresholds and consistency 
------------------------------------------------------------------------------------------------------
