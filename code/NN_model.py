from numpy import genfromtxt
from skmultilearn.problem_transform import LabelPowerset
import pandas as pd
import numpy as np
import numpy.ma as ma
import keras
import keras.backend as K
from keras.models import Sequential
from keras.layers import Dense
import numpy as np
from sklearn.metrics import f1_score,recall_score,precision_score
import timeit
from sklearn.preprocessing import normalize
seed=7
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import train_test_split
np.set_printoptions(suppress=True)
import statistics as s
import warnings
import sys
warnings.simplefilter('ignore')

def read_data(fn):
    df=pd.read_csv(fn,delimiter=",",index_col=0)
    df1 = df.convert_objects(convert_numeric=True)
    mtx=df1.as_matrix()
    print mtx.shape
    return mtx

def main():
        if len(sys.argv) < 12:
            print "Enter files as following:"
            print "NN_model.py LIHC_train_data Tag_matrix Height_file LIHC_test_data liver_tag_lst lncRNA_lst learning_rate epoch_num batch_size hidden_units_num results/log.txt"  
        else:
            fn1= open(sys.argv[1]);
            fn2= open(sys.argv[2]);
            fn3= open(sys.argv[3]);
            fn4= open(sys.argv[4]);
            fn5= open(sys.argv[5]);
            fn6= open(sys.argv[6]);
            lrt = float(sys.argv[7]);
            epc= int(sys.argv[8]);
            bs= int(sys.argv[9]);
            hs= int(sys.argv[10]);
            lg= str(sys.argv[11]);

            name_t_p = sys.argv[2].strip().rsplit(".",1)[0].strip().rsplit("/",1)[-1];

            print "---------------------------------------------------------------"
            print "My training data:"
            print "LIHC data:"
            X=read_data(fn1)
            print "Tag Matrix:"
            y=read_data(fn2)
            print "My Height:"
            Hg=read_data(fn3)
            print "Predictions data:"
            test_lncRNA_data= read_data(fn4)
            liver_tag_lst= fn5.read().splitlines()
            print "Number of tags to predict",len(liver_tag_lst)
            lncRna_lst = fn6.read().splitlines()
            print "Number of lncRNAs to predict",len(lncRna_lst)
            print "Parameters:"
            print "Learning Rate:",lrt
            print "Epoch:",epc
            print "Batch Size:",bs
            print "# Hidden nodes:",hs
            print "---------------------------------------------------------------"

            #fixing training
            X1=np.where(np.isnan(X), ma.array(X, mask=np.isnan(X)).mean(axis=0), X) 
            normX = normalize(X1, axis=0)
            normX = normX.T
            rc_lst=[]
            rc_h_lst=[]
            #splitting to train and test
            lp = LabelPowerset()
            # remember to set n_splits and shuffle!
            kf = StratifiedKFold(n_splits=5, random_state=111, shuffle=True)
            k=1
            for train_index, test_index in kf.split(normX, lp.transform(y)):
                #splitting
                print "processing.."
                X_train = normX[train_index,:]
                y_train = y[train_index,:]
                X_test = normX[test_index,:]
                y_test = y[test_index,:]
                x_size =X_train.shape[1]  
                y_size = y_train.shape[1] 
        
                #start
                start = timeit.default_timer()
                model = Sequential()
    
                #layers
                model.add(Dense(hs, activation='relu', input_dim=x_size))
                model.add(Dense(x_size, activation='relu'))
                model.add(Dense(y_size, activation='sigmoid'))

                #model oprimizer
                ada=keras.optimizers.Adagrad(lr=lrt, epsilon=None)

                #model
                model.compile(loss='binary_crossentropy',optimizer=ada,metrics=['categorical_accuracy'])
                history_callback=model.fit(X_train,y_train, batch_size=bs, epochs=epc,verbose=0)
                loss_history = history_callback.history["loss"]
                numpy_loss_history = np.array(loss_history)
                np.savetxt(lg, numpy_loss_history, delimiter=",")
         
                #time is up!
                stop = timeit.default_timer()
                c_time=stop - start
    
                #predict test data
                preds = model.predict(X_test)
                preds[preds>= 0.50] = 1
                preds[preds< 0.50] = 0
        
                #Recall score
                rc=recall_score(y_test, preds, average="macro")
                rc_lst.append(rc)

            
                #Depth Adjusted Recall score
                numerator = 0.0;
                for i in range(len(preds)):
                    for j in range(len(preds[0])):
                        if preds[i][j]==1 and y_test[i][j] == 1: 
                            numerator += Hg[j]
                dinominator = sum(np.dot(y_test,Hg))
                recall_w_H=numerator /dinominator
                rc_h_lst.append(recall_w_H)
                print "Fold:",k
                print "Recall=",rc
                print "Depth Adjusted Recall=",recall_w_H[0]
                k+=1

        
            #final score for all folds
            print "============================================="
            print "Final Scores:"
            print "Recall=",s.mean(rc_lst)
            rc_h_lst1=np.concatenate(rc_h_lst, axis=0 )
            print "Depth Adjusted Recall=", s.mean(rc_h_lst1)
            print "============================================="

            #test_data lncRNA prediction
            normXtest = normalize(test_lncRNA_data, axis=0)
            preds_test_prob = model.predict(normXtest.T)
            print preds_test_prob.shape

            prediction_df= pd.DataFrame(data=preds_test_prob, index=lncRna_lst, columns=liver_tag_lst)
            print prediction_df.info()
            prediction_df.to_csv("results/prediction_"+name_t_p+".txt",sep=",")
            
main();
