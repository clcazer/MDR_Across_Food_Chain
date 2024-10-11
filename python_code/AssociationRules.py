import pandas as pd
import numpy as np
from mlxtend.frequent_patterns import apriori, association_rules, fpgrowth
from mlxtend.preprocessing import TransactionEncoder
import networkx as nx
import matplotlib.pyplot as plt


#Make class to handle all association mining activities INCLUDING:
# 1) [DONE] run apriori (supply a minsupp)
# 2) [DONE] run fpfgrowth (supply minsupp)
# 3) cut down rules based on several measures
# ^----> need function to calculate measure that are not already built-in?
# ^----> need function to find cutoff values?
# ^----> need function to find out which are the best measures to look at?

class RuleMining():
    def __init__(self):
        pass

    # 1) apriori
    # the the functions from mlxtend already allow freq item selection based on minsupp
    #and then rule selection based on one metric (selectionMetric and threshold)
    #an additional function is needed to make further selections based on multiple metrics
    def runApriori(self, incidenceMatrix, minsupp, maxItemSetLen, selectionMetric, threshold):
        df = pd.read_csv(incidenceMatrix) #read in the incidence matrix (can be produced with the functions in the DataWrangling class)
        df= df.drop('SAMPLE_ID', axis=1) #drop the sample ID's for formatting purposes

        # provides the option to set minsupp to one occurence
        if minsupp == 'one occurence':
            minsupp = 1/len(df)

        #provides the option to set the rule selection threshold to one occurence (if metric is support)
        #might want to do this if you want to save rule selection for the more extensive selection criteria
        if threshold == 'one occurence':
            threshold = 1/len(df)

        #generate frequent itemsets
        frequentItemsdf = apriori(df = df,
                        min_support = minsupp,
                        use_colnames = True,
                        max_len= maxItemSetLen,
                        )
        #mine association rules
        self.df_apriori_ar = association_rules(df= frequentItemsdf,
                                metric = selectionMetric,
                                min_threshold = threshold
                                )
        self.df_apriori_ar = self.df_apriori_ar[self.df_apriori_ar['consequents'].str.len() < 2]
        # return the rules
        return self.df_apriori_ar

    def runFPGrowth(self, incidenceMatrix, minsupp, maxItemSetLen, selectionMetric, threshold):
        df = pd.read_csv(incidenceMatrix) #read in the incidence matrix (can be produced with the functions in the DataWrangling class)
        df= df.drop('SAMPLE_ID', axis=1) #drop the sample ID's for formatting purposes
        
        
        # provides the option to set minsupp to one occurence
        if minsupp == 'one occurence':
            minsupp = 1/len(df)
        
        #provides the option to set the rule selection threshold to one occurence (if metric is support)
        #might want to do this if you want to save rule selection for the more extensive selection criteria
        if threshold == 'one occurence':
            threshold = 1/len(df)
        
        
        #generate frequent itemsets
        frequentItemsdf = fpgrowth(df = df,
                        min_support = minsupp,
                        use_colnames = True,
                        max_len= maxItemSetLen,
                        )
        
        #mine association rules
        self.df_fpg_ar = association_rules(df= frequentItemsdf,
                                metric = selectionMetric,
                                min_threshold = threshold
                                )

        # return the rules
        return self.df_fpg_ar
    


    def visMetricDistibutions(self, data):
        metricList = []
        data = pd.DataFrame(data)
        print(len(data))
        
        data.replace([np.inf, -np.inf], np.nan, inplace=True)
        #data = data.dropna()



        for col in data.columns[2:]:
            metricList.append(col)


        for metric in metricList:
            x = data[metric]
            print(x)
            plt.hist(x, bins=50)
            plt.title(str(metric) + ' distribution')
            plt.xlabel(str(metric) + ' value')
            plt.ylabel('Number of Rules')

            plt.savefig(str(metric)+'Dist.png')
            plt.clf()

        plt.hist2d(data['support'], data['confidence'], bins=[5,20])
        plt.colorbar()
        plt.title('support by confidence')
        plt.xlabel('support value')
        plt.ylabel('confidence value')
        plt.savefig('2dsupp+conf.png')



testRules = RuleMining()
testRules.runApriori( incidenceMatrix='incidenceMatrix.csv',
                     minsupp= 'one occurence',
                     maxItemSetLen= None,
                     selectionMetric= 'support',
                     threshold='one occurence')

testRules.visMetricDistibutions(data= testRules.df_apriori_ar)

df = pd.DataFrame(testRules.df_apriori_ar)
#df= df['consequents']
#df = df[df['consequents'].str.len() < 2]
#print(df)

