import pandas as pd
import numpy as np
import time as time
import operator
from collections import Counter
import matplotlib.pyplot as plt
import glob
import re

#Make a class to handle all data wrangling needs INCLUDING:
# 1) [DONE] select which organism you want to look at (for example for this project we'll be looking at E. coli)
# 2) [DONE] select which host species you want to look at (for example for this project we'll only be looking at samples isolated from cattle)
# 3) [DONE] select which antimicrobials you want to look at
# 4) [DONE] map MIC values to breakpoints in order to discretize the data into resitant (1) and susceptible (0)
# 5) [DONE] create incidence matrix for phenotypic resistance
# 6) create incidence matrix for resistance genes
# 7) [DONE] get specific year
# 8) [DONE] calculate percent missing data per year per antimicrobial and put into table
class DataWrangling():
    #construct by loading in the data and making class attributes then pass to other methods
    def __init__(self, dataset):
       
        self.data = pd.read_excel(str(dataset))#, nrows=5000) #LIMIT ROWS JUST FOR FAST TESTINGS !!!REMEMBER TO CHANGE!!!!        
        #drop all no Growth
        self.data = self.data.drop(self.data[self.data['GROWTH'] == 'NO'].index)



        pass

    
    
    
    
    
    # 1) function to select which organism to look at
    def selectOrganism(self, genus, species):
        
        if genus == 'ALL':
            self.data = self.data
        else:
            self.data = self.data.drop(self.data[self.data['GENUS'] != str(genus)].index)
        if species == 'ALL':
            self.data = self.data
        else:
            self.data = self.data.drop(self.data[self.data['SPECIES'] != str(species)].index)

    
    
    
    
    
    
    # 2) function to select which host species to look at
    def selectHostSpecies(self, host_species):
            self.data = self.data.drop(self.data[self.data['HOST_SPECIES'] != str(host_species)].index)

    
     
    def selectYear(self, year):
            self.data = self.data.drop(self.data[self.data['Year'] != year].index)
         

    
    
    
    
    
    
    # 3) function to select which antimicrobials to look at
    def selectAM(self, AMList, dropMissing, missThresh = None):
         self.MICdf = self.data #make new df because going to want original data df when getting genotype info

         #make everything uppercase to make sure there's no mismatch when checking the AMList against the df
         self.MICdf.columns = map(str.upper, self.MICdf.columns)
         self.AMList = list(map(lambda x: x.upper(), AMList))
         self.AMList= list(set(self.AMList)) #make sure there's no duplicates

         self.AMListSign = list(map(lambda x: x + ' SIGN', self.AMList)) # Make an AMList with sign associated with it
         keepCols = self.AMList + self.AMListSign #concat sign list with AMList
         newKeepCols = ['SAMPLE_ID'] + ['YEAR'] + keepCols #add the sample ID and Year to the columns we want to keep
         self.MICdf = self.MICdf.loc[:, newKeepCols] #define the df as only the columns we want to keep
         

         if dropMissing == 'row-col': #drop missing data by rows first and then remaining missing columns
               self.MICdf = self.MICdf.dropna(thresh=self.MICdf.shape[1]-missThresh, axis=0)# drop rows with 5 or more miss values
               self.MICdf = self.MICdf.dropna(how='all', axis=1) #drop any columns with remaining missing values

               #the dropna above dropped some Antimicrobial columns which were tested sporadically (had lots of missing data)
               #the next two lines just update the AMList and AMListSign to reflect those dropped Antimicrobials
               self.AMList = [AM for AM in self.AMList if AM not in list((Counter(self.AMList) - Counter(self.MICdf.columns[2:])).elements())]
               self.AMListSign = [AM for AM in self.AMListSign if AM not in list((Counter(self.AMListSign) - Counter(self.MICdf.columns[2:])).elements())]
         elif dropMissing == 'col-row': #drop missing data by columns and then by remaining missing rows
               self.MICdf = self.MICdf.dropna(how='all', axis=1) #drop any columns with remaining missing values
               self.MICdf = self.MICdf.dropna(thresh=self.MICdf.shape[1]-missThresh, axis=0)# drop rows with 5 or more miss values
               
               self.AMList = [AM for AM in self.AMList if AM not in list((Counter(self.AMList) - Counter(self.MICdf.columns[2:])).elements())]
               self.AMListSign = [AM for AM in self.AMListSign if AM not in list((Counter(self.AMListSign) - Counter(self.MICdf.columns[2:])).elements())]
         elif dropMissing == 'rows':
               self.MICdf = self.MICdf.dropna(thresh=self.MICdf.shape[1]-missThresh, axis=0)# drop rows with 5 or more miss values
               
               self.AMList = [AM for AM in self.AMList if AM not in list((Counter(self.AMList) - Counter(self.MICdf.columns[2:])).elements())]
               self.AMListSign = [AM for AM in self.AMListSign if AM not in list((Counter(self.AMListSign) - Counter(self.MICdf.columns[2:])).elements())]
         elif dropMissing == 'columns':
               self.MICdf = self.MICdf.dropna(how='all', axis=1) #drop any columns with remaining missing values
              
               self.AMList = [AM for AM in self.AMList if AM not in list((Counter(self.AMList) - Counter(self.MICdf.columns[2:])).elements())]
               self.AMListSign = [AM for AM in self.AMListSign if AM not in list((Counter(self.AMListSign) - Counter(self.MICdf.columns[2:])).elements())]
         else: self.MICdf = self.MICdf
              
         

         # keep copy of wide format for other utilities
         self.MICdfWide = self.MICdf

        #next 6 lines are just reshaping the dataframe into the long format (will be helpful later)
         self.MICdf1 = self.MICdf.melt(id_vars='SAMPLE_ID', value_vars= self.AMList, var_name='AM', value_name="MIC")
         self.MICdf2 = self.MICdf.melt(id_vars='SAMPLE_ID', value_vars= self.AMListSign, var_name='AM SIGN', value_name="SIGN")
         self.MICdf3 = self.MICdf.melt(id_vars='YEAR', value_vars= self.AMListSign, var_name='AM SIGN', value_name="SIGN")
         self.MICdf3= self.MICdf3.drop(['SIGN', 'AM SIGN'], axis=1)
         self.MICdf2= self.MICdf2.drop(['SAMPLE_ID', 'AM SIGN'], axis=1)
         self.MICdf = pd.concat([ self.MICdf3, self.MICdf1, self.MICdf2], axis=1)



         return self.MICdf #return the reshaped dataframe with only the wanted columns
    
   
   
   
   
   
    # 4) function to map MIC values to resistance status
    def getResistanceStatus(self, dataframe, breakpoints):
         self.resistancedf = dataframe #dataframe is the MICdf from the selectAM function
        
        #step 1 make min and max columns to describe the MIC in interval form
         self.resistancedf['MIN'] = 'NULL'
         self.resistancedf['MAX'] = 'NULL'
         self.resistancedf.loc[(self.resistancedf['SIGN'] == '='), ['MAX', 'MIN']] = self.resistancedf['MIC']
         self.resistancedf.loc[(self.resistancedf['SIGN'] == '<='), 'MAX'] = self.resistancedf['MIC']
         self.resistancedf.loc[(self.resistancedf['SIGN'] == '<='), 'MIN'] = 0
         self.resistancedf.loc[(self.resistancedf['SIGN'] == '>='), 'MIN'] = self.resistancedf['MIC']
         self.resistancedf.loc[(self.resistancedf['SIGN'] == '>='), 'MAX'] = 999999
         self.resistancedf.loc[(self.resistancedf['SIGN'] == '>'), 'MAX'] = 999999
         self.resistancedf.loc[(self.resistancedf['SIGN'] == '<'), 'MIN'] = 0
         self.resistancedf.loc[(self.resistancedf['SIGN'] == '>'), 'MIN'] = self.resistancedf['MIC'] + 0.1


         

         # now need to check min and max values to see if susceptible
         self.resistancedf.loc[:, 'isResistant'] = 0 #create column for resistances status
         self.resistancedf.loc[:,'susceptible value'] = '' #create column to hold the NARMS sucseptibility breakpoint
         
         
         

         #take the susceptibility breakpoint from the breakpoints data and add it to the resistance df
         #can probably just compare two separate dataframes but I found this easier
         for index, row in self.resistancedf.iterrows():
                for i, point in breakpoints.iterrows():
                    if row['AM'] == point['abbreviation']:
                         self.resistancedf.loc[index, 'susceptible value'] = point['susceptible value']
         

         #STR needs some special attention because breakpoints change based on the year tested 
         for index, row in self.resistancedf.iterrows():
              if (row['AM'] == 'STR') and (row['YEAR'] < 2014):
                   self.resistancedf.loc[index, 'susceptible value'] = 32
              elif (row['AM'] == 'STR') and (row['YEAR'] > 2014):
                   self.resistancedf.loc[index, 'susceptible value'] = 16
              else: self.resistancedf.loc[index, 'susceptible value'] = row['susceptible value']
         #^(would be nice to solve this problem in a more general purpose way though) 

         # go through and check whether the max of the MIC interval is greater than the sucsecptibility breakpoint
         # if true then that means the isolate is not susceptible (code 1 for not susceptible and 0 for susceptible)
         for index, row in self.resistancedf.iterrows():
              if row['MAX'] > row['susceptible value']:
                   self.resistancedf.loc[index, 'isResistant'] = 1
        
        
         for index, row in self.resistancedf.iterrows():
              if (row['MAX'] > row['susceptible value']) and (row['MIN'] < row['susceptible value']):
                   print("UNINTERPRETABLE")
                   print(" ")
                   print(row)
                   print(" ")



         for index, row in self.resistancedf.iterrows():
              if (row['MAX'] > row['susceptible value']) and (row['MIN'] < row['susceptible value']):
                   self.resistancedf = self.resistancedf.drop(index = index, axis=0)


         for index, row in self.resistancedf.iterrows():
              if (row['MAX'] > row['susceptible value']) and (row['MIN'] < row['susceptible value']):
                   print("UNINTERPRETABLE")
                   print(" ")
                   print(row)
                   print(" ")

         #self.resistancedf['isResistant'] = pd.cut(row['MAX'], bins=[-np.inf, self.resistancedf['susceptible value'], np.inf], labels=[0,1])
         #^I can't figure out why this won't work, but it would be so much faster than looping through the df... so try to figure it out later    
         
         #save to csv so don't have to keep running this (though it doesn't really take that long to run so not totally necessary)
         #self.resistancedf.to_csv('EcoliResistanceStatus.csv', index=False)


         #print(self.resistancedf)

         return self.resistancedf
    


    # 5) make incidence matrix for phenotypic resistance
    def phenotypeIncidence(self, dataframe, title):
         self.incidenceMatrix = pd.DataFrame(dataframe) #takes a dataframe as an argument so that I can do either option: 1) load in df from csv 2) pass the returned df from function number 4
         SAMPLE_IDs_unique = np.array(self.incidenceMatrix['SAMPLE_ID'].unique())


         self.incidenceMatrix=self.incidenceMatrix.pivot(columns='AM', values='isResistant')
         v = self.incidenceMatrix.values
         i = np.arange(v.shape[1])
         a = np.isnan(v).argsort(0, kind='mergesort')
         v[:] = v[a, i]
         self.incidenceMatrix = self.incidenceMatrix.dropna()
         #self.incidenceMatrix['SAMPLE_ID'] = SAMPLE_IDs_unique

         self.incidenceMatrix.to_csv('incidenceMatrices/' + str(title) + '_incidenceMatrix.csv', index=False)
         #print(len(SAMPLE_IDs_unique))

     













    def getPercentMissing(self):
         years = list(self.MICdfWide['YEAR'].unique())

         for year in years:
              self.selectYear(year)
              percent_missing = self.selectedYeardf.isnull().sum() * 100 / len(self.selectedYeardf)
              missing_value_df = pd.DataFrame({'column_name': self.selectedYeardf.columns,
                                               'percent_missing': percent_missing})
              
              #print(missing_value_df)

              fig = plt.figure()
              ax = fig.add_subplot()

              ax.table(cellText = missing_value_df.values,
              #rowLabels = missing_value_df.index,
              colLabels = missing_value_df.columns,
              loc = "center")
              ax.set_title("Percent Missing Data Per Antimicrobial\nNARMS Retail Meats " + "(" + str(year) + ")")

              ax.axis("off")

              plt.savefig('missingData/RetailMeats/' + str(year) + 'missingdata.png')
              plt.clf()
         

         

       




AMs = pd.read_excel('EcoliBreakPoints.xlsx')
breakpoints = pd.read_excel('EcoliBreakPoints.xlsx')
AMList = list(AMs['abbreviation'])


NARMSbyYear = DataWrangling(dataset='RetailMeats\Downloaded101023\CVM-2020-NARMS-RetailMeatsData.xlsx')
NARMSbyYear.selectOrganism(genus='EC', species='coli')
NARMSbyYear.selectHostSpecies(host_species='Cattle')
NARMSbyYear.selectYear(year= 2003)
NARMSbyYear.selectAM(AMList=AMList, dropMissing= 'col-row', missThresh=1)
NARMSbyYear.getResistanceStatus(dataframe=NARMSbyYear.MICdf, breakpoints=breakpoints)
#df = NARMSbyYear.resistancedf
#df = df[df["isResistant"] == 1]
#print(df)
#NARMSbyYear.phenotypeIncidence(dataframe=NARMSbyYear.resistancedf, title=2005)

print('ran!')



def get_prevalence_per_resistance():
     path = "incidenceMatrices/*.csv"
     df_list = []
     for fname in glob.glob(path):
          year = re.sub(r"\D", "", fname)
          df = pd.read_csv(fname)
          df = df.drop("SAMPLE_ID", axis = 1)
          for col in df.columns:
               #df[col] = np.round((df[col].sum()/len(df[col])) * 100, 1)
               df[col] = np.round(df[col].sum(), 1)

               #df[col] = df[col].apply( lambda x : str(x) + '%')


          df["year"] = year
          df = df[['year'] + [ col for col in df.columns if col != 'year' ]]
          df = df.loc[[0]]
          df_list.append(df)
     full_df = pd.concat(df_list, axis = 0, ignore_index=True)
     print(full_df)
     full_df.to_csv("missingData/RetailMeats/individual_AM_resistance_counts.csv")




#get_prevalence_per_resistance()

    

     
        

        
        








