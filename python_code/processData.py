import DataWrangling as dw
import multiprocessing
import pandas as pd
import time


# this data processessing takes a while (like 6 minutes or so) because it's being done separately for each year
# this script sets up some multiprocessing to speed things up
# there are a few things I NEED TO GO BACK AND FIX to make the processes more efficient







# split a list into evenly sized chunks
def chunks(l, n):
    return [l[i:i+n] for i in range(0, len(l), n)]

#wrap all the data processing up into one function
def loopThroughYears(listid, yearlist): #NEED TO make this function take a half-way processed df/obj as an argument so don't need to do all the processing in each loop
     start = time.time()
     

     AMs = pd.read_excel('EcoliBreakPoints.xlsx')
     breakpoints = pd.read_excel('EcoliBreakPoints.xlsx')
     AMList = list(AMs['abbreviation'])


     for year in yearlist:
          NARMSbyYear = dw.DataWrangling(dataset='RetailMeats\Downloaded101023\CVM-2020-NARMS-RetailMeatsData.xlsx')
          NARMSbyYear.selectOrganism(genus='EC', species='coli')
          NARMSbyYear.selectHostSpecies(host_species='Cattle')
          NARMSbyYear.selectYear(year= year)
          NARMSbyYear.selectAM(AMList=AMList, dropMissing= 'col-row', missThresh=1)
          NARMSbyYear.getResistanceStatus(dataframe=NARMSbyYear.MICdf, breakpoints=breakpoints)
          NARMSbyYear.phenotypeIncidence(dataframe=NARMSbyYear.resistancedf, title=year)
     end = time.time()
     print('time=',end-start)

          


#function to split into multiple processes to run simultaneously
def dispatch_jobs(data, job_number): #job number doesn't work out exactly (for example to get 4 jobs job number should be 3)
    total = len(data)
    chunk_size = total / job_number
    slice = chunks(data, int(chunk_size))
    print(slice)
    jobs = []

    for i, s in enumerate(slice):
        print(s)
        j = multiprocessing.Process(target=loopThroughYears, args=(i,s))
        jobs.append(j)
    for j in jobs:
        j.start()
    

     




if __name__ == "__main__":
    
     
     print('started processes!')
     yeardf = pd.read_excel('RetailMeats\Downloaded101023\CVM-2020-NARMS-RetailMeatsData.xlsx', usecols=['Year'])
     print(list(yeardf['Year'].unique()))

     data = list(yeardf['Year'].unique())
     dispatch_jobs(data, 3)
