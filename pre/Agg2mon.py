def agg_mon(yr,mo):
    ''' aggregate all NOAA20 VIIRS 375m fire location data and create a monthly file similar to VNP14IMGML
    '''
    from datetime import date
    import calendar
    import os
    from glob import glob
    import pandas as pd
    
    # set parameters
    strdir = 'D:/shared_data/VIIRS_375m_NOAA20/'
    dirmon = strdir + str(yr)+'/'+str(mo).zfill(2)+'/'
    nd = calendar.monthrange(yr,mo)[1]   # total number of days within a month
    
    # loop over each day in the month
    monst = True
    for d in range(nd):
        print(date(yr,mo,d+1))
        
        # find all text files in a day
        dird = dirmon + str(d+1).zfill(2)
        fnms = glob(dird+'/*.txt')

        # loop over all files in a day and make aggregated daily dataframe dfday
        dayst = True
        for fnm in fnms:

            # read data for each file
            df = pd.read_csv(fnm)
            
            # conform the name of the last column 
            df.columns = [*df.columns[:-1], 'persist_anomaly']

            # aggregate dataframe from all files
            if dayst:
                dfday = df
                dayst = False
            else:
                dfday = dfday.append(df,ignore_index=True)    

        # aggregate over all days in a month and make the monthly dataframe dfmon
        if len(fnms) > 0:
            if monst:
                dfmon = dfday
                monst = False
            else:
                dfmon = dfmon.append(dfday,ignore_index=True)
            print(len(dfmon))

    # save the aggregated data to the yearly directory
    fnmout = strdir + 'monthly/' + str(yr)+'/VJ114IMGML_'+str(yr)+str(mo).zfill(2)+'.txt'
    dfmon.to_csv(fnmout)
    
if __name__ == "__main__":
    for mo in range(1,13):
        print(mo)
        try:
            agg_mon(2018,mo)
        except:
            continue
