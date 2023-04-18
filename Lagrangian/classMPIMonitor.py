import pandas as pd
import sys

class MPIMonitor:

    def __init__(self):

        self.m_DF = pd.DataFrame(columns=['Case', 'Integrating from', 'to', 'simTime'])

    def setDF(self, df):
        self.m_DF = df
        # print(self.m_DF)
        # sys.exit('\n\n\ncheck')

    def updateDF(self, row):
        ind = self.m_DF.index[self.m_DF['Case']==row[0]]
        row = pd.DataFrame(row).T#, columns=['Case', 'Integrating from', 'to', 'simTime'])
        self.m_DF.loc[ind] = row

    def outputDF(self, name):
        self.m_DF.to_csv(name, index=False)
        # sys.exit('check')
