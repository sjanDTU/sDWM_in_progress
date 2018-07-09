import WindFarm as wf
import WindTurbine as wt

def getInputs():
    inputs = dict()
    ##########################################################################
    inputs = {'WD': 222,
              'WS': 9,
              'TI': 0.056,
              'WTcoord': '../data/Lill_rowB.dat',
              'WTG': 'NREL5MW',
              'HH': 65.0,
              'R': 46.5,
              'stab': 'N',
              'accum': 'dominant',
              'optim': True}
    ##########################################################################
    return inputs

def getWT():
    inputs = getInputs()
    WT = wt.WindTurbine('Windturbine', '../WT-data/' + inputs['WTG'] + '/' + inputs['WTG'] + '_PC.dat', inputs['HH'],
                        inputs['R'])
    return WT

def getWF():
    inputs = getInputs()

    WT = getWT()

    WF = wf.WindFarm('Windfarm', inputs['WTcoord'], WT)

    return WF