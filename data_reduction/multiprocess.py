from multiprocessing import Pool

def myProcess(iProcess):
    for i in runParams[iProcess]:
        print('Hello World: this is process ',iProcess,' with parameter ',i)

nRuns = 100
runParams = [range(i) for i in range(nRuns)]

nCores = 32
p = Pool(processes=nCores)
iProcess = range(nRuns)
p.map(myProcess,iProcess)
p.close()
