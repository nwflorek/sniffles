#class for containing read information
class reads:
    #list contianing all isolate ids
    idList = []
    #data dictonary containing all of the run information
    #formatted id: ["read1 path","read2 path"]
    data = {}
    def __init__(self,path):
        #list of reads
        readList = []
        for root,dirs,files in os.walk(path):
            #scan path and look for fastq files and record ids
            for file in files:
                if '.fastq' in file:
                    if '_R1' in file or '_1' in file:
                        id = file.split('_')[0]
                        if id not in self.idList:
                            self.idList.append(id)
                    if '_R2' in file or '_2' in file:
                        id = file.split('_')[0]
                        if id not in self.idList:
                            self.idList.append(id)
                    readList.append(root+'/'+file)
        readList.sort()
        for id in self.idList:
            for read in readList:
                if id in read and '_R1' in read or '_1' in read:
                    self.data[id] = [read]
            for read in readList:
                if id in read and '_R2' in read or '_2' in read:
                    self.data[id].append(read)

    #return a list of id with each paired path as a 3 item sublist
    def retList(self,):
        l = []
        for id in self.idList:
            l.append([id,self.data[id][0],self.data[id][1]])
        return l
