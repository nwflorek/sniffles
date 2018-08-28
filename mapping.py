import os;
import shlex
import subprocess as sub

def mapping(readData,runCFG,threads='1',ids=''):
    #inital parameters
    reference_sequence = os.path.abspath(runCFG['exec']['referenceSequence'])
    reference_sequence_name = os.path.basename(reference_sequence)
    libPath = runCFG['libPath']
    outDir = runCFG['exec']['outdir']

    #use id if provided otherwise get list
    if not ids:
        ids = readData.idList

    #index reference
    if 'bowtieindexed' not in runCFG:
        cmd = f'{libPath}/bin/bowtie2-build {reference_sequence} {reference_sequence_name}'
        cmd = shlex.split(cmd)
        sub.Popen(cmd,cwd=outDir).wait()
        readData.data['bowtieindexed'] = True

    

    #map reads
    for id in ids:
        print(readData.data)
        #determine reads
        if 'trimmed' in readData.data and id not in readData.data['mapProgress']['trimmed']:
            #TODO add status for unpaired read information
            reads = [readData.data['trimm'][id][0],readData.data['trimm'][id][1]]
            readData.data['mapProgress']['trimmed'].append(id)
            interleaved = False

        elif 'normalized' in readData.data and id not in readData.data['mapProgress']['normalized']:
            reads = readData.data['normalized'][id]
            readData.data['mapProgress']['normalized'].append(id)
            interleaved = True

        elif 'consensus' in readData.data and id not in readData.data['mapProgress']['consensus']:
            print('do things for mapping to consensus')

        else:
            print(f'There was an error mapping {id}.')
            exit()

        #determine interleaved or not
        #interleaved cmd
        if interleaved:
            interleaved_cmd = f"{libPath}/bin/bowtie2 -x {reference_sequence_name} --interleaved {reads} -S {id}_remapped.sam -p {threads} --local"
            cmd = shlex.split(interleaved_cmd)
        #split cmd
        else:
            split_cmd = f"{libPath}/bin/bowtie2 -x {reference_sequence_name} -1 {reads[0]} -2 {reads[1]} -S {id}.sam -p {threads} --local"
            cmd = shlex.split(split_cmd)

        #run command
        #sub.Popen(cmd,cwd=outDir).wait()
