import os

def threadbalancer(threads,jobtype=None):
    pass

def checkexists(path):
    path = os.path.abspath(path)
    if not os.path.is_dir(path):
        os.mkdir(path)

def procTitle(title):
    title = ' '+title+' '
    rows,cols = os.popen('stty size','r').read().split()
    print('')
    print(title.center(int(cols),'*'))

def mainTitle():
    print('''
 #####
#     #  #    #  #  ######  ######  #       ######   ####
#        ##   #  #  #       #       #       #       #
 #####   # #  #  #  #####   #####   #       #####    ####
      #  #  # #  #  #       #       #       #            #
#     #  #   ##  #  #       #       #       #       #    #
 #####   #    #  #  #       #       ######  ######   ####
    ''')
    print('\n')
