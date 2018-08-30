import os

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
