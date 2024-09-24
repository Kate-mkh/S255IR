
print('class')
set_of_files = ['1', '2', '4', '6', '8', '10', '12', '14', '16', '18', '20', '1_2', '2_2', '4_2', '6_2', '8_2', '10_2', '12_2', '14_2', '16_2', '18_2', '20_2']

for i in range(len(set_of_files)):
    #print('cd /home/mikh_kate/kalenskii/CASA/s255ch'+str(1+ 2048*i)+'/')
    #print('class')
    print('fits read s255ch' + str(set_of_files[i]) + '_SMA2.fits')
    print('SET OBSERVATORY MAUNA_KEA')
    print('set variab position write')
    print('R%HEAD%POS%SYSTEM = 2')
    print('R%HEAD%POS%LAM = 1.6270837042')
    print('R%HEAD%POS%BET = 3.1399525289e-01')
    print('MODIFY DOPPLER 0')
    print('modi velo 7')
    print('set unit v f')
    print('pl')
    #print('number = ' + str(11+ 1+i) )
    print('number = ' + str(1+i) )   
    
    if i == 0:
        print('file out s255IR_SMA2 mult')
    else:
        print('file out s255IR_SMA2')  
    #print('file out s255IR_SMA1')  #
    print('write')
    #print('exit')
    #print(' ')
