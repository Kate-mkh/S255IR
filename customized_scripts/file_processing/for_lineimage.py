
# Каналы без линий
chans = ["" for x in range(11)]
chans[0] = '500~1000,1100~1400,1700~2000'
chans[1] = '358~1590'
chans[2] = '68~539, 889~1058, 1715~1867'
chans[3] = '1186~1426, 1487~1734'
chans[4] = '437~852, 1281~1600'
#chans[5] = '245~500'
chans[5] = '50~100, 200~230, 550~650'
chans[6] = '1065~1277, 1783~1874'
chans[7] = '107~276, 1866~2034'
chans[8] = '258~480, 1800~2000'
chans[9] = '220~400, 836~1100, 1875~2000'
chans[10] = '550~750, 1177~1400'




for i in range(11):
    print('cd /home/mikh_kate/kalenskii/CASA/s255ch'+ str(1 + i*2048)+'/')
    print('default(\'imcontsub\')')
    print('imagename = \'s255ch'+ str(1 + i*2048)+ '.image\'')
    print('linefile = \'s255ch'+ str(1 + i*2048)+ '.lineimage\'')
    print('contfile = \'s255ch'+ str(1 + i*2048)+ '.contimage\'')
    print('fitorder = 1')
    print('chans = \'' + chans[i] + '\'')
    print('stokes = \'I\'')
    print('imcontsub()')
    print('')
