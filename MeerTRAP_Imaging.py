# L. N. Driessen
# Last updated: 2021.01.25

import numpy as np
import os
import string
from astropy.time import Time
from astropy import units
import struct
from textwrap import wrap
import subprocess as sub
import glob


class schedFiles:
    '''
    Write files required for Sched and run Sched.
    
    Make sure you execute this in the folder you
    are running Sched in, otherwise it won't work.
    
    Args:
    outputpath (str): path to files. Should be the
                      path to the current working
                      directory
    SOURCE (str): a three-letter code for the source name.
                      eg. 'PKS', 'VEL'
    RA (str): the right ascension of the source (J2000)
              eg. '04:08:20.38'
    DEC (str): the declination of the source (J2000)
               eg. '-65:45:09.1'
    YEAR (str): the year of the observation
    MONTH (str): the month of the observation.
                 Include a 0 if the month is under 10.
                 eg. '07' for July, '12' for December
    DAY (str): the day of the observation.
               Include a 0 if the day is only a
               single digit. eg. '01'.
    START (str): the start time of the observation in UT.
                 eg. '09:47:41.835'
    DUR (str): the duration of the observation in seconds.
               As this is for Sched, the duration just needs to
               be longer than the overall observing time of
               the actual observation.
    MAIN_DUR (str): the total duration of the integration in seconds
    numchan (str): the number of channels (an integer)
    centre_freq (float): the actual centre frequency of your overall
                         observing band
    chanwidth (float): the actual width of each observing channel
    meekat_loc_file (str): path to and name of the csv file containing
                           the MeerKAT antenna names and Earth coordinates
    EOP_file (str): path to and name of the EOP.txt file to use for making
                    the v2d file
                    
    Kwargs:
    observing_antennas (list/array): list or array of the antennas
                                     used in the observation.
                                     e.g. np.arange(0, 64, 1) for
                                     all 64 MeerKAT antennas
                                     Default: np.arange(0, 64, 1)
    num_cores (int): number of cores on the machine you're using
                     Default: 1
    num_threads (int): number of threads to use per core
                       Default: 8
    '''
    def __init__(self, outputpath,
                 SOURCE, RA, DEC,
                 YEAR, MONTH, DAY,
                 START, DUR, MAIN_DUR,
                 numchan, centre_freq,
                 chanwidth,
                 meerkat_loc_file,
                 EOP_file,
                 observing_antennas=np.arange(0, 64, 1),
                 num_cores=1, num_threads=8):
        if outputpath[-1] != '/':
            outputpath += '/'
        self.outpath = outputpath
        
        self.source = SOURCE
        self.ra = RA
        self.dec = DEC
        self.year = YEAR
        self.month = MONTH
        self.day = DAY
        self.start = START
        self.dur = DUR
        self.main_dur = MAIN_DUR
        self.numchan_str = numchan
        
        self.cwidth = int(numchan)*chanwidth
        self.le = centre_freq-(self.cwidth/2) + (0.5*0.208984)
        self.lower_edge = str(self.le)
        self.chan_width = str(self.cwidth)
        self.sample_rate = str(2*self.cwidth)
        
        mkt_locs = np.genfromtxt(meerkat_loc_file,
                                 delimiter=',',
                                 dtype=str)

        self.mkt_codes = mkt_locs[:, 0]
        self.mkt_X = mkt_locs[:, 1]
        self.mkt_Y = mkt_locs[:, 2]
        self.mkt_Z = mkt_locs[:, 3]

        self.twoletternames = wrap(''.join(map(''.join,
                                               zip(('L'*26 +
                                                    'M'*26 +
                                                    'N'*(64-(26*2))),
                                                   ((string.ascii_uppercase)*3)[:64]))),
                                   2)
        antstr = np.insert(np.array(self.twoletternames)[observing_antennas],
                           range(1, len(np.array(self.twoletternames)[observing_antennas]), 1),
                           ',')
        antstr = np.insert(antstr, range(16, len(antstr), 16), '\n')
        self.antennastr = ''.join(antstr)

        self.num_ant = str(len(observing_antennas))
        
        self.observing_antennas = observing_antennas
        
        self.EOP = EOP_file
        
        self.num_cores = num_cores
        self.num_threads = num_threads


    def write_machines(self):
        '''
        Writes the machines file for Sched (mkt_1.machines)
        
        Uses num_threads, and outpath
        
        The machines file tells Sched how many machines are
        available to use
        '''
        machines_lines = np.repeat('localhost\n', self.num_threads)
        machinefilename = '{}mkt_1.machines'.format(self.outpath)
        with open(machinefilename, 'w') as f:
            for line in machines_lines:
                f.write(line)
        print('{} written'.format(machinefilename))


    def write_threads(self):
        '''
        Writes the threads file for Sched (mkt_1.threads)
        
        Uses num_cores, num_threads, and outpath
        
        The threads file tells Sched how many threads are
        available for use per machine you're using
        '''
        threads_lines = ['NUMBER OF CORES:    {}\n'.format(self.num_cores),
                         '{}\n'.format(self.num_threads)]
        threadfilename = '{}mkt_1.threads'.format(self.outpath)
        
        with open(threadfilename, 'w') as f:
            for line in threads_lines:
                f.write(line)
        print('{} written'.format(threadfilename))


    def make_station_file(self):
        '''
        Write the MeerKAT station file for Sched (stations_mkt.dat)

        This writes the station file for Sched using
        the MeerKAT antenna locations.
        '''
        station_lines = []
        for i in np.arange(0, 64, 1):
            station_template = ('STATION=MKT'+self.mkt_codes[i]+
                                '\n    STCODE='+self.twoletternames[i]+
                                '\n    CONTROL=VLBA'+
                                '\n    DBNAME='+self.mkt_codes[i]+'-MKT'+
                                '\n        MOUNT=ALTAZ AX1LIM=-90,450 AX2LIM=2.25,90'+
                                '\n        AX1RATE=83.6 AX2RATE=29.0'+
                                '\n        AX1ACC=0.75 AX2ACC=0.25'+
                                '\n        TSETTLE=6 TLEVSET=5 MINSETUP=5 DAR=RDBE2 NBBC=16'+
                                '\n        DISK=MARK5C   MEDIADEF=DISK    TSCAL=CONT'+
                                '\n        HOR_AZ =   0,  5, 10, 15, 25, 30, 40, 45, 70, 75,120,125,130,135,'+
                                '\n                 155,160,185,190,195,220,225,235,240,245,250,255,265,270,'+
                                '\n                 275,300,305,310,315,330,335,340,345,350,360'+
                                '\n        HOR_EL =   2,  2,  3,  2,  2,  3,  3,  4,  4,  5,  5,  4,  4,  3,'+
                                '\n                   3,  2,  2,  3,  4,  4,  3,  3,  4,  4,  5,  6,  6,  5,'+
                                '\n                   6,  6,  5,  6,  5,  5,  4,  4,  3,  2,  2'+
                                '\n        AXISOFF=  0.0'+
                                '\n        X='+self.mkt_X[i]+' Y='+self.mkt_Y[i]+' Z='+self.mkt_Z[i]+
                                '\n        DXDT= 0.0  DYDT=  0.0  DZDT= 0.0  EPOCH=54466'+
                                '\n        FRAME=\'GEO\''+
                                '\n    /'+
                                '\n')
            station_lines.append(station_template)
        stationfilename = '{}stations_mkt.dat'.format(self.outpath)
        
        with open(stationfilename, 'w') as f:
            for line in station_lines:
                f.write(line)

        print('{} written'.format(stationfilename))
    
    def make_freq_file(self):
        '''
        Write the frequency file for Sched (freq_mkt.dat)

        This frequency file now has the L-band
        and UHF frequency info in it.
        '''
        L_LINES = []
        UHF_LINES = []
        for i in np.arange(0, 64, 1):
            L_NAME = 'lband'
            L_STATION = 'MKT'+self.mkt_codes[i]
            L_PRIORITY = '1'
            L_RF1 = '856, 856'
            L_RF2 = '1712, 1712'
            L_IFNAME = 'A, C' # This is what Adam has for ASKAP
            L_FE = '\'20cm\', \'20cm\''
            L_POL = 'RCP, LCP' # Work out what linear polarisation is for here
            L_LO1 = '2100, 2100'
            syn = '2.1'
            l_template = ('Name='+L_NAME+
                          '\n    Station='+L_STATION+
                          '\n    Priority='+L_PRIORITY+
                          '\n    rf1='+L_RF1+
                          '\n    rf2='+L_RF2+
                          '\n    ifname='+L_IFNAME+
                          '\n    fe='+L_FE+
                          '\n    pol='+L_POL+
                          '\n    lo1='+L_LO1+
                          '\n    syn(2)='+syn+
                          '\n/'+
                          '\n')
            L_LINES.append(l_template)
            
            UHF_NAME = 'uhfband'
            UHF_STATION = 'MKT'+self.mkt_codes[i]
            UHF_PRIORITY = '1'
            UHF_RF1 = '580, 580'
            UHF_RF2 = '1015, 1015'
            UHF_IFNAME = 'A, C' # This is what Adam has for ASKAP
            UHF_FE = '\'20cm\', \'20cm\''
            UHF_POL = 'RCP, LCP' # Work out what linear polarisation is for here
            UHF_LO1 = '2100, 2100'
            syn = '2.1'
            UHF_template = ('Name='+UHF_NAME+
                          '\n    Station='+UHF_STATION+
                          '\n    Priority='+UHF_PRIORITY+
                          '\n    rf1='+UHF_RF1+
                          '\n    rf2='+UHF_RF2+
                          '\n    ifname='+UHF_IFNAME+
                          '\n    fe='+UHF_FE+
                          '\n    pol='+UHF_POL+
                          '\n    lo1='+UHF_LO1+
                          '\n    syn(2)='+syn+
                          '\n/'+
                          '\n')
            UHF_LINES.append(UHF_template)

        freqfilename = '{}freq_mkt.dat'.format(self.outpath)
        with open(freqfilename, 'w') as f:
            for line in L_LINES+UHF_LINES:
                f.write(line)
        print('{} written'.format(freqfilename))


    def make_main_file(self):
        '''
        Write the main file for Sched (mkt.key)

        The main file gives Sched the instructions
        for what it needs to do and the files it
        needs to make.
        '''
        # Pre-text about who is doing
        # the observing etc
        # IMPORTANT: obstype must be
        # VLBI
        main1 = ('version  = 1'+'\n'+
                 'expt     = \'meerkat\''+'\n'+
                 'expcode  = \'mkt\''+'\n'+
                 'obstype  = \'VLBI\''+
                 'piname   = \'L. N. Driessen\''+'\n'+
                 'address1 = \'U of Manchester\''+'\n'+
                 'address2 = \'\''+'\n'+
                 'address3 = \'United Kingdom\''+'\n'+
                 'email    = \'Laura.Driessen@postgrad.manchester.ac.uk\''+'\n'+
                 'phone    = \'+44 7480 324 387 (w)\''+'\n'+
                 'obsphone = \'+44 7480 324 387 (w)\''+'\n'+
                 'obsmode  = \'MeerKAT 4096 channel\''+'\n'+
                 'note1    = \'MeerKAT\''+'\n'+
                 '\n')

        # Information about the correlator
        # you'll be using and how averaged
        # the channels will be etc.
        # Just copied from Adam's script,
        # nothing to mess around with here.
        correlator = ('! ================================================================'+'\n'+
                      '!       Correlator section'+'\n'+
                      '! ================================================================'+'\n'+
                      'correl   = \'Socorro\''+'\n'+
                      'coravg   = 2'+'\n'+
                      'corchan  = 16'+'\n'+
                      'cornant  = '+self.num_ant+'\n'+
                      'corpol   = \'on\''+'\n'+
                      'corwtfn  = \'uniform\''+'\n'+
                      'corsrcs  = \'from .sum file only\''+'\n'+
                      'cortape  = \'ftp\''+'\n'+
                      'cornote1 = \'This is special MeerKAT correlation\''+'\n'+
                      '\n')
        
        # These params are for fudging Sched and WILL NOT change
        L_setinit_nchan = '2'
        L_setinit_bbfilt = '128.0' # this is where the BW is coming from!
        L_setinit_netside = 'U'
        L_setinit_bits = '2'
        L_setinit_fLO = '2100.0'
        L_setinit_fREF = '2100.0'
        L_setinit_fOFF = '-608.0, -608.0'
        
        # This is where the main setup
        # and calls to the catalogues are.
        # This also calls to the centre "source"
        # setup earlier.
        # format = 'VDIF' is important here.
        catalogues = ('! ================================================================'+'\n'+
                      '!       Catalogs (special meerkat versions)'+'\n'+
                      '! ================================================================'+'\n'+
                      'stafile  = {}stations_mkt.dat'.format(self.outpath)+'\n'+
                      'freqfile = {}freq_mkt.dat'.format(self.outpath)+'\n'+
                      'overwrite'+'\n'+
                      'srccat /'+'\n'+
                      'EQUINOX = J2000'+'\n'+
                      'SOURCE  =\''+self.source+'\' RA='+self.ra+' DEC='+self.dec+''+'\n'+
                      'REMARKS =\'Beam centre for dumped voltage data\' /'+'\n'+
                      'endcat /'+'\n'+
                      ''+'\n'+
                      'setinit = MK_L.set /'+'\n'+
                      ' dbe      = \'rdbe_ddc\''+'\n'+
                      ' format   = \'vdif\''+'\n'+
                      ' nchan    = '+L_setinit_nchan+'\n'+
                      ' bbfilt   = '+L_setinit_bbfilt+'\n'+
                      ' netside  = '+L_setinit_netside+'\n'+
                      ' bits     = '+L_setinit_bits+'\n'+
                      ' firstlo  = '+L_setinit_fLO+'\n'+
                      ' freqref  = '+L_setinit_fREF+'\n'+
                      ' freqoff  = '+L_setinit_fOFF+'\n'+
                      ' pcal     = \'off\''+'\n'+
                      '   /'+'\n'+
                      'endset /'+'\n'+
                      ''+'\n'+
                      '\n')
        # Setup the source
        # and the observing scan
        source_setup = ('! ================================================================'+'\n'+
                        '!       Source setup'+'\n'+
                        '! ================================================================'+'\n'+
                        'year     = '+self.year+'\n'+
                        'month    = '+self.month+'\n'+
                        'day      = '+self.day+'\n'+
                        'start    = '+self.start+'\n'+
                        'stations = '+self.antennastr+'\n'+
                        'setup    = \'MK_L.set\''+'\n'+
                        'minpause = 5'+'\n'+
                        ''+'\n'+
                        'source = \''+self.source+'\'  dur = '+self.main_dur+'  gap = 0   /'+'\n'+
                        '\n')

        mainfile_text = (main1 +
                         correlator +
                         catalogues +
                         source_setup)
        mainfilename = '{}mkt.key'.format(self.outpath)
        with open(mainfilename, 'w') as f:
            for line in mainfile_text:
                f.write(line)

        print('{} written'.format(mainfilename))
    
    def make_v2d_file(self):
        '''
        Write the v2d file for DiFX (mkt.v2d)

        The v2d file is needed for vex2difx
        '''
        vexfile = '{}mkt.vex'.format(self.outpath)
        
        # Get the format of the
        # antenna string in the vexfile
        try:
            vexantennastr = self.antennastr.replace('\n', '')
        except ValueError:
            vexantennastr = self.antennastr

        ### Get the EOPs ###
        firstday = '%s-%s-%sT%s' % (self.year, '01', '01', '00:00:00.00')
        time = '%s-%s-%sT%s' % (self.year, self.month, self.day, self.start)
        t = Time(time, format='isot', scale='utc')
        ft = Time(firstday, format='isot', scale='utc')

        DOY = str(int(t.mjd -ft.mjd)-2)

        try:            
            with open((self.EOP).format(self.outpath)) as f:
                lines = f.readlines()
            print('Opening the backup EOP.txt file')
        except FileNotFoundError:
            print('Backup EOP not found *gasp*')

        TAI = []
        UT1 = []
        xwobble = []
        ywobble = []
        for l in lines:
            if 'TAI-UTC' in l:
                TAI.append(l)
            if 'ut1-utc' in l:
                UT1.append(l)
            if 'x_wobble' in l:
                xwobble.append(l)
            if 'y_wobble' in l:
                ywobble.append(l)

        eops = []
        for i, sub in enumerate(np.arange(-2, 3, 1)):
            mjd = int(t.mjd)+sub
            ut1 = UT1[i].split()[2]
            tai = TAI[i].split()[1]
            xw = xwobble[i].split()[2]
            yw = ywobble[i].split()[2]

            eop = ('EOP %s { xPole=%s yPole=%s tai_utc=%s ut1_utc=%s }\n' % (mjd, xw, yw, tai, ut1))
            eops.append(eop)

        # The v2d file preamble
        v2d_preamble = ['#  Template v2d file for DiFX correlation of mkt',
                        '\n\n',
                        'vex = '+vexfile,
                        '\n',
                        'minLength=1',
                        '\n',
                        'antennas='+vexantennastr,
                        '\n',
                        'tweakIntTime=True',
                        '\n\n']

        # Set up the basic info for each
        # antenna here. The 'format=VDIF/5032/2'
        # bit ensures that you don't get that
        # INTERLACEDVDIF error thing.
        # The files don't exist (hence
        # file_fake.vdif)
        antenna_info = []

        for ant in self.observing_antennas:
            code = self.twoletternames[ant]
            dummy_name = code+'_fake.vdif'
            info = ('ANTENNA %s {file=%s format=VDIF/5032/2 '
                    'clockOffset=0 clockRate=0 '
                    'clockEpoch=57000.0}\n' % (code, dummy_name))
            antenna_info.append(info)
        antenna_info.append('\n')
        
        # The v2d postamble
        v2d_postamble = [('# The nChan should never be less than 128.'+'\n'+
               '# For numbers of channels < 128, set specAvg so nChan/specAvg'+'\n'+
               '# gives the desired number of channels'+'\n'+
               'SETUP default'+'\n'+
               '{'+'\n'+
               '  tInt =  '+self.dur+'\n'+
               '  nFFTChan =    '+self.numchan_str+'\n'+
               '  nChan =  '+self.numchan_str+'\n'+
               '  doPolar = True # Full stokes'+'\n'+
               '}'+'\n'+
               ''+'\n'+
               '# This, along with SETUP default above, should always be done'+'\n'+
               'RULE default'+'\n'+
               '{'+'\n'+
               '  setup = default'+'\n'+
               '}'+'\n'+
               ''+'\n'+
               '#  SETUP place holders (commented)'+'\n'+
               '# SETUP MK_L.set {}'+'\n'+
               ''+'\n'+
               '# Sources (pointing centers) with recorded data but no offset pointing centers:'+'\n'+
               'SOURCE '+self.source+' { }'+'\n'+'\n')]
        
        v2d_full = (v2d_preamble+antenna_info+
                    v2d_postamble+eops)
        
        v2dfilename = '{}mkt.v2d'.format(self.outpath)
        with open(v2dfilename, 'w') as f:
            for line in v2d_full:
                f.write(line)

        print('{} written'.format(v2dfilename))
        
    def edit_vex(self):
        '''
        Replace dummy values in the vex file with real values

        For Sched to run correctly some of the
        original values have to be fudged. Here
        we change them to the correct values in
        the vex file to run DiFX correctly.
        '''
        os.system(('cp {0}mkt.vex '
                   '{1}mkt_original.vex').format(self.outpath,
                                                 self.outpath))
        
        filein = '{0}mkt.vex'.format(self.outpath)
        fileout = '{0}mkt.vex'.format(self.outpath)
        f = open(filein,'r')
        filedata = f.read()
        f.close()
        newdata = filedata.replace('1492.00MHz',
                                   '{}MHz'.format(self.lower_edge))
        newdata = newdata.replace('1492.00 MHz',
                                  '{}  MHz'.format(self.lower_edge))
        newdata = newdata.replace('2x128MHz',
                                  '{}0.208984MHz'.format(self.numchan_str))
        newdata = newdata.replace('128.00 MHz',
                                  '0.0208984 MHz')
        newdata = newdata.replace('256.000 Ms/sec',
                                  '{} Ms/sec'.format(self.sample_rate))
        
        f = open(fileout,'w')
        f.write(newdata)
        f.close()
        print(('New vex file written with '
               'values replaced: {}').format(fileout))


    def generate(self):
        '''
        Make all the files for Sched, then run Sched,
        vex2difx, and calcif2. This does the correction
        for the geometric delays between antennas.
        
        Here we also check if we're in the correct
        directory.
        
        So that you don't run anything unecessarily, also
        checks for whether the files you need already
        exist. If they don't exist, make them.
        '''
        if os.getcwd()+'/' == self.outpath:
            print('Running in {}'.format(self.outpath))
            if not os.path.exists('{}mkt_1.machines'.format(self.outpath)):
                self.write_machines()
            else:
                print('Machines file exists')

            if not os.path.exists('{}mkt_1.threads'.format(self.outpath)):
                self.write_threads()
            else:
                print('Threads file exists')

            if not os.path.exists('{}stations_mkt.dat'.format(self.outpath)):
                self.make_station_file()
            else:
                print('Station file exists')

            if not os.path.exists('{}freq_mkt.dat'.format(self.outpath)):
                self.make_freq_file()
            else:
                print('Frequency file exists')

            if not os.path.exists('{}mkt.key'.format(self.outpath)):
                self.make_main_file()
            else:
                print('Main file exists')

            if not os.path.exists('{}mkt.v2d'.format(self.outpath)):
                self.make_v2d_file()
            else:
                print('v2d file exists')

            # Run Sched
            print('\nRunning Sched')
            try:
                sub.check_output('sched < mkt.key', shell=True)
            except sub.CalledProcessError as e:
                if '+++ERROR+++ SCHED terminating' in (e.output).decode('utf-8'):
                    print(('-------------------------------\n'
                           '+++ERROR+++ SCHED terminating\n'
                           'Sched really failed\n'
                           '-------------------------------\n'))
                    sched_completed = False
                    print((e.output).decode('utf-8'))
                else:
                    print(('-------------------------------\n'
                           'Sched run successfully\n'
                           '-------------------------------\n'))
                    sched_completed = True
            if sched_completed:
                self.edit_vex()
                # Run vex2difx
                print('Running vex2difx')
                try:
                    sub.check_output(['vex2difx', 'mkt.v2d'])
                    print(('-------------------------------\n'
                           'vex2difx run successfully\n'
                           '-------------------------------\n'))

                    # Run calcif2
                    print('Running calcif2')
                    try:
                        sub.check_output(['calcif2', '-v', 'mkt_1.calc'])
                        print(('-------------------------------\n'
                               'calcif2 run successfully\n'
                               '-------------------------------\n'))
                    except sub.CalledProcessError as c:
                        print(('-----------------------------------\n'
                               'Process killed, calcif2 did not run\n'
                               '-----------------------------------'))
                        print(c.output)
                except sub.CalledProcessError as v:
                    print(('-----------------------------------\n'
                           'Process killed, vex2difx did not run\n'
                           '-----------------------------------'))
                    print(v.output)
        else:
            print('Running in the incorrect directory')
            print('Need to run in the output directory: {}'.format(self.outpath))


class xgpuToDiFX:
    '''
    Convert the XGPU data file into a DiFX file
    and use DiFX to convert that to a fits file.
    
    This class will take the input xGPU correlation
    file and convert that into the correct DiFX
    format. Then it writes that to VDIF format,
    to then run DiFX on it to create the final
    visibilities.
    
    Args:
    xGPU_filename (str): name and path to the xGPU
                         file
    num_chan (str/int): number of channels
    inttime (str/float): integration time for a single
                         timestep in seconds
    date (str): string of format 'YYYY-MM-DD' for the date
                of the observation
    starttime (str): string of format 'HH:MM:SS.ssssssss'
                     which as much precision as possible,
                     the exact start time of the observation
    final_path (str): the path to the final DiFX file
    
    Kwargs:
    num_ant (str/int): number of antennas included
    num_pol (str/int): number of polarisations
    '''
    def __init__(self, xGPU_filename,
                 num_chan, inttime,
                 date, starttime,
                 final_path,
                 num_ant=64,
                 num_pol=4,
                 observing_antennas=np.arange(0, 64, 1)):
        if final_path[-1] != '/':
            final_path += '/'
        
        self.num_ant = int(num_ant)
        print('Number of antennas: {}'.format(self.num_ant))
        self.num_chan = int(num_chan)
        print('Number of channels: {}'.format(self.num_chan))
        self.num_pol = int(num_pol)
        print('Number of polarisations: {}'.format(self.num_pol))
        self.inttime = float(inttime)
        print('Int time: {}'.format(self.inttime))
        
        f = open(xGPU_filename, 'rb')
        f.seek(4096, 0)
        data = np.fromfile(f, dtype=np.complex64)
        data = np.conjugate(data)
        f.close()
        
        self.num_baselines = int(self.num_ant *
                                 (self.num_ant+1)/2)
        self.vals_per_channel = (self.num_baselines *
                                 self.num_pol)
        self.vals_per_timestep = (self.vals_per_channel *
                                  self.num_chan)
        self.num_timesteps = int(len(data)/
                                 self.vals_per_timestep)
        self.total_length = (self.num_baselines *
                             self.num_pol *
                             self.num_chan *
                             self.num_timesteps)

        self.data = np.copy(data[:self.total_length])
        
        self.antenna_order = observing_antennas
        
        data = None # Free up memory

        print('xGPU file information\n---------------------------------')
        print('Number of baselines: {}'.format(self.num_baselines))
        print('Values per channel: {}'.format(self.vals_per_channel))
        print('Values per time step: {}'.format(self.vals_per_timestep))
        print('Number of time steps: {}'.format(self.num_timesteps))
        
        ### Info for WriteDiFX ###

        # The polarisations in the correct
        # byte format required by the DiFX header
        self.pol = [b'RR', b'RL', b'LR', b'LL']

        # These are all parameters needed by the
        # DiFX header, but they're the same everytime
        # and don't matter to us. So I set them once, here.
        self.h1 = (int.to_bytes(4278255360,
                                byteorder='little',
                                length=4) + 
                   int.to_bytes(1,
                                byteorder='little',
                                length=4))
        self.h2 = (int.to_bytes(0,
                                byteorder='little',
                                length=4) +
                   int.to_bytes(0,
                                byteorder='little',
                                length=4) +
                   int.to_bytes(0,
                                byteorder='little',
                                length=4))
        self.h3 = (int.to_bytes(0,
                                byteorder='little',
                                length=4) +
                   struct.pack('<d', 1) +
                   struct.pack('<d', 0) +
                   struct.pack('<d', 0) +
                   struct.pack('<d', 0))
        
        self.date = date
        self.start = starttime

        ### Set up the DiFX output file name by getting the time
        ### and date of the start of the observation.
        time = Time('{0}T{1}'.format(self.date, self.start),
                    format='isot', scale='utc')
        MJD = int(time.mjd)
        self.mjd = MJD
        self.all_dt_seconds = (time +
                          ((self.inttime*(0.5 + (np.arange(0, self.num_timesteps))))*units.s) -
                          Time('{0}T{1}'.format(self.date, '00:00:00.0'),
                               format='isot', scale='utc')).sec
        self.mjd_bytes = int.to_bytes(MJD, byteorder='little',
                                 length=4)

        os.system('mkdir {}mkt_1.difx/'.format(final_path))
        self.filename = (('{0}mkt_1.difx/DIFX_{1}_{2}.'
                          's0000.b0000').format(final_path,
                                                MJD,
                                                int(self.all_dt_seconds[0])))
        print('DiFX output filename: ', self.filename)

    def make_baseline_map(self):
        '''
        Make a map of which baseline matches
        which position in a polarisation triangle.

        This makes a lower triangle that matches the
        positions of the xGPU lower triangle to
        the baseline index.
        '''

        baseline_map = np.zeros((self.num_ant,
                                 self.num_ant),
                                dtype=int) * np.nan
        baseline_num = 0
        for ant1 in range(self.num_ant):
            for ant2 in range(self.num_ant-ant1):
                baseline_map[ant2+ant1, ant1] = baseline_num
                baseline_num += 1

        baseline_index = []
        for a1 in range(self.num_ant):
            for a2 in range(a1+1):
                baseline_index.append(baseline_map[a1, a2])

        return np.array(baseline_index).astype(int), np.array(baseline_map).astype(int)

    def get_baseline_num(self, ant_1, ant_2):
        '''
        Defines the old-school baseline number.

        Takes the number of each antenna and converts
        it into the old-school format for baseline labels,
        a 4 bit little endian integer that is given by
        256*(ant_1+1) + (ant_2+1)
        '''
        baseline = int(256*(ant_1+1) + (ant_2+1))
        return int.to_bytes(baseline,
                            byteorder='little',
                            length=4)

    def generate_difxinput(self):
        '''
        Does the conversion between xGPU and DIFX format
        
        Conversion for the data array only
        '''
        baseline_index, baseline_map = self.make_baseline_map()

        bl_val = baseline_index * self.num_chan * self.num_pol

        bl_level0 = (np.repeat(bl_val, self.num_pol) +
                     np.tile((np.arange(0, self.num_pol, 1)*
                              self.num_chan),
                             self.num_baselines))
        bl_level1 = np.tile(bl_level0, self.num_chan)
        bl_level2 = np.tile(bl_level1, self.num_timesteps)

        ch_level1 = np.repeat(np.arange(0, self.num_chan, 1),
                              (self.num_baselines *
                               self.num_pol))
        ch_level2 = np.tile(ch_level1, self.num_timesteps)

        t_level2 = np.repeat((np.arange(0, self.num_timesteps, 1) *
                              self.vals_per_timestep),
                             self.vals_per_timestep)

        index_map = bl_level2 + ch_level2 + t_level2

        difx_input = np.zeros(self.total_length,
                              dtype=np.complex64)
        difx_input[index_map] = self.data[:self.total_length]
        
        return difx_input
    
    def generate(self):
        '''
        Performs all the steps to convert the xGPU data
        file in the DiFX data file, and then runs difx2fits
        '''
        difx_input = self.generate_difxinput()

        ###### Write out the final DiFX output ######

        # Open up the file as a binary write file
        with open(self.filename, 'wb') as binfile:
            # We need a header for each timestep,
            # baseline, and polarisation. The outer
            # loop being the timestep. So loop through
            # the timesteps, getting the centre time
            # of each one to put in the header.
            for stamp in range(self.num_timesteps):
                seconds = struct.pack('<d', self.all_dt_seconds[stamp])

                # The pointer for where you are in the output
                # array (populated in the giant loop above),
                # and start at baseline 0-0.
                pointer = 0
                a2_start = 0
                # Loop fully through all the antennas as 
                # antenna one
                for a1_index in range(self.num_ant):
                    # Only get non-redundant baselines, so
                    # loop through the second antenna differently
                    a1 = self.antenna_order[a1_index]
                    for a2_index in range(self.num_ant)[a2_start:]:
                        a2 = self.antenna_order[a2_index]
                        # Loop through the four polarisations
                        for polarisation in self.pol:
                            # Get the baseline number according to
                            # the strange old-school notation
                            baselinenum = self.get_baseline_num(a1, a2)

                            # The header in the correct order
                            header = (self.h1+
                                      baselinenum+
                                      self.mjd_bytes+
                                      seconds+
                                      self.h2+
                                      polarisation+
                                      self.h3)
                            # The bit of the data that matches the header
                            standard = stamp*self.num_chan*self.num_baselines*self.num_pol
                            data_chunk = difx_input[(pointer*self.num_chan +
                                                    standard):
                                                    ((pointer+1)*self.num_chan +
                                                    standard)]
                            # Go to the next data chunk
                            pointer += 1
                            # Write the header and data_chunk to the file
                            binfile.write(header)
                            binfile.write(data_chunk)
                    # Skip to the next antenna pair (missing redundant ones)
                    a2_start += 1
        print('Running difx2fits')
        try:
            sub.check_output('difx2fits mkt_1 >\'mkt_difx2fits.log\'', shell=True)
            
            print(('-------------------------------\n'
                   'difx2fits run successfully\n'
                   '-------------------------------\n'))
        except sub.CalledProcessError as d:
            print(('###########################\n'
                   'difx2fits failed'
                   '###########################\n'))
            print(d.output)


class createMeasurementSet:
    '''
    Uses CASA to convert the final difx2fits file
    into a measurement set
    
    The measurement set can then be cleaned by WSClean
    
    Args:
    flag_antennas (str): string of the format 'LA&LI&NJ'
                         listing all of the antennas missing
                         (or really, fudged) in the xGPU
                         correlation
    measurementset (str): the name and path to the final
                          final measurement set
    path (str): path to the MKT.0.bin0000.source0000.FITS
                file made by difx2fits
    '''
    def __init__(self,
                 flag_antennas,
                 measurementset,
                 path):
        if path[-1] != '/':
            path += '/'
        self.fitsname = '{}MKT.0.bin0000.source0000.FITS'.format(path)
        self.ms = measurementset
        self.flagant = flag_antennas
        
    def process(self):
        '''
        Convert the FITS file to a measurement set
        and flag the antennas that are missing
        '''
        if self.flagant == '':
            casa_script = [(('importfitsidi(fitsidifile=\'{0}\', '
                             'vis=\'{1}\')\n').format(self.fitsname,
                                                      self.ms)),
                           (('vishead(vis=\'{0}\', '
                             'mode=\'put\', hdkey=\'telescope\', '
                             'hdvalue=\'MeerKAT\')\n').format(self.ms))] 
        else:
            casa_script = [(('importfitsidi(fitsidifile=\'{0}\', '
                             'vis=\'{1}\')\n').format(self.fitsname,
                                                      self.ms)),
                           (('vishead(vis=\'{0}\', '
                             'mode=\'put\', hdkey=\'telescope\', '
                             'hdvalue=\'MeerKAT\')\n').format(self.ms)),
                           (('flagdata(vis=\'{0}\', '
                             'antenna=\'{1}\', '
                             'action=\'apply\')\n').format(self.ms,
                                                           self.flagant))]
        with open(('casa_script.py'), 'w') as f:
            for line in casa_script:
                f.write(line)
        print('Running casa script')
        os.system('casa -c \'casa_script.py\'')
        print('Casa script complete')


class getObsInformation:
    '''
    Read the xGPU header and manually entered information
    into the correct formats
    
    Args:
    dat_file (str): filename and path to the xGPU file
    source (str): capitalised three letter name for the source
                  e.g. 'VEL' for Vela
    ra (str): right ascension of the source in the format
              'hh:mm:ss.ssss'
    dec (str): declination of the source in the format
               'dd:mm:ss.ssss'
    time_step_size (float): the size of each individual time
                            step of the observation (as in the
                            xGPU file)
                            
    kwargs:
    current_folder (str): the path to the foler
                          that you're working in
                          Default: '/raid/driessen/Correlation/'
    casa_path (str): the path where you'd like CASA to work
                     Default: /raid/driessen/Correlation/'
    output_path (str): the path where you would like the
                       final products saved to
                       Default: '/raid/driessen/Correlation/'
    mkt_loc_file (str): the path to the MeerKAT telescope
                        location file
                        Default: ('/raid/driessen/'
                                  'Correlation/'
                                  'meerkat_dish_positions.csv')
    print_info (bool): if True, will print information about the
                       observation
                       Default: True
    antenna_labels (str): The path to the file containing the antenna
                          label codes
                          Default: ('/raid/driessen/'
                                    'Correlation/'
                                    'antenna_labels.csv')
    antennas_used (array): the names (e.g. 'm001') or numbers (e.g. 1)
                           of the antennas used in the observation.
                           If set to 'None', will read the antennas directly
                           from the xGPU header
                           Default: 'None'
    antenna_order (array): the order that the antennas_used are in. Will also
                           apply this order to the antennas read from the xGPU
                           head if antennas_used='None'. If antenna_order='None',
                           no re-ordering will be applied.
                           Default: 'None'
    
    '''
    def __init__(self, dat_file, source, ra, dec, time_step_size,
                 current_folder='/raid/driessen/Correlation/',
                 casa_path='/raid/driessen/Correlation/',
                 output_path='/raid/driessen/Correlation/',
                 mkt_loc_file=('/raid/driessen/'
                               'Correlation/'
                               'meerkat_dish_positions.csv'),
                 print_info=True,
                 antenna_labels=('/raid/driessen/'
                                 'Correlation/'
                                 'antenna_labels.csv'),
                 antennas_used='None',
                 antenna_order='None'):
        
        antenna_dict = np.genfromtxt(antenna_labels,
                                     delimiter=',',
                                     dtype=str)
        antenna_dictionary = dict()
        for row in antenna_dict:
            try:
                key = eval(row[0])
            except NameError:
                key = row[0]
            antenna_dictionary[key] = row[1]
        self.antenna_dictionary = antenna_dictionary

        self.dat_file = dat_file
        self.filename_base = self.dat_file.split('/')[-1].split('.')[0]
        self.measurementset = (('{0}{1}.ms').format(current_folder,
                                                    self.filename_base))
        self.SOURCE = source
        self.RA = ra
        self.DEC = dec
        self.time_step_size = time_step_size
        
        self.casa_path = casa_path
        self.output_path = output_path
        self.mkt_loc_file = mkt_loc_file

        with open(self.dat_file, 'rb') as f:
            count = 0
            while count < 49:
                line = f.readline().decode("utf-8")
                count += 1

                if 'IDX2_LIST' in line:
                    line = line.strip()
                    line = line.split()
                    antennas = np.array(line[1].split(',')).astype(int)
                if 'IDX3_LIST' in line:
                    line = line.strip()
                    line = line.split()
                    dt = line[-1]
                    dt = dt.split('-')
                    self.YEAR = dt[0]
                    self.MONTH = dt[1]
                    self.DAY = dt[2]
                    self.START = dt[3]
                if ('FREQ' in line) and ('REFFREQ' not in line):
                    line = line.strip()
                    line = line.split()
                    self.centre_freq = float(line[1])*1e-6
                if 'NCHANS' in line:
                    if 'TOTAL' in line:
                        line = line.strip()
                        line = line.split()
                        self.total_nchans = int(line[1])
                    else:
                        line = line.strip()
                        line = line.split()
                        self.numchan = line[1]
            self.chanwidth = 856. / self.total_nchans

        if print_info:
            if antennas_used != 'None':
                print('Antennas included: {}\n'.format(antennas_used))
                print('(Antennas manually listed)')
            else:
                print('Antennas included: {}\n'.format(antennas))
                print('(Antennas read from xGPU header)')
            print(('YEAR {0} MONTH {1} '
                   'DAY {2} START {3}\n').format(self.YEAR,
                                                 self.MONTH,
                                                 self.DAY,
                                                 self.START))
            print('Centre freq {} MHz\n'.format(self.centre_freq))
            print('Total number of channels in obs: {0}\n'.format(self.total_nchans))
            print('Number of channels in this file: {}\n'.format(self.numchan))
            print('Channel width: {} MHz\n'.format(self.chanwidth))
        

        self.total_antennas = np.arange(0, 64, 1)

        if antennas_used != 'None':
            if isinstance(antennas_used[0], str):
                antennas_used_number = []
                for ant in antennas_used:
                    ant1 = ant.strip('m')
                    antennas_used_number.append(int(ant1))
                antennas = antennas_used_number
            else:
                antennas = list(antennas_used)
        else:
            antennas = list(antennas)
            
        if antenna_order != 'None':
            antennas = np.array(antennas)[antenna_order]
            antennas = list(antennas)

        missing_antennas = ''
        for ante in self.total_antennas:
            if ante not in antennas:
                missing_antennas += antenna_dictionary[ante]
                missing_antennas += ','
                antennas.append(ante)
        missing_antennas = missing_antennas[:-1]
        self.missing_antennas = missing_antennas
        if print_info:
            print(('Antennas missing '
                   'from obs: {}').format(self.missing_antennas))
            print(('Final antenna list: {}').format(antennas))
        self.antennas = np.array(antennas)

        # Number of antennas
        self.numant = 64

        # The real integration time of each time step
        # in the observation
        self.DUR = str(time_step_size)
        if print_info:
            print('Time step size: {}'.format(self.DUR))

        # The fudged total integration time for Sched
        # This meeds to be at least as long as the actual
        # full integration time. It's good if it can
        # actually be the real length of the observations
        self.MAIN_DUR = str(8000.*0.00489988785046528)

        # Number of polarisations
        self.numpol = 4
        # Date for xGPU file writing
        self.date = '{0}-{1}-{2}'.format(self.YEAR,
                                         self.MONTH,
                                         self.DAY)

        # Antennas to flag
        self.flag_antennas = missing_antennas
        
        start_split = self.START.split(':')
        start_sec = int(float(start_split[2]))
        if start_sec < 10:
            start_sec = '0{0}'.format(start_sec)
        else:
            start_sec = str(start_sec)
        start_rd = ('{0}:{1}:{2}').format(start_split[0],
                                          start_split[1],
                                          start_sec)
        self.START_rd = start_rd
        print('Start rounded down: {}'.format(start_rd))


if __name__ in '__main__:
    backup_EOP = ('/path/to/backupEOPfile/'
                  'backup_EOP.txt')
    dat_files = sorted(glob.glob('/path/to/data/files/*.corr'))

    for dat_file in dat_files:
        bf = dat_file.split('/')[-1].split('.')[0]
        folder_name = '{}_Results'.format(bf)
        folder_exist = glob.glob(folder_name)

        if len(folder_exist) == 0:
            print('###########################\n')
            print('Running on {}\n'.format(bf))
            print('###########################\n')
            obs_info = getObsInformation(dat_file,
                                         'SPF',
                                         '04:52:34.11',
                                         '-17:59:23.4',
                                         0.00128,
                                         current_folder='/path/to/folder/',
                                         casa_path='/path/to/folder/',
                                         output_path='/path/to/folder/',
                                         mkt_loc_file=('/path/to/folder/'
                                                       'meerkat_dish_positions.csv'),
                                         print_info=True,
                                         antenna_labels='antenna_labels.csv',
                                         antennas_used = ['m001', 'm005', 'm008', 'm009',
                                                          'm010', 'm012', 'm015', 'm016',
                                                          'm017', 'm018', 'm019', 'm020',
                                                          'm022', 'm024', 'm026', 'm027',
                                                          'm029', 'm031', 'm032', 'm033',
                                                          'm034', 'm035', 'm036', 'm038',
                                                          'm039', 'm040', 'm041', 'm042',
                                                          'm044', 'm046', 'm047', 'm048',
                                                          'm050', 'm051', 'm055', 'm056',
                                                          'm057', 'm058', 'm060', 'm061',
                                                          'm062', 'm063'],
                                         antenna_order=[32,33,34,35,36,37,
                                                        38,39,40,41,19,28,
                                                        29,31,0,1,2,3,4,5,
                                                        6,7,8,9,10,11,12,13,
                                                        14,15,16,17,18,20,21,
                                                        22,23,24,25,26,27,30])

            print('###########################\n')
            print('Running Sched')
            print('###########################\n')
            schedFiles(obs_info.output_path,
                       obs_info.SOURCE,
                       obs_info.RA,
                       obs_info.DEC,
                       obs_info.YEAR, 
                       obs_info.MONTH,
                       obs_info.DAY,
                       obs_info.START_rd,
                       obs_info.DUR,
                       obs_info.MAIN_DUR,
                       obs_info.numchan,
                       obs_info.centre_freq,
                       obs_info.chanwidth,
                       obs_info.mkt_loc_file,
                       backup_EOP).generate()

            print('###########################\n')
            print('Running xGPU to DiFX')
            print('###########################\n')
            # observing_antennas should be the actual order
            # of the antennas
            xgpuToDiFX(obs_info.dat_file,
                       obs_info.numchan,
                       obs_info.DUR,
                       obs_info.date,
                       obs_info.START, 
                       obs_info.output_path,
                       observing_antennas=obs_info.antennas).generate()

            print('###########################\n')
            print('Running Casa to make the MS')
            print('###########################\n')
            createMeasurementSet(obs_info.flag_antennas,
                                 obs_info.measurementset,
                                 obs_info.casa_path).process()

        else:
            print('\n###########################\n')
            print(('Not running on {}, already exists'
                   '\n').format(obs_info.filename_base))
            print('###########################\n')
