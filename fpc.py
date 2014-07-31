# -*- coding: utf-8 -*-
"""
Created on Mon Dec 02 13:23:48 2013

@author: rjs3
"""
from __future__ import division, print_function
import csv, re, os
import numpy as np
from xml.etree import ElementTree as et
from Tkinter import Tk
from tkFileDialog import askopenfilename,asksaveasfilename
import matplotlib.pyplot as plt

def get_files():
    '''
    Open a dialog and return a set of files to parse. 
    '''
    # we don't want a full GUI, so keep the root window from appearing
    Tk().withdraw() 
    
    # show an "Open" dialog box and return the paths to the selected files
    fullpaths = askopenfilename(multiple=1,defaultextension='.xrdml', 
                    filetypes=[('XRDML','.xrdml'),('CSV','.csv')])
    fullpaths = re.findall('\{(.*?)\}', fullpaths)
    if len(fullpaths):
        print('User opened:', *fullpaths, sep='\n')
    else:
        print('No files selected')
        raise RuntimeError('Exiting, not an error')
    return fullpaths
    
def load_csv(fullpath):
    '''
    Parse Philips CSV files and return the arrays "angle", "cps", and "dcps".
    
    Don't use this yet!
    '''
#    raise NotImplementedError
    # fullpath=r'C:\Users\rjs3\SkyDriveNISTjnc\data\xrr\2014-01-14\fw 100m 160c washed retry crit.csv'
    dialect_dict={'lineterminator': '\r\n', 'skipinitialspace': True, 
             'quoting': 0, 'delimiter': ',', 'quotechar': '"',
             'doublequote': False}
    
    
    csv.register_dialect('philips-csv',dialect_dict)
    ''' check dialect
    with open(fullpath, 'rb') as csvfile:
        sniffdialect = csv.Sniffer().sniff(csvfile.read(),delimiters=',')  
        dialect_dict=vars(sniffdialect)
        didictviews=dialect_dict.viewkeys()
        badviews=[]
        for view in didictviews:
            if '_' in view:
                print('deleting '+view)
                badviews.append(view)
            
        for view in badviews:
                del dialect_dict[view]
             
        csv.register_dialect('philips-csv2',dialect_dict)
        for row in reader:
            print(row)
    '''
    headers={}
    with open(fullpath,'r') as csvfile:
#        csvfile.seek(0)
        reader = csv.reader(csvfile, 'philips-csv')
        for row in reader:
            csvheader=row[0]
            if csvheader=='Sample identification':
                headers['title']=row[1]
            elif csvheader=='K-Alpha1 wavelength':
                headers['wavelength']=row[1]
            elif csvheader=='File date and time':
                headers['date']=row[1]
            elif csvheader=='Time per step':
                time=float(row[1])
            elif csvheader=='Angle':
                break
        #Generate a list of [angle, count] "datapoint" lists
        data=[[float(row[0]),float(row[1])] for row in reader]
        
    # split them into individual lists, then convert to arrays
    angle, counts = tuple(zip(*data))
    angle, counts = np.array(angle), np.array(counts)
    
    cps=counts/time
    dcps=np.sqrt(counts)/time
    keepers = counts>0
    
    return headers, angle[keepers], cps[keepers], dcps[keepers]
        
def load_xrdml(fullpath):
    '''
    Parse Philips XRDML and return the arrays "angle", "cps", and "dcps".
    '''
    headers={}
    xrdtree=et.parse(fullpath)
    xrdroot=xrdtree.getroot()
    
    # I suppose hardcoded index access is "bad", but there would be tons of
    # pointless iteration on hardcoded data otherwise
    headers['title']=xrdroot[0][0].text
    headers['wavelength']=xrdroot[1][1][0].text
    headers['instrument']='X-ray'
    headers['date']=xrdroot[1][5][0][0].text
    
    startangle=float(xrdroot[1][5][1][1][0].text)
    stopangle=float(xrdroot[1][5][1][1][1].text)
    time=float(xrdroot[1][5][1][3].text)
    counts=xrdroot[1][5][1][4].text.split(' ')
    
    counts=np.array([int(count) for count in counts])
    cps=counts/time
    dcps=np.sqrt(counts)/time
    angle=np.linspace(startangle,stopangle,len(cps))
    keepers = counts>0
    
    return headers, angle[keepers], cps[keepers], dcps[keepers]
        
def stitch_data(fullpaths):
    '''
    Load data into list, zip it up, test it, and stitch/sort if necessary
    '''
    data=[load_xrdml(fullpath)
        for fullpath in fullpaths if fullpath.endswith('.xrdml')]
            
    data+=[load_csv(fullpath) 
        for fullpath in fullpaths if fullpath.endswith('.csv')]
    
    # Cut out early if we only load one file
    if len(data)==0:
        print('no supported types selected')
        raise RuntimeError('Exiting, not an error')
    if len(data)==1:
        return data[0]
        
    #Zipped is way more convenient here
    data=zip(*data)

    # Test if our data come from the same sample.
    sentinel=None
    for datum in data[0]:
        if sentinel is None:
            sentinel = datum['title']
        elif sentinel != datum['title']:
            raise ValueError("Data do not appear to be from the same sample."
                             " Try rewriting the sample IDs to match.",
                             sentinel,datum['title'])
                             
    # Stitch away if we get this far 
    headers = data[0][0]
    angle = np.concatenate(data[1])
    cps = np.concatenate(data[2])
    dcps = np.concatenate(data[3])
    
    # test if our data are already sorted
    if np.array_equal(np.sort(angle),angle):
        return headers, angle, cps, dcps
        
    # then sort based on angle
    sortind = np.argsort(angle)
    angle = angle[sortind]
    cps = cps[sortind]
    dcps = dcps[sortind]
    
    return headers, angle, cps, dcps
    
def mymessage(text):
    print(text)
    plt.title(text,fontsize=16)
    plt.draw()
    
def plot_refl_footprint(q,cps):
    '''
    plot, determine axis range from data, then return handle
    '''
    fig = plt.figure()
    plt.clf()
    plt.plot(q,cps,'k-')
    plt.axis([min(q),q[np.argmax(cps)]*2,0,1.1*max(cps)])
    return fig
    
def write_refl(headers, q, refl, drefl, path):
    '''
    Open a file where the previous was located, and drop a refl with the same 
    name as default.
    '''
    # we don't want a full GUI, so keep the root window from appearing
    Tk().withdraw()
    
    fullpath=asksaveasfilename(initialdir=path, 
                initialfile=headers['title']+'.refl',
                defaultextension='.refl', 
                filetypes=[('reflectivity','.refl')])

#    fullpath = re.findall('\{(.*?)\}', s)
    if len(fullpath)==0:
        raise RuntimeError('Exiting, not an error')
        
    textlist = ['#pythonfootprintcorrect 1 1 2014-07-30']
    for header in headers:
        textlist.append('#'+header+' '+headers[header]) 
    textlist.append('#columns Qz R dR')
    
    for point in tuple(zip(q,refl,drefl)):
        textlist.append('{:.12g} {:.12g} {:.12g}'.format(*point))
#    print('\n'.join(textlist))
    with open(fullpath, 'w') as reflfile:
        reflfile.writelines('\n'.join(textlist))
        
    print('Saved successfully as:', fullpath)
    
def footprintCorrect(fullpaths=None,save=True):
    '''
    As a script, ask for filenames, load/stitch data, plot data, 
    collect footprint range, calc footprint, and write .refl files.
    '''
    if fullpaths is None:
        fullpaths = get_files()
    path=os.path.split(fullpaths[0])[0]
    headers,angle,cps,dcps=stitch_data(fullpaths)
    
    q=4*np.pi/float(headers['wavelength'])*np.sin(np.pi/180*angle)
    fig=plot_refl_footprint(q,cps)
    
    mymessage('Click the beginning and end of the footprint region')
    points=plt.ginput(2,timeout=-1)
#    plt.close(fig)
    
    start,stop=zip(*points)[0]
    keepers=q > start
    logical_range = np.logical_and(keepers, q < stop)
    fp_q = q[logical_range]
    fp_cps = cps[logical_range]
    p = np.polyfit(fp_q,fp_cps,1)
    
    q_correct=q[keepers]
    footprint = p[0]*q_correct+p[1]
    plt.plot(q_correct,footprint,'r--')
    plt.draw()
    refl_correct = cps[keepers]/footprint
    refl_correct = refl_correct/max(refl_correct)
    drefl_correct = dcps[keepers]/footprint/max(refl_correct)
    
    plt.figure()
    plt.semilogy(q,cps/max(cps),'k--')
    plt.semilogy(q_correct,refl_correct,'k-')
    mymessage('Resulting reflectivity')
    plt.draw()
    plt.show(block=False)
    
    if save:
        write_refl(headers,q_correct,refl_correct,drefl_correct,path)    
    return q_correct,refl_correct,drefl_correct
           
if __name__ == '__main__':
    footprintCorrect() # We don't need no stinking commandline arguments
    