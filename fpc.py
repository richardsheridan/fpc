# -*- coding: utf-8 -*-
"""
Created on Mon Dec 02 13:23:48 2013

@author: rjs3
"""
from __future__ import division, print_function
import csv, os, sys
ver = sys.version_info.major
import numpy as np
import matplotlib.pyplot as plt
from xml.etree import ElementTree as et
if ver >= 3:
    from tkinter import Tk
    from tkinter.filedialog import askopenfilename,asksaveasfilename
    from tkinter.simpledialog import askfloat
else:
    from Tkinter import Tk
    from tkFileDialog import askopenfilename,asksaveasfilename
    from tkSimpleDialog import askfloat

def get_files():
    '''
    Open a dialog and return a set of files to parse.
    '''
    # we don't want a full GUI, so keep the root window from appearing
    Tk().withdraw()

    # show an "Open" dialog box and return the paths to the selected files
    fullpaths = askopenfilename(multiple=1,#defaultextension='.xrdml',
                    filetypes=[('XRDML','.xrdml'),('CSV','.csv'),('All files','.*')])

    if len(fullpaths):
        print('User opened:', *fullpaths, sep='\n')
    else:
        print('No files selected')
        raise ExitException('Exiting, not an error')
    return fullpaths

def crop_noise(counts):
    keepers = np.zeros_like(counts,dtype=bool)
    for i,count in enumerate(counts):
        if count:
            keepers[i] = True
        else:
            break

    if not any(keepers):
        raise RuntimeError('No nonzero data in this file')

    return keepers

def load_csv(fullpath):
    '''
    Parse Philips CSV files and return the arrays "angle", "cps", and "dcps".

    '''
    headers={}
    with open(fullpath,'r') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            csvheader=row[0]
            if csvheader=='Sample identification':
                headers['title']=row[1].strip()
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
    keepers = crop_noise(counts)

    return headers, angle[keepers], cps[keepers], dcps[keepers]

def load_xrdml(fullpath):
    '''
    Parse Philips XRDML and return the arrays "angle", "cps", and "dcps".
    '''
    headers={}
    xrdtree=et.parse(fullpath)
    ns='/{http://www.xrdml.com/XRDMeasurement/1.0}'

    headers['title']=xrdtree.findtext(ns.join(('.','sample','id'))).strip()
    headers['wavelength']=xrdtree.findtext(ns.join(('.','xrdMeasurement',
        'usedWavelength','kAlpha1')))
    headers['instrument']='X-ray'
    headers['date']=xrdtree.findtext(ns.join(('.','xrdMeasurement','scan',
        'header','startTimeStamp')))

    dataPoints=xrdtree.find(ns.join(('.','xrdMeasurement','scan','dataPoints')))

    for child in dataPoints:
        tag=child.tag
        if tag == ns[1:]+'positions' and child.get('axis')=='Omega':
            startangle=float(child.findtext('.'+ns+'startPosition'))
            stopangle=float(child.findtext('.'+ns+'endPosition'))
        elif tag == ns[1:]+'commonCountingTime':
            time=float(child.text)
        elif tag == ns[1:]+'intensities':
            counts=child.text

    counts=np.fromstring(counts,sep=' ',dtype=int)
    cps=counts/time
    dcps=np.sqrt(counts)/time
    angle=np.linspace(startangle,stopangle,len(cps))
    keepers = crop_noise(counts)

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
        raise ExitException('Exiting, not an error')
    if len(data)==1:
        return data[0]

    # Zipped is way more convenient here
    data=list(zip(*data))

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
    return headers, angle[sortind], cps[sortind], dcps[sortind]

def mymessage(text):
    print(text)
    plt.title(text,fontsize=16)
    plt.draw()

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
        print('Results not saved')
        raise ExitException('Exiting, not an error')

    textlist = ['#pythonfootprintcorrect 1 1 2014-09-11']
    for header in headers:
        textlist.append('#'+header+' '+headers[header])
    textlist.append('#columns Qz R dR')

    for point in tuple(zip(q,refl,drefl)):
        textlist.append('{:.12g} {:.12g} {:.12g}'.format(*point))
#    print('\n'.join(textlist))
    with open(fullpath, 'w') as reflfile:
        reflfile.writelines('\n'.join(textlist))

    print('Saved successfully as:', fullpath)

def q_from_angle(angle,wavelength=1.5405980):
    return 4*np.pi/wavelength*np.sin(np.pi/180*angle)

def footprintCorrect(fullpaths=None):
    '''
    As a script, ask for filenames, load/stitch data, plot data,
    collect footprint range, calc footprint, and write .refl files.
    '''
    if fullpaths is None:
        fullpaths = get_files()
    path=os.path.split(fullpaths[0])[0]
    headers,angle,cps,dcps=stitch_data(fullpaths)

    wavelength=float(headers['wavelength'])
    q=q_from_angle(angle,wavelength)

    plt.plot(q,cps,'k-')
    plt.axis([min(q),q[np.argmax(cps)]*2,0,1.1*max(cps)])

    mymessage('Click the beginning and end of the footprint region')
    points=plt.ginput(2,timeout=100)
    if len(points)<2:
        raise ExitException('Point selection timeout, not an error')

#    cutoff=askfloat('Cutoff angle','A footprint cutoff angle in degrees',
#             initialvalue=1.0)
    Tk().withdraw()
    width=askfloat('Sample width','Average sample width in mm',initialvalue=24.5)
    if width is None:
        raise ExitException('Exiting, not an error')
    cutoff=np.arcsin(.4/width)/np.pi*180
    q_cut=q_from_angle(cutoff,wavelength)

    start,stop = points[0][0], points[1][0]
    keepers = q > start
    fit_range = np.logical_and(keepers, q < stop)
    p = np.polyfit(q[fit_range],cps[fit_range],1)

    q_correct = q[keepers]
    footprint = np.polyval(p, q_correct)
    fp_region = q_correct<q_cut
    footprint[np.logical_not(fp_region)]=np.polyval(p,q_cut)

    plt.plot(q_correct,footprint,'r--')
    plt.draw()

    refl_correct = cps[keepers]/footprint
    refl_correct = refl_correct/max(refl_correct)
    drefl_correct = dcps[keepers]/footprint/max(refl_correct)

    plt.figure()
    plt.semilogy(q,cps/max(cps),'k--')
    plt.semilogy(q_correct,footprint/max(cps),'r:')
    plt.semilogy(q_correct,refl_correct,'k-')
    mymessage('Resulting reflectivity')
    plt.draw()

    plt.show(block=False)

    write_refl(headers,q_correct,refl_correct,drefl_correct,path)
    return q_correct,refl_correct,drefl_correct

class ExitException(Exception):
     def __init__(self, value):
         self.value = value
     def __str__(self):
         return repr(self.value)

if __name__ == '__main__':
    try:
        footprintCorrect() # We don't need no stinking commandline arguments
    except ExitException:
        pass