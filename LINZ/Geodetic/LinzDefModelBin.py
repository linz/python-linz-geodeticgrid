# Imports to support python 3 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

from LINZ.DeformationModel import Model, Time
from collections import namedtuple
from datetime import datetime
import os.path
from .LinzGrid import LinzGrid
import struct
import re

formats={
    'LINZDEF2B':
     { 'signature':  "LINZ deformation model v2.0B\r\n\x1A",
       'gridformat': 'GRID2B',
       'bigendian': True,
       },

    'LINZDEF2L': 
     { 'signature':  "LINZ deformation model v2.0L\r\n\x1A",
       'gridformat': 'GRID2L',
       'bigendian': False,
       },
}

defaultFormat='LINZDEF2L'
defaultStartDate=Time.Time(datetime(1800,1,1))
defaultEndDate=Time.Time(datetime(2200,1,1))

deformation_resolution=0.0001
velocity_resolution=0.000001


class _bbox( object ):

    def __init__( self,ymin=None,ymax=None,xmin=None,xmax=None):
        self.ymin=ymin
        self.ymax=ymax
        self.xmin=xmin
        self.xmax=xmax

    def add( self, other ):
        self.ymin = other.ymin if self.ymin is None or self.ymin > other.ymin else self.ymin
        self.ymax = other.ymax if self.ymax is None or self.ymax < other.ymax else self.ymax
        self.xmin = other.xmin if self.xmin is None or self.xmin > other.xmin else self.xmin
        self.xmax = other.xmax if self.xmax is None or self.xmax < other.xmax else self.xmax

class _packer( object ):

    def __init__( self, bigendian=False ):
        endian=">" if bigendian else "<"
        self.packschar=struct.Struct(endian+'b').pack
        self.packshort=struct.Struct(endian+'h').pack
        self.packlong=struct.Struct(endian+'l').pack
        self.packdouble=struct.Struct(endian+'d').pack

    def writeschar( self, fh, value ): fh.write(self.packschar(value))
    def writeshort( self, fh, value ): fh.write(self.packshort(value))
    def writelong( self, fh, value ): fh.write(self.packlong(value))
    def writedouble( self, fh, value ): fh.write(self.packdouble(value))

    def writestring( self, fh, text ):
        encoded=text.encode('ascii')
        fh.write(self.packshort(len(encoded)+1))
        fh.write(encoded)
        fh.write('\x00')

    def writedate( self, fh, time ):
        date=time.asDateTime()
        for dp in (date.year,date.month,date.day,date.hour,date.minute,date.second):
            fh.write(self.packshort(dp))

    def writebbox( self, binfile, bbox ):
        if bbox.ymin==None or bbox.ymax==None or bbox.xmin==None or bbox.xmax==None:
            raise RuntimeError('Cannot write uninitialized bbox')
        self.writedouble(binfile,bbox.ymin)
        self.writedouble(binfile,bbox.ymax)
        self.writedouble(binfile,bbox.xmin)
        self.writedouble(binfile,bbox.xmax)

class LinzDefModelBin( object ):
    '''
    The LinzGrid class is used for reading/writing the binary grid format used
    by LINZ software such as SNAP and Landonline.
    '''

    def __init__( self, model, format='LINZDEF2L', description=None, verbose=False):
        '''
        Create a LinzDefModelBin object.  

        Required parameters are:

        model: the LINZ.DeformationModel.Model.Model to be uploaded
        
        Optional parameters are:

        format: one of the format definition strings used by the grid format (eg LINZDEF2L, LINZDEF2B)
        '''
        if format not in formats:
            raise RuntimeError('Invalid grid format {0}'.format(format))
        
        self.format=format
        self.formatdef=formats[format]
        self.model=model
        self.verbose=verbose
        self.analyzeModel()

    def analyzeModel( self ):
        model=self.model
        refdate=Time.Time(model.metadata('datum_epoch'))
        self.datum_epoch=refdate
        self.datum_code=model.metadata('datum_code')
        

        TimeStep=namedtuple('TimeStep','mtype t0 f0 t1 f1')
        DefSeq=namedtuple('DefSeq','component description dimension zerobeyond steps grids timefuncs bbox')
        DefComp=namedtuple('DefComp','date factor before after')
        SeqComp=namedtuple('SeqComp','time factor before after nested')
        TimeFunc=namedtuple('TimeFunc','type params')
        GridDef=namedtuple('GridDef','name description dimension isvelocity bbox')
        small=0.00001

        class TimeEvent:
            def __init__(self,time,f0,time1,f1):
                self.time0=time
                self.time1=time1
                self.days=time.daysAfter(refdate)
                self.f0=f0
                self.f1=f1

            def __str__(self):
                return 'Time:{0} f0:{1} {2} f1:{3}'.format(self.time0,self.f0,self.time1,self.f1)

            def __cmp__( self, other ):
                return cmp( (self.time0,self.time1), (other.time0,other.time1) )

            def __repr__(self):
                return str(self)

        sequences=[]
        gridfiles={}

        for c in model.components():
            if self.verbose: print("Analyzing component:",c.name)
            component=c.submodel
            tm=c.timeFunction
            mtype=tm.time_function
            if mtype not in ('velocity','step','ramp'):
                raise RuntimeError('Cannot handle temporal model type '+mtype)

            grids=[]
            dimension=0
            zerobeyond=True
            for m in c.spatialModel.models():
                if m.spatial_model != 'llgrid':
                    raise RuntimeError('Cannot handle spatial model type '+spatial_model)
                gridfile=m.model().gridFile()
                grids.append(gridfile)
                dimension=len(m.columns)
                if not m.spatial_complete:
                    zerobeyond=False
                if gridfile not in gridfiles:
                    gridfiles[gridfile]=GridDef(gridfile,getattr(m,'description',''),dimension,mtype == 'velocity',_bbox())

            # Reverse grids so that contained grids occur before containing grids..
            grids.reverse()
            step=TimeStep(mtype,tm.time0,tm.factor0,tm.time1,tm.factor1)

            found = False
            for s in sequences:
                if (s.component == component and
                    s.zerobeyond == zerobeyond and
                    s.grids == grids ):
                    s.steps.append(step)
                    found=True
                    break

            if not found:
                sequences.append(DefSeq(component,c.description,dimension,zerobeyond,[step],grids,[],_bbox()))

        for sequence in sequences:
            compname=sequence.component
            if self.verbose: print("Analyzing sequence:",compname)

            timefuncs = []
            events=[]
            for s in sequence.steps:
                if s.mtype == 'velocity':
                    timefuncs.append(TimeFunc('VELOCITY',[s.t0]))
                elif s.mtype == 'step':
                    events.append(TimeEvent(s.t0,s.f0,s.t0,s.f1))
                elif s.mtype == 'ramp':
                    events.append(TimeEvent(s.t0,s.f0,s.t1,s.f1))
                else:
                    raise RuntimeError('Cannot handle time model type '+s.mtype)

            pwm=''
            if len(events) > 0:
                time_model=[]
                events.sort()
                e0=None
                for e in events:
                    if e0 and e0.time1.daysAfter(e.time0) > 0.001:
                        raise RuntimeError('Cannot handle overlapping time events in series')
                    v0 = 0.0 if len(time_model) == 0 else time_model[-1][1]
                    for t in time_model:
                        t[1] += e.f0
                    time_model.append([e.time0,e.f0+v0])
                    time_model.append([e.time1,e.f1+v0])
                for i in reversed(range(len(time_model)-1)):
                     t0=time_model[i]
                     t1=time_model[i+1]
                     if abs(t0[0].daysAfter(t1[0])) < 0.001 and abs(t0[1]-t1[1]) < 0.00001:
                         time_model[i:i+1]=[]
                steps=[time_model[0][1]]
                i0=1 if time_model[1][0].daysAfter(time_model[0][0]) < 0.001 else 0
                for t in time_model[i0:]:
                    steps.append(t[0])
                    steps.append(t[1])

                timefuncs.append(TimeFunc('PIECEWISE_LINEAR',steps))
            sequence.timefuncs[:]=timefuncs

        self.sequences=sequences
        self.gridfiles=gridfiles

    def write(self,binfile):
        '''
        Write the model to a file handle. 
        '''
        packer=_packer(self.formatdef['bigendian'])
        model=self.model
        binfile.write(self.formatdef['signature'].encode('ascii'))
        # Create a pointer to the file index data, which is written immediately 
        # after the pointer in this case (unlike previous perl code)
        indexptrloc=binfile.tell()
        packer.writelong(binfile,0)
        # Write each of the grids used and record its location
        gridfiles=self.gridfiles
        gridloc={}
        bbox=_bbox()
        for g in sorted(gridfiles):
            gridloc[g]=binfile.tell()
            if self.verbose: print("Writing grid {0} at {1}".format(g,gridloc[g]))
            gdef=gridfiles[g]
            gf=LinzGrid(
                format=self.formatdef['gridformat'],
                coordsys=self.datum_code,
                description=[g,gdef.description],
                csvfile=self.model.getFileName(g),
                crdcols=['lon','lat'],
                datacols=['du'] if gdef.dimension==1 else ['de','dn'] if gdef.dimension==2 else ['de','dn','du'],
                resolution=velocity_resolution if gdef.isvelocity else deformation_resolution
                )
            gf.write(binfile)
            gridbbox=_bbox(gf.ymin,gf.ymax,gf.xmin,gf.xmax)
            gdef.bbox.add(gridbbox)
            bbox.add(gridbbox)
        # compile sequence bbox and total bbox
        for s in self.sequences:
            for g in s.grids:
                s.bbox.add(gridfiles[g].bbox)

        indexloc=binfile.tell()
        packer.writestring(binfile,model.metadata('model_name'))
        packer.writestring(binfile,model.version())
        packer.writestring(binfile,self.datum_code)
        packer.writestring(binfile,model.metadata('description'))
        packer.writedate(binfile,model.versionInfo(model.version()).release_date)
        packer.writedate(binfile,defaultStartDate)
        packer.writedate(binfile,defaultEndDate)
        packer.writebbox(binfile,bbox)
        # Coords are lat/lon flag - always true for LINZ deformation model
        packer.writeshort(binfile,1)

        defseq=[]
        for sequence in self.sequences:
            name=sequence.component
            if len(sequence.timefuncs):
                name=name+'_f{0}'
            for i,timefunc in enumerate(sequence.timefuncs):
                defseq.append((name.format(i),sequence,timefunc))

        packer.writeshort(binfile,len(defseq))

        for name,sequence,timefunc in defseq:
            packer.writestring(binfile,name)
            packer.writestring(binfile,sequence.description)
            packer.writedate(binfile,defaultStartDate)
            packer.writedate(binfile,defaultEndDate)
            packer.writebbox(binfile,sequence.bbox)
            packer.writeshort(binfile,sequence.dimension)
            packer.writeshort(binfile,1 if sequence.zerobeyond else 0)
            # Nested sequence - always true for implementation
            packer.writeshort(binfile,1)
            packer.writeshort(binfile,len(sequence.grids))

            for gridfile in sorted(sequence.grids):
                gridname=name+'_'+os.path.basename(gridfile)
                if timefunc.type == 'VELOCITY':
                    tref=timefunc.params[0].asYear()
                    t0=defaultStartDate.asYear()
                    t1=defaultEndDate.asYear()
                    timemodel=[t0-tref,defaultStartDate,t0-tref,defaultEndDate,t1-tref]
                elif timefunc.type == 'PIECEWISE_LINEAR':
                    timemodel=timefunc.params
                else:
                    raise RuntimeError('Invalid timefunc type {0}'.format(timefunc.type))
                packer.writestring(binfile,gridname)
                packer.writedate(binfile,timemodel[1])
                packer.writebbox(binfile,gridfiles[gridfile].bbox)
                # Time model type - always piecewise linear
                packer.writeshort(binfile,1)
                nstep=int((len(timemodel)-1)/2)
                packer.writeshort(binfile,nstep)
                packer.writedouble(binfile,timemodel[0])
                for ns in range(nstep):
                    packer.writedate(binfile,timemodel[ns*2+1])
                    packer.writedouble(binfile,timemodel[ns*2+2])
                # Spatial model always grid
                packer.writeshort(binfile,0)
                packer.writelong(binfile,gridloc[gridfile])
        endloc=binfile.tell()
        binfile.seek(indexptrloc)
        packer.writelong(binfile,indexloc)
        binfile.seek(endloc)

            
    def writefile( self, filename ):
        '''
        Write the grid to a file
        '''
        with open(filename,'wb') as binfile:
            self.write( binfile )
        if self.verbose:
            print("Deformation model written to {0}".format(filename))

    @staticmethod
    def main():
        import argparse
        formatcodes=sorted(formats.keys())
        argparser=argparse.ArgumentParser('Build a LINZDEF deformation file from a published deformation model')
        argparser.add_argument('model_dir',help='Base directory of model')
        argparser.add_argument('linzdef_file',help='Name of final model')
        argparser.add_argument('-f','--format',choices=formatcodes,default=defaultFormat,
                               help='Format of binary file')
        argparser.add_argument('-v','--verbose',action='storetrue',help='Generate more output')

        args = argparser.parse_args()

        model_dir=args.model_dir
        if not isdir(model_dir):
            raise RuntimeError("Invalid model directory: "+md)

        model=Model.Model(model_dir)
        binfiledef=LinzDefModelBin(model,format=args.format,verbose=args.verbose)
        binfiledef.writefile(args.linzdef_file)
        
if __name__ == "__main__":
    LinzDefModelBin.main()
