
from LINZ.DeformationModel import Model, Time
from collections import namedtuple
from datetime import datetime
import LinzGrid
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
defaultStartDate=datetime(1800,1,1)
defaultEndDate=datetime(2200,1,1)

deformation_resolution=0.0001
velocity_resolution=0.000001

class _packer( object ):

    def __init__( self, bigendian=False ):
        endian=">" if bigendian else "<"
        self.packschar=struct.Struct(endian+'b').pack
        self.packshort=struct.Struct(endian+'h').pack
        self.packlong=struct.Struct(endian+'l').pack
        self.packdouble=struct.Struct(endian+'d').pack

    def writestring( self, fh, text ):
        encoded=text.encode('ascii')
        fh.write(self.packshort(len(encoded)+1))
        fh.write(encoded)
        fh.write('\x00')

    def writedate( self, fh, date ):
        for dp in (date.year,date.month,date.day,date.hour,date.minute,date.second):
            fh.write(self.packshort(dp))

class _range( object ):

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


    def write( self, packer, binfile ):
        if self.ymin==None or self.ymax==None or self.xmin==None or self.xmax==None:
            raise RuntimeError('Cannot write uninitialized range')
        binfile.write(packer.packdouble(self.ymin))
        binfile.write(packer.packdouble(self.ymax))
        binfile.write(packer.packdouble(self.xmin))
        binfile.write(packer.packdouble(self.xmax))

class LinzDefModelBin( object ):
    '''
    The LinzGrid class is used for reading/writing the binary grid format used
    by LINZ software such as SNAP and Landonline.
    '''

    def __init__( self, model, format='LINZDEF2L', description=None, coordsys='',islatlon=True,resolution=None,csvfile=None):
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
        self.analyzeModel()

    def analyzeModel( self ):
        model=self.model
        refdate=Time.Time(model.metadata('datum_epoch'))
        self.datum_epoch=refdate
        self.datum_code=model.metadata('datum_code')
        

        TimeStep=namedtuple('TimeStep','mtype t0 f0 t1 f1')
        DefSeq=namedtuple('DefSeq','component dimension zerobeyond steps grids subseq range')
        DefComp=namedtuple('DefComp','date factor before after')
        SeqComp=namedtuple('SeqComp','time factor before after nested')
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
            print "Analyzing component:",c.name
            component=c.submodel
            # print component
            tm=c.timeFunction
            mtype=tm.time_function
            if mtype not in ('velocity','step','ramp'):
                raise RuntimeError('Cannot handle temporal model type '+mtype)
            # print type(c.timeComponent.model())
            # print mtype,tm.factor0,tm.factor1,tm.time0,tm.time1

            grids=[]
            dimension=0
            zerobeyond=True
            for m in c.spatialModel.models():
                if m.spatial_model != 'llgrid':
                    raise RuntimeError('Cannot handle spatial model type '+spatial_model)
                # print m.model().gridSpec()
                #print '    ',m.model().gridFile()
                #print '    ',m.columns
                gridfile=m.model().gridFile()
                grids.append(gridfile)
                dimension=len(m.columns)
                if not m.spatial_complete:
                    zerobeyond=False
                if gridfile not in gridfiles:
                    gridfiles[gridfile]={
                        'name': gridfile,
                        'floc': 0,
                        'dimension': dimension,
                        'description': getattr(m,'description',''),
                        'isvelocity': mtype == 'velocity',
                        'range': _range(),
                        }

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
                sequences.append(DefSeq(component,dimension,zerobeyond,[step],grids,[],_range()))

        for sequence in sequences:
            compname=sequence.component
            print "Analyzing sequence:",compname

            subsequences = []
            events=[]
            for s in sequence.steps:
                if s.mtype == 'velocity':
                    subsequences.append([s.t0,'VELOCITY'])
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
                pwm='PIECEWISE_LINEAR {0}'.format(time_model[0][1])
                i0=1 if time_model[1][0].daysAfter(time_model[0][0]) < 0.001 else 0
                for t in time_model[i0:]:
                    pwm = pwm+" {0} {1}".format( t[0].strftime('%d-%b-%Y'),t[1] )

                subsequences.append([refdate,pwm])
            sequence.subseq[:]=subsequences

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
        binfile.write(packer.packlong(0))
        # Write each of the grids used and record its location
        gridfiles=self.gridfiles
        range=_range()
        for g in sorted(gridfiles):
            print("Writing grid "+g)
            gridfiles[g]['floc']=binfile.tell()
            gf=LinzGrid.LinzGrid(
                format=self.formatdef['gridformat'],
                coordsys=self.datum_code,
                description=[g,gridfiles[g]['description']],
                csvfile=self.model.getFileName(g),
                resolution=velocity_resolution if gridfiles[g]['isvelocity'] else deformation_resolution
                )
            gf.write(binfile)
            gridrange=_range(gf.ymin,gf.ymax,gf.xmin,gf.xmax)
            gridfiles[g]['range']=gridrange
            range.add(gridrange)
        # compile sequence ranges and total range
        for s in self.sequences:
            for g in s.grids:
                s.range.add(gridfiles[g]['range'])

        indexloc=binfile.tell()
        packer.writestring(binfile,model.metadata('model_name'))
        packer.writestring(binfile,model.version())
        packer.writestring(binfile,self.datum_code)
        packer.writestring(binfile,model.metadata('description'))
        packer.writedate(binfile,model.versionInfo(model.version()).release_date)
        packer.writedate(binfile,defaultStartDate)
        packer.writedate(binfile,defaultEndDate)
        range.write(packer,binfile)
        # Coords are lat/lon flag - always true for LINZ deformation model
        binfile.write(packer.packshort(1))

        defseq=[]
        for sequence in self.sequences:
            name=sequence.component.name
            if len(sequence.subsequences):
                name=name+'_{0}'
            for i,subsequence in enumerate(sequence.subsequences):
                defseq.append(name.format(i),sequence,subsequence):

        binfile.write(packer.packshort(len(defseq)))

        for name,sequence,subsequence in enumerate(self.sequences):
            print(name)
            print(sequence.component.name)
            print(sequence.component.description)
            packer.writestring(binfile,sequence.component.name)
            packer.writestring(binfile,sequence.component.description)
            packer.writedate(binfile,defaultStartDate)
            packer.writedate(binfile,defaultEndDate)
            sequence.range.write(packer,binfile)
            binfile.write(packer.packshort(sequence.dimension))
            binfile.write(packer.packshort(1 if sequence.zerobeyondrange else 0))
            # Nested sequence - always true for implementation
            binfile.write(packer.packshort(1))
            binfile.write(packer.packshort(len(sequence.grids)))
            for gridfile in sequence.grids:
                gridname=name+os.path.basename(gridfile)


            
    def writefile( self, filename ):
        '''
        Write the grid to a file
        '''
        with open(filename,'wb') as binfile:
            self.write( binfile )

    @staticmethod
    def main():
        import argparse
        formatcodes=sorted(formats.keys())
        argparser=argparse.ArgumentParser('Build a LINZDEF deformation file from a published deformation model')
        argparser.add_argument('model_dir',help='Base directory of model, containing model and tools directories')
        argparser.add_argument('linzdef_file',help='Name of final model')
        argparser.add_argument('-f','--format',choices=formatcodes,default=defaultFormat,
                               help='Format of binary file')

        args = argparser.parse_args()

        md=args.model_dir
        if not isdir(md) or not isdir(joinpath(md,'model')) or not isdir(joinpath(md,'tools')):
            raise RuntimeError("Invalid model directory: "+md)

        bd=args.target_dir
        if not isdir(bd):
            if os.path.exists(bd):
                raise RuntimeError("Invalid target build directory: "+bd)
            else:
                os.makedirs(bd)

        defname=args.linzdef_file
        deffile=open(joinpath(bd,defname+'.def'),"w")

        sys.path.append(joinpath(md,'tools'))
        from LINZ.DeformationModel import Model, Time

        m=Model.Model(joinpath(md,'model'))
        
if __name__ == "__main__":
    # LinzDefModelBin.main()
    model=Model.Model('/home/ccrook/projects/deformation/models/published/model')
    binfiledef=LinzDefModelBin(model)
    binfiledef.writefile('garbage_model.bin')
