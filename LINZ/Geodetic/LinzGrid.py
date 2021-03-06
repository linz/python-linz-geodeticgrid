import numpy as np
import csv
import struct
import re

formats = {
    "GEOID": "SNAP geoid binary file\r\n\x1A",
    "GRID1L": "SNAP grid binary v1.0 \r\n\x1A",
    "GRID1B": "CRS grid binary v1.0  \r\n\x1A",
    "GRID2L": "SNAP grid binary v2.0 \r\n\x1A",
    "GRID2B": "CRS grid binary v2.0  \r\n\x1A",
}

bigendian_formats = ("GRID1B", "GRID2B")
version2_formats = ("GRID2L", "GRID2B")

mincompress = 16  # Minimum length of array to pack with differences
maxsubset = 1024  # Maximum number of elements in a subset

# Useful constants
max1byte = 2 ** 7 - 2
max2byte = 2 ** 15 - 2
max4byte = 2 ** 31 - 2
missing = 2 ** 31 - 1


class _packer(object):
    def __init__(self, bigendian=False):
        endian = ">" if bigendian else "<"
        self.packschar = struct.Struct(endian + "b").pack
        self.packshort = struct.Struct(endian + "h").pack
        self.packlong = struct.Struct(endian + "l").pack
        self.packdouble = struct.Struct(endian + "d").pack

    def _writestring(self, text, fh):
        encoded = text.encode("ascii")
        fh.write(self.packshort(len(encoded) + 1))
        fh.write(encoded)
        fh.write("\x00")


class LinzGrid(object):
    """
    The LinzGrid class is used for reading/writing the binary grid format used
    by LINZ software such as SNAP and Landonline.
    """

    def __init__(
        self,
        format="GRID2L",
        description=None,
        coordsys="",
        islatlon=True,
        resolution=None,
        csvfile=None,
        crdcols=None,
        datacols=None,
        separator=None,
        colnames=None,
    ):
        """
        Create a LinzGrid object.  Optional parameters are:

        format: one of the format definition strings used by the grid format (eg GRID2L)
        description: header text describing the grid (see setHeaders)
        coordsys: code identifying the coordinate system of the lat/lon grid
        islatlon: identifies that the grid is latitude and longitude (ie not projection coords)
        resolution: precision of the grid data values (values are integer multiples of resolution)
        csvfile: optional name of a csvfile from which grid data is read
        crdcols: optional name of east and north coordinate columns (default [lon, lat])
        datacols: optional names of data columns, default is not coord columns
        """
        if format not in formats:
            raise RuntimeError("Invalid grid format {0}".format(format))

        self.format = format
        self.setDescription(description or "")
        self.coordsys = coordsys
        self.islatlon = islatlon
        self.ngridx = 0
        self.ngridy = 0
        self.xmin = 0
        self.xmax = 0
        self.ymin = 0
        self.ymax = 0
        self.ndim = 0
        self.resolution = resolution
        self.data = None
        if csvfile:
            self.loadFile(
                csvfile, resolution=resolution, crdcols=crdcols, datacols=datacols, separator=separator, colnames=colnames
            )

    def setDescription(self, *description):
        """
        Set the descriptive headers for the model.  The model may have
        up to three headers, each of which is a string. The final string
        may include new line characters and is the detailed description. 
        Headers may be entered as a single string, a list of strings, 
        or a set of string parameters.
        """
        lines = []
        description = list(description)
        try:
            while description:
                h = description.pop(0)
                if isinstance(h, str):
                    lines.extend(h.split("\n"))
                else:
                    description[0:0] = h
        except:
            raise
            raise RuntimeError("Invalid description")
        while len(lines) < 3:
            lines.append("")
        if len(lines) > 3:
            lines[2:] = ["\n".join(lines[2:])]
        self.description = lines

    def _csvIterator(self, csvfile):
        with open(csvfile) as csvfh:
            csvr = csv.reader(csvfh)
            for r in csvr:
                yield r

    def _spaceIterator(self, csvfile):
        with open(csvfile) as csvfh:
            for l in csvfh:
                yield l.split()

    def _fileIterator(self, csvfile, crdcols, datacols, separator="csv", colnames=None):
        iterator = None
        if separator == "csv":
            iterator = self._csvIterator(csvfile)
        elif separator == "space":
            iterator = self._spaceIterator(csvfile)
        else:
            raise RuntimeError("Invalid grid file separator {0}".format(separator))
        nrow = 0
        if colnames is None:
            columns = next(iterator)
            nrow += 1
        else:
            columns = re.findall(r"\w+", colnames)
        columns = [c.lower() for c in columns]
        if len(columns) < 3:
            raise RuntimeError("Grid file {0} must have at least 3 columns".format(csvfile))
        if crdcols is None:
            if "longitude" in columns and "latitude" in columns:
                crdcols = ("longitude", "latitude")
            else:
                crdcols = ("lon", "lat")
        elif len(crdcols) == 2:
            crdcols = tuple(c.lower() for c in crdcols[:2])
        else:
            raise RuntimeError("Grid file {0} called with invalid coordinate columns".format(csvfile))
        crdcolno = tuple(columns.index(c) if c in columns else -1 for c in crdcols)
        if -1 in crdcolno:
            raise RuntimeError("Grid file {0} must have columns {1} and {2}".format(csvfile, crdcols[0], crdcols[1]))
        datacolno = []
        if datacols is None:
            datacolno = tuple(i for i in range(len(columns)) if i not in crdcolno)
        else:
            datacols = [d.lower() for d in datacols]
            for c in datacols:
                if c in columns:
                    datacolno.append(columns.index(c))
                else:
                    raise RuntimeError("Grid file {0} must have data columns {1}".format(csvfile, c))

        ncolmin = crdcolno[0] if crdcolno[0] < crdcolno[1] else crdcolno[1]
        for cn in datacolno:
            if cn > ncolmin:
                ncolmin = cn
        notblank = re.compile(r"\S")

        for row in iterator:
            nrow += 1
            if len(row) <= ncolmin:
                continue
            try:
                ll = tuple(float(row[x]) for x in crdcolno)
                data = tuple(float(row[x]) if notblank.match(row[x]) else None for x in datacolno)
            except:
                raise RuntimeError("Error reading grid file {0} at line {1}".format(csvfile, nrow))
            yield ll, data

    def loadFile(self, csvfile, resolution=None, crdcols=None, datacols=None, separator="csv", colnames=None):
        """
        Load grid data from a csv file.  The csv file is assumed to have as 
        single header line holding the column names.  East and north coordinates
        are read from the columns named in crdcols parameter, and data from the columns listed
        in datacols parameter. If crdcols is None then either ('lon','lat'), or 
        ('longitude','latitude') will be used.  If datacols is None then all none
        coordinate columns will be used.
        The grid must be a regular longitude/latitude grid (ie aligned with 
        the lon/lat axes and with equally spaced columns/rows).  The order 
        of elements must be entered sequentially along rows then columns or vice versa,
        No assumptions are made about the direction that rows and columns are entered,
        or the order.
        Data values are converted to integers by dividing by resolution.
        """

        ndata = 0
        ncol = 0
        ll0 = None
        ll1 = None
        ll2 = None
        dmin = None
        dmax = None
        reader = lambda: self._fileIterator(csvfile, crdcols, datacols, separator, colnames)
        for ll, data in reader():
            ndim = len(data)
            ndata += 1
            if ndata == 1:
                ll0 = ll
            elif ndata == 2:
                dll0 = (ll[0] - ll0[0], ll[1] - ll0[1])
            elif ll1 is None:
                dll1 = (ll[0] - ll2[0], ll[1] - ll2[1])
                if dll1[0] * dll0[0] + dll1[1] * dll0[1] < 0:
                    ll1 = ll2
                    ncol = ndata - 1
            ll2 = ll
            for f in data:
                if f is not None:
                    if dmin is None:
                        dmin = dmax = f
                    elif f < dmin:
                        dmin = f
                    elif f > dmax:
                        dmax = f

        if ncol == 0 or ndata % ncol != 0:
            raise RuntimeError(
                "Grid file {0} number of elements {1} not multiple of row size {2}".format(csvfile, ndata, ncol)
            )
        nrow = int(ndata / ncol)
        tolerance = 0.000001
        gdata = np.full((nrow, ncol, ndim), missing, dtype=np.int32)
        lon0 = ll0[0]
        lat0 = ll0[1]
        dlonc = ll1[0] - lon0
        dlatc = ll1[1] - lat0
        dlonr = ll2[0] - ll1[0]
        dlatr = ll2[1] - ll1[1]
        if abs(dlonc) > tolerance and abs(dlatc) > tolerance:
            raise RuntimeError("Grid file {0} - grid not aligned with lat/lon axes".format(csvfile))
        if abs(dlonr) > tolerance and abs(dlatr) > tolerance:
            raise RuntimeError("Grid file {0} - grid not aligned with lat/lon axes".format(csvfile))

        dlonc /= ncol - 1
        dlatc /= ncol - 1
        dlonr /= nrow - 1
        dlatr /= nrow - 1

        resolution = resolution if resolution is not None else self.resolution
        if resolution is None or resolution == 0.0:
            resolution = 1.0
            vmax = abs(dmin) if abs(dmin) > abs(dmax) else abs(dmax)
            # Limits defining default resolution
            while resolution > 1.0e-9 and vmax < 1.0e7:
                resolution /= 10.0
                vmax *= 10.0

        irow = 0
        icol = -1
        ndata = 1
        for ll, data in reader():
            ndata += 1
            icol += 1
            if icol >= ncol:
                icol = 0
                irow += 1
            lonerr = abs(ll[0] - ll0[0] - dlonc * icol - dlonr * irow)
            laterr = abs(ll[1] - ll0[1] - dlatc * icol - dlatr * irow)
            if lonerr > tolerance or laterr > tolerance:
                raise RuntimeError("Grid file {0} - grid not regular grid at record {1}".format(csvfile, ndata))
            for i, f in enumerate(data):
                if f is not None:
                    gdata[irow, icol, i] = int(round(f / resolution))

        if abs(dlatc) > abs(dlonc):
            gdata = gdata.transpose(1, 0, 2)
        if ll0[0] > ll2[0]:
            gdata = gdata[:, ::-1, :]
        if ll0[1] > ll2[1]:
            gdata = gdata[::-1, :, :]

        self.data = gdata
        self.ngridx = gdata.shape[1]
        self.ngridy = gdata.shape[0]
        self.xmin = ll0[0] if ll0[0] < ll2[0] else ll2[0]
        self.xmax = ll2[0] if ll0[0] < ll2[0] else ll0[0]
        self.ymin = ll0[1] if ll0[1] < ll2[1] else ll2[1]
        self.ymax = ll2[1] if ll0[1] < ll2[1] else ll0[1]
        self.ndim = ndim
        self.resolution = resolution

    def _writestring(self, text, gridfile, packer):
        encoded = text.encode("ascii")
        gridfile.write(packer.packshort(len(encoded) + 1))
        gridfile.write(encoded)
        gridfile.write("\x00".encode("ascii"))

    def _writerowv1(self, irow, gridfile, packer):
        rowloc = gridfile.tell()
        data = self.data[irow]
        if np.amax(np.absolute(data[data != missing])) > max2byte:
            raise RuntimeError("Data cannot be stored at required resolution in this format")
        for r in data:
            for v in r:
                if v == missing:
                    gridfile.write(packer.packshort(0x7FFF))
                else:
                    gridfile.write(packer.packshort(v))
        return rowloc

    def _writerowv2a(self, irow, idim, gridfile, packer):
        # The data for each row is written as a series of subsets.  Each is
        # defined by a short value identifying format with bits
        #    0      subset
        #    1      cont
        #    2,3    difference level (0,1,2)
        #    4,..   bytes per value (1,2,4)
        #
        # If subset is set then this is followed by two short values defining
        # the start and end of the subset.
        # If cont is non-zero then the data is followed by another subset
        # Differencing may be used to compress the data, either single or double differencing
        # Data is compressed to 1, 2, or 4 bytes per value depending on the data variation
        # For single and double differences the data start with 1 or 2 4 byte values followed
        # by the differences using the specified number of bytes
        # Missing values are represented by 2**(nbytes*8)-1

        # Split out subsets between missing data
        rowdata = self.data[irow, :, idim]
        usecols = np.where(rowdata != missing)[0]
        subset = usecols * 0
        usesubset = len(usecols) < len(rowdata)
        if usesubset:
            # Successive columns differ by more than one if
            # missing values between
            subset[1:] = usecols[1:] + usecols[:-1] - 1
            subset = np.cumsum(subset)
        subsetids = np.unique(subset)
        levels = set()
        for subsetid in subsetids:
            if usesubset:
                subsetcols = usecols[subset == subsetid]
                sdata = rowdata[subsetcols]
            else:
                sdata = rowdata

            maxval = np.amax(np.absolute(sdata))
            maxdif1 = max4byte
            maxdif2 = max4byte
            if len(sdata) >= mincompress:
                dif1 = sdata[1:] - sdata[:-1]
                maxdif1 = np.amax(np.absolute(dif1))
                dif2 = dif1[1:] - dif1[:-1]
                maxdif2 = np.amax(np.absolute(dif2))
            nbytes = 1 if maxval < max1byte else 2 if maxval < max2byte else 4
            ndif1bytes = 1 if maxdif1 < max1byte else 2 if maxdif1 < max2byte else 4
            ndif2bytes = 1 if maxdif2 < max1byte else 2 if maxdif2 < max2byte else 4
            diflevel = 0
            wdata = sdata
            startdata = []
            if ndif1bytes < nbytes:
                diflevel = 1
                wdata = dif1
                nbytes = ndif1bytes
                startdata = [sdata[0]]
            if ndif2bytes < nbytes:
                diflevel = 2
                wdata = dif2
                nbytes = ndif1bytes
                startdata = [sdata[0], dif1[0]]
            contflag = 1 if subsetid != subsetids[-1] else 0
            subsetflag = 1 if usesubset else 0
            formatflags = (nbytes << 4) + (diflevel << 2) + (contflag << 1) + subsetflag
            gridfile.write(packer.packshort(formatflags))
            if usesubset:
                gridfile.write(packer.packshort(subsetcols[0]))
                gridfile.write(packer.packshort(subsetcols[-1]))
            for v in startdata:
                gridfile.write(packer.packlong(v))
            packfunc = packer.packschar if nbytes == 1 else packer.packshort if nbytes == 2 else packer.packlong
            for i in range(diflevel):
                gridfile.write(packfunc(0))
            for v in wdata:
                gridfile.write(packfunc(v))

    def _writerowv2(self, irow, gridfile, packer):
        rowloc = gridfile.tell()

        for idim in range(self.ndim):
            self._writerowv2a(irow, idim, gridfile, packer)
        return rowloc

    def write(self, gridfile):
        """
        Write the grid to a file handle. 
        """
        packer = _packer(self.format in bigendian_formats)
        offset = gridfile.tell()
        writerow = self._writerowv2 if self.format in version2_formats else self._writerowv1
        gridfile.write(formats[self.format].encode("ascii"))
        # Create a pointer to the file index data, which is written immediately
        # after the pointer in this case (unlike previous perl code)
        indexptrloc = gridfile.tell() + 4 - offset
        gridfile.write(packer.packlong(indexptrloc))
        gridfile.write(packer.packdouble(self.ymin))
        gridfile.write(packer.packdouble(self.ymax))
        gridfile.write(packer.packdouble(self.xmin))
        gridfile.write(packer.packdouble(self.xmax))
        gridfile.write(packer.packdouble(self.resolution))
        gridfile.write(packer.packshort(self.ngridy))
        gridfile.write(packer.packshort(self.ngridx))
        if self.format != "GEOID":
            gridfile.write(packer.packshort(self.ndim))
            gridfile.write(packer.packshort(1 if self.islatlon else 0))
        self._writestring(self.description[0], gridfile, packer)
        self._writestring(self.description[1], gridfile, packer)
        self._writestring(self.description[2], gridfile, packer)
        self._writestring(self.coordsys, gridfile, packer)
        # Save location of row index, and build dummy row index with
        # pointers to 0
        indexptrloc = gridfile.tell()
        for i in range(self.ngridy):
            gridfile.write(packer.packlong(0))
        # Write the rows out and save the location of each row
        rowlocs = [writerow(i, gridfile, packer) - offset for i in range(self.ngridy)]
        # Update the row index in the binary file with the actual locations
        endloc = gridfile.tell()
        gridfile.seek(indexptrloc)
        for rowloc in rowlocs:
            gridfile.write(packer.packlong(rowloc))
        gridfile.seek(endloc)

    def writefile(self, filename):
        """
        Write the grid to a file
        """
        with open(filename, "wb") as gridfile:
            self.write(gridfile)


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Convert CSV grid file to LINZ binary format grid file")
    parser.add_argument("csv_grid_file", help="CSV grid input file")
    parser.add_argument("linz_grid_file", help="LINZ binary grid output file")
    parser.add_argument("-x", "--coordcols", help="Names of lon/lat columns separated by :")
    parser.add_argument("-v", "--datacols", help="Names of data columns separated by :")
    parser.add_argument("-r", "--resolution", type=float, help="Precision of data values")
    parser.add_argument("-f", "--format", choices=sorted(formats), default="GRID2L", help="Binary file format")
    parser.add_argument("-c", "--coordsys", default="NZGD2000", help="Coordinate system code")
    parser.add_argument(
        "-d", "--description", action="append", help="Binary file descriptive header (up to 3 may be specified)"
    )
    parser.add_argument(
        "-s", "--separator", choices=["csv", "space"], default="csv", help="Separator for input file.  Options csv, space"
    )
    parser.add_argument("-n", "--names", help="Column names separated by non alphanumeric - doesn't user header")
    args = parser.parse_args()
    linzgrid = LinzGrid(
        description=args.description,
        coordsys=args.coordsys,
        format=args.format,
        resolution=args.resolution,
        csvfile=args.csv_grid_file,
        crdcols=args.coordcols.split(":") if args.coordcols else None,
        datacols=args.datacols.split(":") if args.datacols else None,
        separator=args.separator,
        colnames=args.names,
    )
    linzgrid.writefile(args.linz_grid_file)


if __name__ == "__main__":
    main()
