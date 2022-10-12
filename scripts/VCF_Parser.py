from collections import OrderedDict
import datetime
import sys


# VCF file function
def VariantFile(file_name, mode='r', header=None, progName = None, progVersion = None, progCmd = None):
    if mode not in ('r', 'w'):
        raise Exception(f"Unrecognized mode {mode}")
    elif mode == 'r':
        return VariantReader(file_name)
    elif mode == 'w':
        if header is None:
            raise Exception("No header provided for VariantWriter creation")
        else:
            return VariantWriter(file_name, header, progName, progVersion, progCmd)


# VCF parser class
class VariantReader:
    def __init__(self, file_name):
        self.source = open(file_name, 'r')
        header = []
        line = next(self.source)
        while line[:2] == '##':
            header.append(line)
            line = next(self.source)
        if line[0] != '#':
            raise Exception("Labels line missing")
        header.append(line)
        self.header = VariantHeader(header)
        linebins = line.strip().split()
        self.samps = linebins[9:]

    def __iter__(self):
        return self

    def __next__(self):
        try:
            return VariantRecord(next(self.source), self.samps)
        except StopIteration:
            raise StopIteration

    next = __next__

    def close(self):
        self.source.close()


# VCF writer class
class VariantWriter:
    def __init__(self, file_name, header, progName = None, progVersion = None, progCmd = None):
        self.write_file = open(file_name, 'w')
        myHeader = header
        if not(progName is None or progVersion is None or progCmd is None):
            myHeader.addLine("prog", progName, vNum=progVersion, description=progCmd)
        myHeader.UpdateDateline()
        self.write_file.write(str(myHeader))

    def writeline(self, vcfLine):
        self.write_file.write(str(vcfLine))

    def close(self):
        self.write_file.close()


class VariantHeader:
    def __init__(self, in_lines):
        self.versionLine = []
        self.dateline = []
        self.progLines = []
        self.formatLines = []
        self.contigLines = []
        self.infoLines = []
        self.filterLines = []
        self.altLines = []
        self.otherLines = []
        self.labelLine = []

        for line in in_lines:
            if "##fileformat" in line:
                self.versionLine.append(line)
            elif "##filedate" in line:
                self.dateline.append(line)
            elif "##prog" in line:
                self.progLines.append(line)
            elif "##contig" in line:
                self.contigLines.append(line)
            elif "##INFO" in line:
                self.infoLines.append(line)
            elif "##FORMAT" in line:
                self.formatLines.append(line)
            elif "##FILTER" in line:
                self.filterLines.append(line)
            elif "##ALT" in line:
                self.altLines.append(line)
            elif "#CHROM" in line:
                self.labelLine.append(line)
            else:
                self.otherLines.append(line)
                sys.stderr.write(f"WARN: Unrecognized header line type for line {line}.\n")

    def addLine(self, lineType, label, number='.', Type="String", description="", source="", vNum=""):
        if lineType.upper() == "FORMAT":
            self.formatLines.append(
                f'##FORMAT=<ID={label},Number={number},Type={Type},Description="{description}">\n')
        elif lineType.upper() == "INFO":
            self.infoLines.append(
                f'##INFO=<ID={label},Number={number},Type={Type},Description="{description}",Source="{source}",Version="{vNum}">\n')
        elif lineType.upper() == "FILTER":
            self.filterLines.append(
                f'##FILTER=<ID={label},Description="{description}">\n')
        elif lineType.upper() == "ALT":
            self.altLines.append(
                f'##ALT=<ID={label},Description="{description}">\n')
        elif lineType == "prog":
            self.progLines.append(
                f'##prog=<ID={label},Version={vNum},Command={description}>\n')
        else:
            raise Exception(f"Unrecognized header line type {lineType}.\n")
            
    def UpdateDateline(self):
        mydate = datetime.date.today()
        self.dateline = [f"##filedate={mydate.year}{mydate.month:02d}{mydate.day:02d}\n"]
    

    def __str__(self):
        headLines = []
        headLines.extend(self.versionLine)
        headLines.extend(self.dateline)
        headLines.extend(self.progLines)
        headLines.extend(self.contigLines)
        headLines.extend(self.infoLines)
        headLines.extend(self.formatLines)
        headLines.extend(self.filterLines)
        headLines.extend(self.altLines)
        headLines.extend(self.otherLines)
        headLines.extend(self.labelLine)
        return ''.join(headLines)


# VCF line class
class VariantRecord:
    def __init__(self, in_line, sampleNames):
        linebins = in_line.strip().split()
        self.chrom = linebins[0]
        self.pos = int(linebins[1])
        self.id = linebins[2]
        self.ref = linebins[3]
        self.alts = [x for x in linebins[4].split(',')]
        self.qual = int(linebins[5]) if linebins[5] is not '.' else None
        if linebins[6] == '.':
            self.filter = []
        else:
            self.filter = linebins[6].split(';')
        if linebins[7] == '.':
            self.info = {}
        else:
            self.info = {x.split('=')[0]: x.split('=')[1] for x in linebins[7].split(";")}
        self.format = [x for x in linebins[8].split(':')]
        sampBins = [x.split(':') for x in linebins[9:]]
        self.samples = OrderedDict(
            (sampleNames[x], OrderedDict(
                (self.format[y], sampBins[x][y]) for y in range(len(self.format))
            )) for x in range(len(sampleNames))
        )

    def __str__(self):
        outStr = [
            str(self.chrom),
            str(self.pos),
            str(self.id),
            str(self.ref),
            ",".join(self.alts),
            str(self.qual) if self.qual is not None else '.',
            ';'.join(sorted(self.filter)) if len(self.filter) != 0 else '.',
            ";".join([f"{x}={self.info[x]}" for x in self.info]) if len(self.info) != 0 else '.',
            ':'.join(self.format)
        ]
        for samp in self.samples:
            outStr.append(':'.join(self.samples[samp][x] for x in self.samples[samp]))
        retStr = '\t'.join(outStr)
        return f"{retStr}\n"

    def add_filter(self, new_filter):
        if new_filter in self.filter:
            return False
        else:
            self.filter.append(new_filter)
            return True

    def get_filters(self):
        return self.filter

    def has_filter(self, test_filter):
        return test_filter in self.filter

    def remove_filter(self, test_filter):
        if test_filter in self.filter:
            self.filter.remove(test_filter)
            return True
        else:
            return False

    def add_info(self, new_info_tag, new_info_value):
        if new_info_tag in self.info:
            raise Exception
        else:
            self.info[new_info_tag] = new_info_value
            return True
            
def SamHeaderToVcfHeader(inSamHeader, sampName, progName, progVersion, progCmd):
    mydate = datetime.date.today()
    headLines = [
        "##fileformat=VCFv4.3\n", 
        f"##filedate={mydate.year}{mydate.month:02d}{mydate.day:02d}\n"
        ]
    contigBlock = []
    programBlock = []
    
    for line in str(inSamHeader).split('\n'):
        linebins = line.strip().split('\t')
        if linebins[0] == '@HD':
            pass
        elif linebins[0] == '@SQ':
            contigBlock.append(
                f"##contig=<ID={linebins[1].split(':')[1]},length={linebins[2].split(':')[1]}>\n"
            )
        elif linebins[0] == '@PG':
            progHead = None
            progVersion = None
            progCL = None
            for binIter in range(len(linebins)):
                if 'ID' in linebins[binIter]:
                    progHead = f"{linebins[binIter].split(':')[1].split()[0]}"
                elif 'VN' in linebins[binIter]:
                    progVersion = linebins[binIter].split(':')[1]
                elif 'CL' in linebins[binIter]:
                    progCL = " :".join(linebins[binIter:]).split(':')[1]
            if progHead is not None:
                if progVersion is not None:
                    programBlock.append(f'##prog=<ID={progHead},'
                                        f'Version={progVersion},'
                                        f'Command="{progCL}">\n')
    headLines.extend(programBlock)
    headLines.append(f'##prog=<ID={progName},Version={progVersion},Command="{progCmd}">\n')
    headLines.extend(contigBlock)
    headLines.append(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sampName}\n')
    return VariantHeader(headLines)