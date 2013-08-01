'''
    The MobiDB offline consensus generator is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

from Bio import SeqIO
import math

SOURCES = ['disprot','pdb-xray','pdb-nmr','espritz-xray','espritz-nmr','espritz-disprot','iupred-short','iupred-long','ideal','disembl-465','disembl-HL']

class DisConsensus:
    '''
    Class to calculate disorder consensus sequences from MobiDB datasets.
    '''
    def __init__(self,resthres=5,exclusions=[],inclusions=[]):
        self.refseqlen = None
        self.resthres = resthres
        self.exclusions = exclusions
        self.inclusions = inclusions
        self.levels = {}
        self.scores = {}
        self.annfilepath = None
        self.savememory = False

    def getweight(self,annsrc,res=None,thres=None):
        ''' 
       Method that returns the weight for a given annotation source. Fixed weights
       for all sources except pdb-xray, which "awards" higher resolution structures.
       Returns a float containing the weight.

       annsrc (string): name of the annotation source ('pdb-xray', 'espritz-nmr', etc.)
       res (float): resolution of the structure, used only if annsrc='pdb-xray'
       thres (float): resolution threshold, used only if annsrc='pdb-xray'
       '''
        if annsrc == 'pdb-xray':
            weight = 1-math.log(res)/math.log(thres)
            if weight < 0:
                weight = 0
            elif weight < 0.2:
                weight = 0.2
        elif annsrc == 'pdb-nmr':
            weight = 0.2
        elif annsrc == 'disprot' or annsrc == 'ideal':
            weight = 3
        else:
            weight = 0.05

        return(weight)

    def getannsdict(self,annfilepath):
        '''
        Method that parses "annotation.fasta" files from MobiDB datasets and builds a dictionary
        of the form {'refseqid':{'annsrc1':[['annseq1id','annseq1str'],[...]],'annsrc2':[[...]]}}.
        If the annsrc = 'pdb-xray', the leaf lists have a third member for the resolution.
        Returns the dictionary.

        annfilepath (string): path to the annotations fasta file
        '''

        if not annfilepath:
            raise ValueError("annfilepath not set")

        records = SeqIO.index(annfilepath, "fasta")

        rseqs = {}

        for r in records:
            r = records[r]

            header = r.id.split("|")
            refid = header[0]
            annid = header[2]
            annsrc = header[3]
            anntype = header[4]

            if annsrc not in SOURCES or annsrc in self.exclusions:
                continue

            if len(self.inclusions) > 0 and annsrc not in self.inclusions:
                continue

            if annsrc == 'pdb-xray':
                try:
                    res = header[5]
                    entry = [annid,str(r.seq),float(res)]
                except IndexError:
                    print("nope: ",annid)
            else:
                entry = [annid,str(r.seq)]
            
            try:
                rseqs[refid][annsrc].append(entry)
            except KeyError:
                try:
                    rseqs[refid][annsrc] = [entry]
                except KeyError:
                    rseqs[refid] = {annsrc:[entry]}

        return(rseqs)

    def calculate(self):
        '''
        Method to calculate the disorder consensus for each reference sequence in the annotations file from MobiDB.
        Refer to the MobiDB documentation for details on the calculations.
        Returns nothing. Fills the "levels" and "scores" dictionary properties.
        '''

        if not self.annfilepath:
            raise ValueError("annfilepath property not set")

        rseqs = self.getannsdict(self.annfilepath) #generate the dictionary containing all annotations for each ref. sequence

        for rseqid in rseqs:
            annsrckey = rseqs[rseqid].keys()[0]
            self.refseqlen = len(rseqs[rseqid][annsrckey][0][1]) #get the first aligned sequence's length, which is the same as the ref seq's

            levels = [0] * self.refseqlen #disorder levels for each position in the refseq
            dscores = [0] * self.refseqlen #disorder score for each position
            sscores = [0] * self.refseqlen #structure score for each position
            scores = [0] * self.refseqlen #disorder (if dscores[i]>=sscores[i]) or structure (otherwise) score for each position in the refseq

            for annsrcid in rseqs[rseqid]: #for each annotation source
                for seq in rseqs[rseqid][annsrcid]: #for each annotating sequence
                    seqid = seq[0]
                    seqstr = seq[1]
                    seqres = seq[2] if annsrcid == 'pdb-xray' else None

                    if annsrcid == 'pdb-xray' and seqres > self.resthres: #avoid pdb xray structures with a resolution higher than the threshold
                        continue

                    w = self.getweight(annsrcid,seqres,self.resthres) #calculate the weight for the annotating sequence

                    for i in xrange(0,self.refseqlen): #assign the weight to each position, depending on if it's voted to be d or s
                        if seqstr[i] == '1': #disordered
                            dscores[i] += w
                        elif seqstr[i] == '0': #structured
                            sscores[i] += w

            for i in xrange(0,self.refseqlen): #calculate the levels and scores for each position
                if dscores[i] + sscores[i] == 0: #if no scores at all
                    levels[i] = None 
                else:
                    levels[i] = int(round(min(10*float(dscores[i])/(dscores[i]+sscores[i]),9)))
                    scores[i] = round(max(sscores[i],dscores[i]),2)

            if self.savememory:#output the data without storing it as properties, to save memory. Useful for large annotation files
                print(">{0}\n{1}".format(rseqid,"".join([str(c).replace("None","-") for c in levels])))
            else:
                self.levels[rseqid] = levels
                self.scores[rseqid] = scores

def parsearguments():
    import argparse

    description = """MobiDB offline consensus disorder calculator"""
    epilog = """For obtaining datasets and documentation, see http://mobidb.bio.unipd.it. """
    parser = argparse.ArgumentParser(description=description,epilog=epilog,formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('annfile',metavar="ANNFILE", help='the MobiDB annotations fasta file from which to extract the annotations')
    parser.add_argument('-r','--resolution-threshold', metavar='N', help='the maximum resolution for a pdb-xray structure to be used in the calculation', default=5, type=float)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--include',metavar='SOURCE',help='list of the sources to (exclusively) include in the calculation',nargs="*",choices=SOURCES,default=[])
    group.add_argument('--exclude',metavar='SOURCE',help='list of the sources to exclude from the calculation',nargs="*", choices=SOURCES,default=[])

    return(vars(parser.parse_args()))

if __name__ == "__main__":
    args = parsearguments()
    annsf = args['annfile']
    rt = args['resolution_threshold']

    dc = DisConsensus()
    dc.annfilepath = annsf
    dc.resthres = rt
    dc.inclusions = args['include']
    dc.exclusions = args['exclude']
    dc.savememory = True

    dc.calculate()
