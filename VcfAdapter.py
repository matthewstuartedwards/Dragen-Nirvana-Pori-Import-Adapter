import json
import sys
from tempfile import NamedTemporaryFile
import ijson
from NirvanaJsonAdapter import NirvanaJsonAdapter
from jsonStructure import perform2ndPass
from jsonConstants import cnvConsequencePriorityList

class Transcript:
    """
    Class to represent a transcript.
    This class is used to store information about a transcript, including its ID, source, and consequences.
    """
    hgnc = None
    hgvsp = None
    hgvsc = None
    isCanonical = False
    source = None
    consequences = []
    transcript = None
    
    def __init__( self ):
        self.consequences = []
    
    def putSource(self, source):
        self.source = source
    
    def putConsequence(self, consequence):
        self.consequences.append(consequence)
        
    def putHgnc( self, hgnc ):
        self.hgnc = hgnc
        
    def putHgvsp(self, hgvsp):
        self.hgvsp = hgvsp
    
    def putHgvsc(self, hgvsc):
        self.hgvsc = hgvsc
        
    def putCanonical(self, isCanonical):
        self.isCanonical = isCanonical
        
    def putTranscript(self, transcript):
        self.transcript = transcript
        
    def getConsequences(self):
        return self.consequences
    
    def getSource(self):
        return self.source
    
    def getGene(self):
        return self.hgnc
    
    def getCanonical(self):
        return self.isCanonical
    
    def getTranscript(self):
        return self.transcript
    
    def getHgvsProtein(self):
        return self.hgvsp
    
    def getHgvsCds(self):
        return self.hgvsc


class VcfAdapter(NirvanaJsonAdapter):
    """
    Adapter class for CNV data import.
    This class is responsible for handling the import of CNV data from Nirvana Pori.
    """
    
    # private class members
    iterator= 0
    transcriptEvents = []
    transcripts = []
    currentTranscript = None
    
    def __init__(self, output_handle):
        """
        Initialize the VcfAdapter with the given VNC file path.

        :param cnv_file: Path to the input JSON file with VNC information.
        """
        super().__init__()
        self.context['positions'] = [{}]
        self.setOutputHandle(output_handle)
        
        # Add mappings
        self.addSimpleMapping(('positions.item.chromosome', 'string'), 'chromosome')
        self.addSimpleMapping(('positions.item.variants.item.begin', 'number'), 'startPosition')
        self.addSimpleMapping(('positions.item.variants.item.end', 'number'), 'endPosition')
        self.addSimpleMapping(('positions.item.variants.item.refAllele', 'string'), 'refSeq')
        self.addSimpleMapping(('positions.item.variants.item.altAllele', 'string'), 'altSeq')
        self.addSimpleMapping(('positions.item.variants.item.transcripts.item.hgvsp', 'string'), 'hgvsp')
        self.addSimpleMapping(('positions.item.variants.item.transcripts.item.hgvsc', 'string'), 'hgvsc')
        self.addSimpleMapping(('positions.item.variants.item.transcripts.item.isCanonical', 'boolean'), 'isCanonical')
        self.addSimpleMapping(('positions.item.variants.item.transcripts.item.source', 'string'), 'source')
        
        self.addComplexMapping(('positions.item', 'end_map'), self.handle_end_map_positions_item)
        #self.addComplexMapping(('positions.item.variants.item', 'start_map'), self.handle_start_map_variants_item)
        self.addComplexMapping(('positions.item.variants.item.transcripts.item.consequence.item', 'string'), self.handleTranscriptConsequence)
        self.addComplexMapping(('positions.item.variants.item.transcripts.item', 'start_map'), self.handleNewTranscript)
        self.addComplexMapping(('positions.item.variants.item.transcripts.item', 'end_map'), self.addFinishedTranscript)
        
    
    # When a position ends, it's time to figure out if we need to print the position.
    # We don't print a position if it doesn't have the right filter item, or if the gene has already been printed.
    # PORI only maps gene names to the database, so we can't have multiple genes being output.
    def handle_end_map_positions_item(self, value):
        
        position = self.context['positions'][0]
        
        if position.get('filters') and self.passFilter in position['filters']:
            printPosition = self.massagePosition( position, self.transcripts )
            
            
            # Check if the position has a gene as PORI will expect one.
            if 'gene' not in printPosition or printPosition['gene'] is None:
                self.context['positions'] = [{}]  # Reset positions to avoid printing empty objects
                return
            
            self.printComma() # Function handles printing a comma if needed for array or map purposes
            print(json.dumps(printPosition, indent=4), file=self.output_handle)
        
        self.context['positions'] = [{}]
        
    # This function handles the start of a new transcript item
    # If this is the first new transcript, it initializes the list and the context
    def handleNewTranscript(self, value):
        self.currentTranscript = Transcript()
        
    def addFinishedTranscript(self, value):
        self.currentTranscript.putSource( self.context['positions.item.variants.item.transcripts.item.source'] )
        self.currentTranscript.putHgnc( self.context['positions.item.variants.item.transcripts.item.hgnc'] )
        self.currentTranscript.putHgvsp( self.context['positions.item.variants.item.transcripts.item.hgvsp'] )
        self.currentTranscript.putHgvsc( self.context['positions.item.variants.item.transcripts.item.hgvsc'] )
        self.currentTranscript.putCanonical( self.context['positions.item.variants.item.transcripts.item.isCanonical'] )
        self.currentTranscript.putTranscript( self.context['positions.item.variants.item.transcripts.item.transcript'] )
        
        
        self.transcripts.append(self.currentTranscript)

#    def handle_start_map_variants_item(self, value):
#        self.addArrayToContext( ['positions', 'variants'] )
        
    def setOutputHandle(self, handle):
        return super().setOutputHandle(handle)

    def printCNVHeader( self ):
        print( "\t\"copyVariants\":", file=self.output_handle )
        
    def getOutputHandle(self):
        return self.output_handle if hasattr(self, 'output_handle') else sys.stdout
        
    def readCnvFile(self, cnvJsonFile ):
        self.iterator= 0
        self.context['positions'] = [{}]
        
        originalOutputHandle = self.getOutputHandle()
        # Read the JSON file
        with open( cnvJsonFile, 'r') as f:
            parser = ijson.parse( f )
            #position = None
            self.printCNVHeader()
            
            with NamedTemporaryFile( mode='w+', delete=True) as tempfile:
                self.setOutputHandle( tempfile )
                print( "[", file=tempfile )
                for prefix, event, value in parser:
                    self.processEvents( prefix, event, value )
                print( "]", file=tempfile )
                
                # Go back to the start of the temporary file and read it for the 2nd pass.
                tempfile.seek(0)
                self.setOutputHandle( originalOutputHandle ) # Not sure if this is needed
                perform2ndPass( tempfile, originalOutputHandle )

    def massagePosition(self, position, transcripts):
        """
        Massage the position data to prepare it for output.
        This includes converting chromosome band, handling variants, and removing unnecessary fields.
        """
        printPosition = position.copy()
    
        
        if transcripts:
            bestTranscript = self.getBestTranscript(transcripts)
            self.processTranscript( printPosition, bestTranscript )
        if 'variants' in printPosition: # Situation where variants are still around because there were no transcripts
            printPosition.pop('variants', None)
        if 'samples' in printPosition:
            self.processSample( printPosition )
        
        return printPosition
    
    def convertCytogeneticBand(self, chromosome, cytogeneticBand):
        return self.getChromosomeNumberOnly(chromosome) + ":" + cytogeneticBand if cytogeneticBand else None
    
    def getChromosomeNumberOnly(self, chromosome):
        if 'chr' in chromosome:
            return chromosome.replace('chr', '')
        return chromosome

    def processSample(self, position):
        if len(position.get('samples')) > 1:
            print("More than one sample found, only processing the first one.", file=sys.stderr)
        
        
        sample = position.get('samples')[0]
        
        copyNumber = sample.get('copyNumber', None)
        minorHaplotype = sample.get('minorHaplotypeCopyNumber', None)
        genotype = sample.get('genotype', None)
        
        #cna = copyNumber / minorHaplotype 
        copyChange = (copyNumber - 2)
        #log2Cna = math.log2(cna)
        lohState = sample.get('lossOfHeterozygosity', None)
        if lohState is not None:
            position['lohState'] = "LOH"
        #else:
        #    position['lohState'] = "HET"
            
        #position['cna'] = cna
        position['copyChange'] = copyChange
        #position['log2Cna'] = log2Cna
        

        position.pop('samples', None)  # Remove samples as we are only interested in the first one for now


    def processTranscript(self, position, transcript):
        position['gene'] = transcript.getHgnc()
        position['source'] = transcript.getSource()
        position['isCanonical'] = transcript.getCanonical()
        position['transcript'] = transcript.getTranscript()
        position['hgvsProtein'] = transcript.getHgvsp()
        position['hgvsCds'] = transcript.getHgvsc()

    def getBestTranscript(self, transcripts):
        """
        Get the best transcript from the list of transcripts.
        The best transcript is the one that is canonical and has the lowest consequence rank
        TODO: Check if the best transcript should also be a completeOverlap transcript.
        """
        if not transcripts:
            return {}
        
        consequenceRank = []
        for transcript in transcripts:
            consequenceRank.append( (self.getBestConsequence( transcript.getConsequence(), cnvConsequencePriorityList), transcript.getSource() ) )
        bestRank = len(cnvConsequencePriorityList)
        bestIndex = -1
        i = 0
        for (rank, source) in consequenceRank:
            
            canonical = transcripts[i].getCanonical()
            bestCanonical = transcripts[bestIndex].getCanonical() if bestIndex >= 0 else False
            # Check for the canonical transcript first.  This should be one of the best transcripts to use.
            # If the i indexed transcript is canonical and the best one is not, then we should use this one.
            if canonical and not bestCanonical: 
                bestRank = rank
                bestIndex = i
            elif rank < bestRank: # standard checking for best ranking
                bestRank = rank
                bestIndex = i
            elif rank == bestRank and source == 'RefSeq': # This can be an elif. Check if the source is our preferred source
                bestRank = rank
                bestIndex = i
            i += 1
        
        return transcripts[bestIndex]
    
    def handleTranscriptConsequence(self, value):
        self.currentTranscript.putConsequence(value)
        
    def printHeader(self):
        print("\t\"smallMutations\":", file=self.output_handle)

