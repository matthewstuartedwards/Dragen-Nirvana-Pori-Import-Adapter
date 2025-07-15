import json
import sys
from tempfile import NamedTemporaryFile
import ijson
from NirvanaJsonAdapter import NirvanaJsonAdapter
from jsonStructure import perform2ndPass
from jsonConstants import cnvConsequencePriorityList, variantTypeKBCategoryMap

class CnvAdapter(NirvanaJsonAdapter):
    """
    Adapter class for CNV data import.
    This class is responsible for handling the import of CNV data from Nirvana Pori.
    """
    
    # private class members
    iterator= 0
    positions = [{}]
    printedGenes = set()  # Set to keep track of printed genes to avoid duplicates
    currentTranscript = None
    
    def __init__(self, cnv_file):
        """
        Initialize the CnvAdapter with the given CNV file path.

        :param cnv_file: Path to the input JSON file with CNV information.
        """
        super().__init__()
        self.cnvFile = cnv_file
        

        
        # Add mappings
        self.addSimpleMapping(('positions.item.chromosome', 'string'), 'chromosome')
        self.addSimpleMapping(('positions.item.position', 'number'), 'start')
        self.addSimpleMapping(('positions.item.svEnd', 'number'), 'end')
        self.addSimpleMapping(('positions.item.cytogeneticBand', 'string'), 'cytogeneticBand')
        self.addSimpleMapping(('positions.item.variants.item.variantType', 'string'), 'kbCategory')
        self.addSimpleMapping(('positions.item.variants.item.transcripts.item.transcript', 'string'), 'transcript')
        self.addSimpleMapping(('positions.item.variants.item.transcripts.item.source', 'string'), 'source')
        self.addSimpleMapping(('positions.item.variants.item.transcripts.item.hgnc', 'string'), 'gene')
        self.addSimpleMapping(('positions.item.samples.item.copyNumber', 'number'), 'copyNumber')
        self.addSimpleMapping(('positions.item.samples.item.minorHaplotypeCopyNumber', 'number'), 'minorHaplotype')
        self.addSimpleMapping(('positions.item.samples.item.lossOfHeterozygosity', 'number'), 'lossOfHeterozygosity')
        self.addSimpleMapping(('positions.item.samples.item.genotype', 'string'), 'genotype')
        self.addSimpleMapping(('positions.item.variants.item.transcripts.item.isCanonical', 'boolean'), 'isCanonical')
        self.addSimpleMapping(('positions.item.variants.item.transcripts.item.completeOverlap', 'boolean'), 'completeOverlap')
        self.addSimpleMapping(('positions.item.variants.item.transcripts.item.consequence.item', 'string'), 'consequence')
        self.addSimpleMapping(('positions.item.filters.item', 'string'), 'filters')
        
        self.addComplexMapping(('positions.item', 'end_map'), self.handle_end_map_positions_item)
        self.addComplexMapping(('positions.item.variants.item', 'start_map'), self.handle_start_map_variants_item)
        self.addComplexMapping(('positions.item.variants.item.transcripts.item.consequence.item', 'start_array'), self.handleTranscriptConsequence)
        self.addComplexMapping(('positions.item.variants.item.transcripts.item', 'start_map'), self.handleNewTranscript)
    
    # When a position ends, it's time to figure out if we need to print the position.
    # We don't print a position if it doesn't have the right filter item, or if the gene has already been printed.
    # PORI only maps gene names to the database, so we can't have multiple genes being output.
    def handle_end_map_positions_item(self, value):
        position = self.positions[0]
        
        if position.get('filters') and self.passFilter in position['filters']:
            printPosition = self.massagePosition( position )
            
            # Check if the position has a gene as PORI will expect one.
            if 'gene' not in printPosition or printPosition['gene'] is None:
                self.positions = [{}]  # Reset positions to avoid printing empty objects
                return
            
            self.printComma() # Function handles printing a comma if needed for array or map purposes
            print(json.dumps(printPosition, indent=4), file=self.outputFile)
        
        self.positions = [{}]
        
    # This function handles the start of a new transcript item
    # If this is the first new transcript, it initializes the list and the context
    def handleNewTranscript(self, value):
        self.addArrayToContext(['positions','variants','transcripts'])

    def handle_start_map_variants_item(self, value):
        self.addArrayToContext( ['positions', 'variants'] )
        
    def setOutputHandle(self, handle):
        return super().setOutputHandle(handle)

    def printCNVHeader( self ):
        print( "\t\"copyVariants\": [" )
        
    def getOutputHandle(self):
        return self.outputFile if hasattr(self, 'outputFile') else sys.stdout
        
    def readCnvFile(self, cnvJsonFile ):
        self.iterator= 0
        self.positions = [{}]
        
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
                self.setOutputHandle( originalOutputHandle )
                perform2ndPass( tempfile )
                
    def massagePosition(self, position):
        """
        Massage the position data to prepare it for output.
        This includes converting chromosome band, handling variants, and removing unnecessary fields.
        """
        printPosition = position.copy()
        
        # Convert chromosome band to a more readable format
        printPosition['chromosomeBand'] = self.convertCytogeneticBand(printPosition.get('chromosome'), printPosition.pop('cytogeneticBand', None))
        
        if 'transcripts' in printPosition['variants'][0]:
            self.processVariant( printPosition )
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


    def processVariant(self, position):
        #for variant in position.get('variants'):
        variant = position.get('variants')[0] #TODO, handle multiple variants
        transcript = self.getBestTranscript( variant.get('transcripts') )
        position['transcript'] = transcript['transcript']
        position['kbCategory'] = variant['variantType']
        position['kbCategory'] = variantTypeKBCategoryMap.get(position.get('kbCategory'), 'unknown')
        position['gene'] = transcript['hgnc']
        position['source'] = transcript['source']
        
        position.pop('variants', None)  # Remove variants as we are only interested in the first one for now

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
            consequenceRank.append( (self.getBestConsequence( transcript.get('consequence'), cnvConsequencePriorityList), transcript.get('source') ) )
        bestRank = len(cnvConsequencePriorityList)
        bestIndex = -1
        i = 0
        for (rank, source) in consequenceRank:
            
            canonical = transcripts[i].get('isCanonical', False)
            bestCanonical = transcripts[bestIndex].get('isCanonical', False) if bestIndex >= 0 else False
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
    
    def getBestConsequence(self, consequences, consequencePriorityList):
        """
        Get the best consequence from the list of transcript events.
        """
        bestRank = len(consequencePriorityList)
        indexOfBestRank = -1
        index = 0
        for event in consequences:
            #print( "Event: " + str(event), file=sys.stderr )
            
            if event in consequencePriorityList:
                rank = consequencePriorityList.index( event )
                if rank < bestRank:
                    bestRank = rank
                    indexOfBestRank = index
            index+=1
        return indexOfBestRank
    
    def handleTranscriptConsequence(self, value):
        if not self.currentTranscript:
            self.currentTranscript = { 'consequence': [value]}
        else:
            self.currentTranscript['consequence'].append(value)
        
        #transcript = self.currentTranscript
        #if 'consequence' not in transcript:
            #transcript['consequence'] = [value]
        #else:
            #transcript['consequence'].append(value)
        self.currentTranscript['consequenceRank'] = cnvConsequencePriorityList.index(value) if value in cnvConsequencePriorityList else len(cnvConsequencePriorityList)
