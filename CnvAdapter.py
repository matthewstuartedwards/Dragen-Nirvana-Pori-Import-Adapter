import json
from NirvanaJsonAdapter import NirvanaJsonAdapter

class CnvAdapter(NirvanaJsonAdapter):
    """
    Adapter class for CNV data import.
    This class is responsible for handling the import of CNV data from Nirvana Pori.
    """
    
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
        position = self.context['positions'][0]
        
        printedGenes = self.context.get('printedGenes', set())
        if position.get('filters') and self.passFilter in position['filters']:
            printPosition = self.massagePosition( position )
            
            # Check if the position has a gene as PORI will expect one.
            if 'gene' not in printPosition or printPosition['gene'] is None:
                self.context['positions'] = [{}]  # Reset positions to avoid printing empty objects
                return
            
            # Check if the gene has already been printed
            # COMMENTING THIS OUT FOR NOW.  GOING TO USE A 2ND PASS TO REMOVE DUPLICATE GENES.
            #
            #if printPosition['gene'] not in printedGenes:
            #    printedGenes.add(printPosition['gene'])
            #    context['printedGenes'] = printedGenes
            #else:
            #    context['positions'] = [{}]  # Reset positions to avoid printing empty objects
            #    return
            
            self.printComma() # Function handles printing a comma if needed for array or map purposes
            print(json.dumps(printPosition, indent=4), file=outputFile)
        
        self.context['positions'] = [{}]
        
    # This function handles the start of a new transcript item
    # If this is the first new transcript, it initializes the list and the context
    def handleNewTranscript(self, value):
        self.addArrayToContext(['positions','variants','transcripts'])

    def handle_start_map_variants_item(self, value):
        self.addArrayToContext( ['positions', 'variants'] )
        