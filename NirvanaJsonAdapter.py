import sys
from tempfile import NamedTemporaryFile
import jsonConstants
import jsonStructure
import ijson

class NirvanaJsonAdapter:
    """
    Adapter class for Nirvana Pori import.
    This class is responsible for handling the import of data from Nirvana Pori.
    """
    
    def __init__( self, output_handle=None ):
        """
        Initialize the NirvanaJsonAdapter with an optional output file path.

        :param outputFile: Path to the output JSON file. If None, no output file is set.
        """
        self.output_handle = output_handle
        self.context = {}
        self.simpleMapping = {}
        self.complex_handlers = {}
        self.passFilter = jsonConstants.passFilter
        
    def printOutputHeader(self, patientID, diseaseName, projectName, template="genomic"):
        """
        Print the header for the output JSON.

        :param patientID: Patient ID. Can be "ANONYMOUS" or a specific patient ID.
        :param diseaseName: Disease name for kbDiseaseMatch to match the cancer type.
        :param projectName: Project name.
        :param template: Template for the Pori import. See https://bcgsc.github.io/pori/ipr/templates/".
        """
        print("{", file=self.output_handle)
        print(f'\t"patientId": "{patientID}",', file=self.output_handle)
        print(f'\t"kbDiseaseMatch": "{diseaseName}",', file=self.output_handle)
        print(f'\t"project": "{projectName}",', file=self.output_handle)
        print(f'\t"template": "{template}",', file=self.output_handle)

    def printOutputFooter(self):
        """
        Print the footer for the output JSON.
        """
        print("}", file=self.output_handle)
        
    def setOutputHandle( self, handle ):
        """
        Set the output file for the JSON output.
        This is used to redirect the output to a file or to stdout.
        """
        
        if handle is not None:
            self.output_handle = handle
        else:
            self.output_handle = sys.stdout
            
    # Function to process events
    def processEvents(self, prefix, event, value, ):
        key = (prefix, event)
        path = [x for x in prefix.split('.') if x != 'item'] # Remove all 'item' from the path to make parsing easier
        if key in self.simpleMapping:
            self.handleMapping(path, value)
        elif key in self.complex_handlers:
            self.complex_handlers[key](value)
            

    def handleMapping(self, path, value):
        pathEnd = path[-1]
        container = self.ensureContainerExists( path[:-1] )
        if pathEnd in container: # Check if the key already exists in the container.
            # This is possibly an array item if the item already exists.  Convert it if not already an array.
            if not isinstance(container.get(pathEnd), list):
                oldItem = container.pop(pathEnd)
                container[pathEnd] = [ oldItem ]
            container[pathEnd].append( value )
        else:
            container[pathEnd] = value
            
    def addSimpleMapping(self, key, value):
        """
        Add a simple mapping to the context.
        This is used for key-value pairs that do not require complex handling.
        
        :param key: The JSON event to look for paired with the type. eg ('positions.item.chromosome', 'string')
        :param value: The name that will be associated with the key.  For simple mappings this 
        will generally be the name you want to use in the output JSON.
        """
        self.simpleMapping[key] = value
        
    def addComplexMapping(self, key, handler):
        """
        Add a complex handler for a specific key.
        This is used for keys that require special handling beyond simple mapping.
        Generally this is used for events that require processings such as the start of arrays and maps.
        
        :param key: The JSON event to look for paired with the type. eg ('positions.item', 'end_map')
        :param handler: The function that will handle the complex event.
        """
        self.complex_handlers[key] = handler

    def ensureContainerExists(self, path) -> dict:
        """
        Ensure that the container for the given path exists in the context.
        If it does not exist, create it as an empty dictionary.
        """           
        container = self.context
        for part in path:
            if part not in container:
                container[part] = [{}]
            container = container[part][-1]
        return container

    def handleSimpleMapping(self, key, value):
        collection = key[0].split('.')[-3]
        if collection not in self.context:
            self.context[collection] = {}
        collectionMap = self.context[collection]
        collectionMap.update( {key[0].split('.')[-1]:  value})
        self.context[collection] = collectionMap

    def printComma( self ):
        """
        Print a comma if the iterator is greater than 0.
        This is used to separate JSON objects in an array.
        """
        iterator = self.context.get('iterator', 0)
        if iterator > 0:
            print(",", end=' ', file=self.output_handle)
        self.context['iterator'] = iterator + 1
    
    def addArrayToContext(self, path):
        jsonStructure.addArrayToContext(self.context, path)

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
    
    def getOutputHandle(self):
        return self.output_handle if hasattr(self, 'output_handle') else sys.stdout
    
    def printHeader(self):
        pass
    
    def readJsonFile(self, jsonFile ):
        self.iterator= 0
        self.context['positions'] = [{}]
        
        originalOutputHandle = self.getOutputHandle()
        # Read the JSON file
        with open( jsonFile, 'r') as f:
            parser = ijson.parse( f )
            #position = None
            self.printHeader()
            
            with NamedTemporaryFile( mode='w+', delete=True) as tempfile:
                self.setOutputHandle( tempfile )
                print( "[", file=tempfile )
                for prefix, event, value in parser:
                    self.processEvents( prefix, event, value )
                print( "]", file=tempfile )
                
                # Go back to the start of the temporary file and read it for the 2nd pass.
                tempfile.seek(0)
                self.setOutputHandle( originalOutputHandle ) # Not sure if this is needed
                jsonStructure.perform2ndPass( tempfile, originalOutputHandle )