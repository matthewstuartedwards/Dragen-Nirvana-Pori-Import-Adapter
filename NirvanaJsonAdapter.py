import sys
import jsonConstants
import jsonStructure

class NirvanaJsonAdapter:
    """
    Adapter class for Nirvana Pori import.
    This class is responsible for handling the import of data from Nirvana Pori.
    """
    
    def __init__( self, outputFile=None ):
        """
        Initialize the NirvanaJsonAdapter with an optional output file path.

        :param outputFile: Path to the output JSON file. If None, no output file is set.
        """
        self.outputFile = outputFile
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
        print("{")
        print(f'\t"patient_id": "{patientID}",')
        print(f'\t"disease_name": "{diseaseName}",')
        print(f'\t"project_name": "{projectName}",')
        print(f'\t"template": "{template}",')

    def printOutputFooter(self):
        """
        Print the footer for the output JSON.
        """
        print("}")
        
    def setOutputHandle( self, handle ):
        """
        Set the output file for the JSON output.
        This is used to redirect the output to a file or to stdout.
        """
        
        if handle is not None:
            self.outputFile = handle
        else:
            self.outputFile = sys.stdout
            
    # Function to process events
    def processEvents(self, context, prefix, event, value, ):
        key = (prefix, event)
        path = [x for x in prefix.split('.') if x != 'item'] # Remove all 'item' from the path to make parsing easier
        if key in self.simpleMapping:
            self.handleMapping(path, value)
        elif key in self.complex_handlers:
            self.complex_handlers[key](context, value)
            

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
            print(",", end=' ', file=self.outputFile)
        self.context['iterator'] = iterator + 1
    
    def addArrayToContext(self, path):
        jsonStructure.addArrayToContext(self.context, path)