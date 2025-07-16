import json
from collections import defaultdict

def addArrayToContext(context, path):
    container = context
    pathEnd = path[-1]
    for part in path[:-1]:
        container = container[part][-1]
    
    if pathEnd in container:
        container[pathEnd].append({})
    else:
        container[pathEnd] = [{}]
        
def perform2ndPass(input_handle, output_handle):
    # Load your JSON data (replace this with loading from a file if needed)
    #with open(output_file, 'r') as input:
    data = json.load(input_handle)

    # Group entries by gene
    gene_entries = defaultdict(list)
    for entry in data:
        gene_entries[entry['gene']].append(entry)
        
    # Select the best entry for each gene
    selected_entries = []
    for gene, entries in gene_entries.items():
        # Priority 1: kbCategory is 'deep deletion'
        deep_deletion = [e for e in entries if e.get('kbCategory') == 'deep deletion']
        if deep_deletion:
            selected_entries.append(deep_deletion[0])
            continue
        
        # Priority 2: source is 'RefSeq' and isCanonical is True
        refseq_canonical = [e for e in entries if e.get('source') == 'RefSeq' and e.get('isCanonical') is True]
        if refseq_canonical:
            selected_entries.append(refseq_canonical[0])
            continue
        
        # Priority 3: source is 'RefSeq'
        refseq = [e for e in entries if e.get('source') == 'RefSeq']
        if refseq:
            selected_entries.append(refseq[0])
            continue
        
        # Fallback: first available entry
        selected_entries.append(entries[0])

    # Output the selected entries
    #print( "\t\"smallMutations\": ", end = "" )
    print( json.dumps(selected_entries, indent=4), end = "", file=output_handle )
    #print( "," )


def ensureContainerExists(context, path):
    """
    Ensure that the container for the given path exists in the context.
    If it does not exist, create it as an empty dictionary.
    """
    container = context
    for part in path:
        if part not in container:
            container[part] = [{}]
        container = container[part][-1]
    return container

def handleMapping(context, path, value):
    pathEnd = path[-1]
    container = ensureContainerExists( context, path[:-1] )
    if pathEnd in container: # Check if the key already exists in the container.
        # This is possibly an array item if the item already exists.  Convert it if not already an array.
        if not isinstance(container.get(pathEnd), list):
            oldItem = container.pop(pathEnd)
            container[pathEnd] = [ oldItem ]
        container.get( pathEnd ).append( value )
    else:
        container[pathEnd] = value

def handleSimpleMapping(context, key, value):
    collection = key[0].split('.')[-3]
    if collection not in context:
        context[collection] = {}
    collectionMap = context[collection]
    collectionMap.update( {key[0].split('.')[-1]:  value})
    context[collection] = collectionMap