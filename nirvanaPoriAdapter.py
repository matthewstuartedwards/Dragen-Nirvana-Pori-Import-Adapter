import argparse
from NirvanaJsonAdapter import NirvanaJsonAdapter
from CnvAdapter import CnvAdapter
#from VcfAdapter import VcfAdapter
#from ExpressionAdapter import ExpressionAdapter

    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Nirvana Pori Import Adapter")
    parser.add_argument('--cnv', metavar='c', type=str, required=False, help='Path to the input JSON file with CNV information')
    parser.add_argument('--vcf', metavar='v', type=str, required=False, help='Path to the input JSON file with VCF information')
    parser.add_argument('--diseaseZscores', metavar='z', type=str, required=False, help='Path to the input TSV file with disease Z-scores')
    parser.add_argument('--biopsyZscores', metavar='b', type=str, required=False, help='Path to the input TSV file with biopsy Z-scores')
    parser.add_argument('--outputFile', metavar='o', type=str, required=True, help='Path to the output JSON file' )
    parser.add_argument('--patientID', metavar='p', type=str, required=False, default="ANONYMOUS", help='Patient ID')
    parser.add_argument('--diseaseName', metavar='d', type=str, required=True, help='Disease name for kbDiseaseMatch and is used to populate the matchedCancer flag. eg: sarcoma, colorectal cancer')
    parser.add_argument('--projectName', metavar='j', type=str, required=False, default="NoProjectName", help='Project name for Pori')
    parser.add_argument('--template', metavar='t', type=str, required=False, default="genomic", help='Template for the Pori import. Default is "genomic".')
    args = parser.parse_args()
    
    # Because many objects will be writing to the output file, I'm opening it here.
    # Could refactor this so that each object opens the file and writes when needed with a lock
    # That would allow multithreaded processing until printing is needed.
    with open(args.outputFile, 'w') as output_handle:
        mainAdapter = NirvanaJsonAdapter( output_handle )
        mainAdapter.printOutputHeader(args.patientID, args.diseaseName, args.projectName, args.template)
        
        
        if args.cnv:
            cnvAdapter = CnvAdapter( output_handle, args.cnv )
            cnvAdapter.readCnvFile( args.cnv )
        
        # if args.vcf:
        #     adapter = VcfAdapter(args.vcf)
            
        # if args.diseaseZscores and args.biopsyZscores:
        #     adapter = ExpressionAdapter(args.diseaseZscores, args.biopsyZscores)
            
        
        mainAdapter.printOutputFooter()
    
    
    
    