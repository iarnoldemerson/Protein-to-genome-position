# Protein-to-genome-position
Perl program to convert the protein position into genomic position


1) Required perl modules:

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Transcript;



2) Input format
(five column details : Ensemble_id,	domain_start,	domain_end,	Pfam_id,	Pfam_name)

ENST00000390396.1	20	115	PF07686.12	V-set
ENST00000390400.2	20	114	PF07686.12	V-set
ENST00000621184.1	20	115	PF07686.12	V-set
ENST00000390372.3	20	114	PF07686.12	V-set
ENST00000390369.2	20	115	PF07686.12	V-set




3) How to run the program
(type at the command or terminal prompt)

perl protein_to_genome.pl



4) Output format
(tow column details : Ensemble_id,	Chr_no:start-end)

ENST00000390396.1	7:142646178-142646465	
ENST00000390400.2	7:142720874-142721158	
ENST00000621184.1	7:142581138-142581425	
ENST00000390372.3	7:142482734-142483018	
ENST00000390369.2	7:142455346-142455633	


