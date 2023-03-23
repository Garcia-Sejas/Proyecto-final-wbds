# Proyecto-final-wbds
Proyecto final
A. FINALIDAD 
El surgimiento de una subpoblación de M. Tuberculosis extremadamente resistentes XDR-TB al tratamiento antibioterápico múltiple, es diferente al de las subpoblaciones catalogadas como M.Tuberculosis multidrogo resistente o MDR- TB debido  a su pronóstico de sobrevida pobre sobre los pacientes. Un solo caso de XDR-TB y el estudio de sus contactos deben ser enfocados como una emergencia sanitaria. Ya que el desarrollo de la XDR-TB revela debilitamiento de los servicios asistenciales en el primer nivel de atención. Siendo dos posibles factores de riesgo más fuertemente asociados con la XDR-TB: 1) Fracaso a un régimen antituberculoso que contiene drogas de segunda línea que incluye un inyectable y una fluoroquinolona y 2) Contacto estrecho con un individuo con XDR-TB.
Debido a este problema, es necesario una investigación no solo la identificación de cepas XDR-TB,  y que mejor que conocer dominios funcionales en el subset de proteínas de interés y  predecir secuencias codificantes.
B. DESCRIBIR LOS PASOS PARA SOLUCIONAR EL PROBLEMA (CÓDIGO Y COMENTARIOS).

Para preparar nuestro entorno, necesitamos de las siguientes librerias:
pip3 install matplotlib
pip3 install pandas
pip3 install pycirclize
pip3 install pyrodigal
pip3 install requests
pip3 install seaborn
pip3 install biopython
import credentials
import matplotlib.pyplot as plt
import numpy  as np
import pandas as pd
import pyrodigal
import requests
import seaborn as sns
import subprocess
import sys
from Bio import SeqIO
from Bio import Entrez
from io                 import StringIO
from matplotlib.patches import Patch
from pycirclize         import Circos
from pycirclize.parser  import Gff
from requests.adapters  import HTTPAdapter, Retry
________________________________________
1. Obtención de una secuencia genómica
________________________________________
Nuestro genoma elegido fue seleccionado del sitio web de NCBI - Genome Mycobacterium tuberculosis del  laboratorio Alamos National 
Los Alamos National Lab,  Genome assembly ASM73579v1 y buscaremos genes de producción de antimicrobianos en el cromosoma con número de acceso CP008959.1
accession = "CP042279.1"
genome = Entrez.efetch(db="nucleotide",
                       id=accession,
                       format="gb",
                       rettype="text")
record = SeqIO.read(genome, "genbank")
genome_length = len(record.seq)
________________________________________
2. Predicción por pyrodigal
________________________________________

Con este código encontramos los genes codificantes de proteínas en nuestro genoma 
________________________________________
[ ]
orf_finder = pyrodigal.OrfFinder()
orf_finder.train(bytes(record.seq))
orf_genes  = orf_finder.find_genes(bytes(record.seq))
________________________________________
Con el siguiente código, logramos almacenar las secuencias aminoacídicas de los genes predichos en un nuevo archivo CP008959.1.faa
Para identificar las secuencias aminoacídicas del genoma de Mycobaterium tuberculosis usamos un prefijo como MicoisaXDR________________________________________
aa_file = accession + ".faa"
prefix  = " MicoisaXDR "
with open(aa_file, "w") as orf_gene:
    orf_genes.write_translations(orf_gene,sequence_id=prefix)
________________________________________
Y almacenamos las cordenadas con este nombre CP008959.1.gff
________________________________________
gff_file = accession + ".gff"
prefix  = " MicoisaXDR "
with open(gff_file, "w") as orf_gene:
    orf_genes.write_gff(orf_gene,sequence_id=prefix)
________________________________________
3. Obtención de un set de secuencias de referencia
________________________________________
Descargamos secuencias de UniProt con el objeto uniprot_ref_seqs
uniprot_api_url  = "https://rest.uniprot.org/uniprotkb/stream"
uniprot_api_args = {"compressed" : "false",
                    "format"     : "fasta",
                    "query"      : "(antibiotic biosynthesis) AND (reviewed:true)"}
uniprot_ref_seqs = requests.get(uniprot_api_url,params=uniprot_api_args).text
________________________________________
Y convertimos las secuencias en formato fasta a un archivo (uniprot_sequences.fasta) en nuestro disco duro
________________________________________
uniprot_seqs_file = open("uniprot_sequences.fasta", "wt")
uniprot_seqs_file.write(uniprot_ref_seqs)
uniprot_seqs_file.close()
________________________________________
4. Creación de una base de datos tipo BLAST
________________________________________
A partir de las secuencias aminoacídicas de las predicciones de pyrodigal, creamos una base de datos, para ello descargamos previamente la carpeta acorde a nuestro software:
https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
copiamos el URL de donde se instalo, usarlo en los siguientes pasos:
________________________________________
makeblastdb_path = “/Users/dell/Desktop/WBDS/ncbi-blast-2.13.0+/bin/makeblastdb"
makeblastdb_command = [makeblastdb_path,'-in',aa_file,'-dbtype','prot']
subprocess.call(makeblastdb_command)
________________________________________
5. Obtención de secuencias de interés en el genoma analizado
blastp_path       = "/Users/dell/Desktop/WBDS/ncbi-blast-2.13.0+/bin/blastp"
blastp_out_format = "6 qseqid sseqid qlen slen qstart sstart qend send score evalue length positive"
blastp_out_file   = accession + ".blast.tsv"
blastp_command    = [blastp_path,
                     "-db",          aa_file,
                     "-query",       "uniprot_sequences.fasta",
                     "-evalue",      "1e-9",
                     "-out",         blastp_out_file,
                     "-outfmt",      blastp_out_format,
                     "-num_threads", "12"]
subprocess.call(blastp_command)
________________________________________
6. Examinación de los resultados de la búsqueda tipo BLAST
________________________________________
 Logramos importar una tabla a un dataframe de pandas y asignamos los nombres de las columnas usando la variable blastp_out_format que definimos anteriormente
________________________________________
blastp_column_names = blastp_out_format.split(" ")[1:]
blastp_df = pd.read_csv(blastp_out_file,sep="\t",names=blastp_column_names)
________________________________________
En el genoma de interés encontramos 629 genes potencialmente asociadas con biosíntesis de antibióticos.
________________________________________
candidate_genes=blastp_df["sseqid"].unique().tolist()
len(candidate_genes)
________________________________________
# Visualización preliminar de los datos
________________________________________
Transformamos la información a un formato legible para pyCirclize
________________________________________
# Transformar nuestro archivo gff en un dataframe
________________________________________
gff_columns     = ["chr","source","feature_type","start","end","score","strand","phase","info"]
gff_df          = pd.read_csv(gff_file,sep="\t",comment="#",header=None,names=gff_columns)
gff_df["start"] = gff_df["start"].astype(int)
gff_df["end"]   = gff_df["end"].astype(int)
________________________________________
# Obtención de información adicional del dataframe gff_df
________________________________________

def get_gff_info(info_str):
    out_dict = {}
    info_arr = info_str.split(";")
    for line in info_arr:
        if "=" in line:
            line_arr    = line.split("=")
            field_name  = line_arr[0]
            field_value = line_arr[1]
            out_dict[field_name] = field_value
    return out_dict
________________________________________

gff_df["annotation"] = gff_df["info"].apply(lambda x: get_gff_info(x))
________________________________________
# Filtramos de datos para incluir solamente los genes identificados como asociados a la biosíntesis de antibióticos
________________________________________

gff_df["candidate"] = gff_df["annotation"].apply(lambda x: "include" if x["ID"] in candidate_genes else "exclude")
________________________________________
# Almacenamos en un nuevo archivo gff para que pyCirclize visualice unicamente los genes de interés
________________________________________
candidate_df = gff_df.copy()
candidate_df = candidate_df[candidate_df["candidate"]=="include"][gff_columns]
candidate_df.to_csv("candidates.gff",sep="\t",header=False,index=False)
________________________________________
# Visualización de los datos con pyCirclize
________________________________________
Con el siguiente código construimos un mapa circular donde se encuentran los operones potenciaes en el genoma.
________________________________________

circos = Circos(sectors={accession: genome_length})
circos.text("Streptomyces sp.")
circos_gff = Gff(gff_file="candidates.gff")
sector = circos.get_sector(accession)
sector = circos.sectors[0]
cds_track = sector.add_track((80, 100))
cds_track.axis(fc="#EEEEEE", ec="none")
cds_track.genomic_features(circos_gff.extract_features("CDS", target_strand =  1), r_lim=(90, 100),fc="red" )
cds_track.genomic_features(circos_gff.extract_features("CDS", target_strand = -1), r_lim=(80,  90),fc="blue")
pos_list, labels = [], []
cds_track.xticks_by_interval(
    interval=500000,
    label_formatter=lambda label_value: f"{label_value/ 1000000:.1f} Mb",
    label_orientation="vertical")
fig = circos.plotfig().set_figwidth(5)
________________________________________
# para una mejor visualización de los datos lo analizamos por seaborn
________________________________________
num_bins = 25
counter_1 = 0
counter_2 = 0
fig, axes = plt.subplots(5,5,figsize=(30,30))
bin_len  = (genome_length - (genome_length % (num_bins - 1))) / (num_bins)
for bin_num in range(num_bins):
    start_pos = bin_num * bin_len
    end_pos   = (bin_num + 1) * bin_len
    mb_df = gff_df.copy()
    mb_df = mb_df[(mb_df["start"]>start_pos) & (mb_df["end"]<=end_pos)]
    sns.swarmplot(ax = axes[counter_1,counter_2],data = mb_df,y="candidate",x="start",hue="strand",dodge=True,order=["exclude","include"],hue_order=["+","-"])
    axes[counter_1,counter_2].set(ylabel=None)
    counter_2 += 1
    if (counter_2%5 == 0):
        counter_2 = 0
        counter_1 += 1
plt.show()
________________________________________
7. Examinación detallada del operón seleccionado
________________________________________
De los genes de dicha región a través del servicio de búsqueda de dominios conservados de InterProScan.
________________________________________
operon_df = gff_df.copy()
operon_df = operon_df[(operon_df["start"]     >= 3350000) &
                      (operon_df["end"]       <= 3450000) &
                      (operon_df["strand"]    == "-")     &
                      (operon_df["candidate"] == "include")]
operon_df.reset_index(drop=True, inplace=True)
________________________________________
Obtuvimos 2 operones…
 len(operon_df)
________________________________________
operon_gene_list = []
for index in operon_df.index.tolist():
    gene_id = operon_df["annotation"][index]["ID"]
    operon_gene_list.append(gene_id)
________________________________________
Los cuales  fueron:
['Micoisaxdr_3252', 'Micoisaxdr_3254']
________________________________________
#ENVIO INTERPROSCAN
Eliminamos los asteriscos.

query_str = ""
________________________________________

for record in SeqIO.parse(aa_file, "fasta"):
    seq_id  = record.id
    if(seq_id in operon_gene_list):
        seq_str = str(record.seq)
        query_str+=">"+seq_id+"\n"+seq_str+"\n"
query_str = query_str.replace("*","")
________________________________________
Agregamos los links específicos para cada paso: 
submit_url   = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/run"
progress_url = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/status"
results_url  = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/result"
________________________________________
De igual manera, cada etapa acepta headers específicos
________________________________________
submit_headers   = {"Accept":"text/plain"}
progress_headers = {"Accept":"text/plain"}
results_headers  = {"Accept":"text/tab-separated-values"}
________________________________________
# Enviamos las secuencias
________________________________________
En esta etapa, construiremos un diccionario de python que adjuntaremos a requests para buscar los dominios funcionales
________________________________________
[ ]
submit_data = {"email":"vflorelo@gmail.com",
               "title":"operon_335_345",
               "goterms":"false",
               "pathways":"false",
               "stype":"p",
               "sequence":query_str}
________________________________________
submit_request = requests.post(submit_url,data=submit_data,headers=submit_headers)
________________________________________
submit_status_code = submit_request.status_code
submit_job_id      = submit_request.text
________________________________________
print(submit_status_code)
print(submit_job_id)
________________________________________
# Podemos ver el Progreso de nuestro análisis:
________________________________________
progress_request     = requests.get(progress_url+"/"+submit_job_id,headers=progress_headers)
progress_status_code = progress_request.status_code
progress_status      = progress_request.text
print(progress_status_code)
print(progress_status)
200
RUNNING todo esta OK!!!

# Obtención de los resultados
•	log para acceder al reporte de texto del programa
•	tsv para acceder al reporte tabular programa
________________________________________
results_log_request = requests.get(results_url+"/"+submit_job_id+"/log",headers=results_headers)
results_tsv_request = requests.get(results_url+"/"+submit_job_id+"/tsv",headers=results_headers)
________________________________________

print(results_log_request.text)
________________________________________
# Examinamos los dominios conservados con StringIO
________________________________________
results_tsv_str = StringIO(results_tsv_request.text)
results_column_names = ["sequence","md5","length","database","accession","description","start","end","evalue","post_processed","date","entry","name"]
results_df = pd.read_csv(results_tsv_str,sep="\t",names=results_column_names)
________________________________________
results_df
________________________________________


C. SOFTWARE Y BASES DE DATOS ESPECÍFICOS.
Trabajamos con Jupyter NOtebook desde un Windoms: recurrimos a este genoma en particular ya que no tenia ningua anotación
Annotation details
GenBank
Provider	NCBI
Name	NCBI Prokaryotic Genome Annotation Pipeline
Date	Jul 18, 2014
Genes	3979
Protein-coding	3816
Non-coding	54
D. RESUME LOS RESULTADOS QUE OBTUVISTE EXPLICANDO SI TU APROXIMACIÓN FUE ADECUADA O NO. EN ESTA SECCIÓN DE LA DOCUMENTACIÓN QUE INCLUYA LAS CONCLUSIONES Y OBSERVACIONES QUE PUDISTE HACER GRACIAS A ESTE ANÁLISIS
En resumen, de todas las secuecias aminoacidicas se logro identificar  2 genes asociados a la biosíntesis de antibióticos.
