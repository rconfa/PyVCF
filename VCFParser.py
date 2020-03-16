"""
    TODO: file compresso
    'fileDate', 'reference'


    #QUESTION
    # contig ... come è formato perchè non si capisce


"""

import re # import per l'utilizzo delle regex
import sys # import per chiudere il programma
from urllib.parse import urlparse  # import per controllare l'url
from pathlib import Path # import per prendere il path del file nell'url
from dateutil.parser import parse # per controllare se una stringa è una data



class VCFParser:
    # parametri da inizializzare 
    filename = ""
    zipped = True
    
    # Definisco i pattern per le regex comuni come privati 
    __general_match = "[a-zA-Z0-9_]+"
    __id_pattern = "ID=" + __general_match
    __desc_pattern = "Description=\"[^\"]+\""
    __num_patter = "Number=(?P<num>([0-9]+|A|R|G|.))"
    
    # mi salvo le righe in cui trovo fileformat (deve essere la prima),
    # fileDate, reference, pedigree/pedigreeDB
    __lineFileFormat = None
    __lineFileDate = None
    __lineReference = None
    __linePedigree = None
    
    
    
    __actualLine = 0 # riga in cui sono arrivato a leggere (per print error)
    
    # inizializzo la classe, il parametro zipped non è richiesto, di defaul è False
    def __init__(self, filename, zipped = False):
        self.filename = filename
        self.zipped = zipped


    """
        funzione che legge tutto il contenuto di un file e richiama il parser   
        filename = nome del file
        zipped = specifica se il file è in formato compresso oppure no
    """
    def parseFile(self):
        # TODO: gestione della compressione
        
        # apro e leggo tutto il contenuto del file
        with open(self.filename, 'r') as vcf_file:
            lines = vcf_file.readlines()
    
        vcf_file.close() # chiudo il file
        self.__parser(lines) # richiamo il parser
    
    
    # parso tutto il file riga per riga! se incontro un errore il programma
    # stampa il messaggio e si blocca
    def __parser(self, allLines):
        for line in allLines:
            metadataAllows = True # variabile booleana che mi dice se sono ancora nella sezione ## oppure sono nella sezione dati
            
            # guardo che tipo di riga sto leggendo metadati o no
            isMetadata = re.match("##", line) 
            # Se ho trovato una stringa di metadati e sono nella sezione metadati allora:
            # tolgo i ## iniziali e passo la stringa al parser ad hoc per i metadati
            if (isMetadata and metadataAllows):
                self.__parseMetadata(line[2:])
                self.__actualLine = self.__actualLine + 1

        print("GREAT! VCF IS WELL FORMED!!!")


    # Parsa la riga di metadati in ingresso
    def __parseMetadata(self, line):
        """
            1) cerco il primo = dato che ogni stringa è una coppia key-value
            2) faccio uno "switch" per capire in quale metadato mi trovo
        """
        indexEqual = line.find("=") 
        
        info = line[:indexEqual]
        # metto la riga in maiuscolo, cosi matcho anche con le righe in minuscolo
        info = info.upper()
        if (info == "FILEFORMAT"):
            self.__myFileformat(line)
        if (info == "FILEDATE"):
            self.__myFiledate(line)
        elif (info == "FILTER"):
            self.__myFilter(line)
        elif (info == "INFO"):
            self.__myInfo(line)
        elif (info == "FORMAT"):
            self.__myFormat(line)
        elif (info == "ALT"):
            self.__myAlt(line)
        elif (info == "ASSEMBLY"):
            self.__myAssembly(line)
        elif (info == "CONTIG"):
            self.__myContig(line)
        elif (info == "SAMPLE"):
            self.__mySample(line)
        elif (info == "PEDIGREEDB"):
            self.__myPedigreedb(line)
        elif (info == "PEDIGREE"):
            self.__myPedigree(line)        
    
    
        
        
        
    """
       funzione che controlla che la prima riga sia del tipo 
       ##fileformat=VCFvNumber.Number
    """
    def __myFileformat(self, line):
        # controllo che sia la prima volta che incontro la riga fileFormat
        if(self.__lineFileFormat is not None):
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", fileformat already define!")
        elif(self.__actualLine != 0):
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", fileformat must be defined at line 0")
        else:
            # controllo che sia well-formed tramite regex
            response =  re.match("fileformat=VCFv[0-9].[0-9]", line, re.IGNORECASE) 
            
            # se well-formed aggiorno __lineFileFormat altrimenti stampo errore ed esco
            if (response):
                self.__lineFileFormat = self.__actualLine
            else:
                sys.exit("Error at line: " + str(self.__actualLine + 1) + ", fileformat malformed!")
        
        
    """
       controllo la riga FILEDATI 
       ##filedate=date
    """
    def __myFiledate(self, line):
        if(self.__lineFileDate is not None):
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", filedate already define!")
            
        # prendo tutto quello che c'è dopo = e controllo sia una data
        indexEqual = line.find("=") 
        
        date = line[indexEqual+1:]

        result = self.__is_date(date)
        if(result):
            self.__lineFileDate = self.__actualLine
                
       
    """
       controllo la riga FILTER
       ##filter=<ID=ID,Description="description">
    """ 
    def __myFilter(self,line):
        response = re.match("FILTER=<" + self.__id_pattern + ",\s*" + self.__desc_pattern + ">" ,line, re.IGNORECASE)
        
        # se non è well-formed stampo errore
        if (not response):
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", FILTER malformed!")
    
    
    """
       controllo la riga info
       ##info=<ID=ID,Number=number,Type=type,Description="description",Source="source",Version="version">
    """ 
    def __myInfo(self,line):
        response = re.match("INFO=<" + self.__id_pattern + ",\s*" + self.__num_patter + ",\s*" +
                         "Type=(?P<type>(Integer|Float|Flag|Character|String)),\s*" + 
                         self.__desc_pattern + "(,\s*Source=(\")*[a-zA-Z0-9_\s]+(\")*)*" +
                         "(,\s*Version=(\")*[0-9]+(\")*)*>" ,line, re.IGNORECASE)
        
        # se well-formed stampo eventuale warning flag-num
        if (response):
            # controllo che se type = flag allora number = 0.. E' UNO SHOULD BE, QUINDI NON OBBLIGATORIO
            if (response.group("type") == "Flag" and response.group("num") != "0"):
                print("Warning: number should be zero with Flag type at line: " + str(self.__actualLine + 1))
        else:
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", INFO malformed!")
            
    
    
    """
       controllo la riga format
       ##info=<ID=ID,Number=number,Type=type,Description="description">
    """ 
    def __myFormat(self,line):
        response = re.match("FORMAT=<" + self.__id_pattern + ",\s*" + self.__num_patter + ",\s*" +
                         "Type=(?P<type>(Integer|Float|Character|String)),\s*" + 
                         self.__desc_pattern + ">", line, re.IGNORECASE)
        
        # se non è well-formed stampo errore
        if (not response):
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", FORMAT malformed!")
    
    
    """
       controllo la riga alt
       ##info=<ID=ID;ID,Description="description">
    """ 
    def __myAlt(self,line):
        # aggiungo una seconda parte al patter ID per accettare una colon-separated list
        response = re.match("ALT=<" + self.__id_pattern +  "(:" + self.__general_match + ")*"+ ",\s*" + self.__desc_pattern + ">", line, re.IGNORECASE)
        
        # se non è well-formed stampo errore
        if (not response):
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", ALT malformed!")
    
    
    """
       controllo la riga assembly
       ##assembly=url to fasta file
    """ 
    def __myAssembly(self,line):
        
        # sono sicuro che inizia con assembly= perchè lo matcho nella chiamata di funzione
        # quindi prendo tutto quello che c'è dopo l'= e controllo sia una url
        indexEqual = line.find("=") 
        url = line[indexEqual+1:] # +1 per escludere l'uguale
        
        # uso il metodo urlparse della libreria urllib.parse 
        getUrlResponse = urlparse(url)
        
        # prendo il nome del file, voglio che controllare che sia un fasta
        urlFile = Path(getUrlResponse.path).name
        
        # accetto qualsiasi cosa che termini con .fasta, sono sicuro che
        # va bene perchè ho già controllato che sia un url
        response = re.match("(.*).fasta", urlFile, re.IGNORECASE)
        
        # se non è well-formed stampo errore
        if (not response):
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", ASSEMBLY malformed!")
    
    
    """
       controllo la riga contig
       ##contig=<ID=ctg1,URL=ftp://somewhere.org/assembly.fa,...>
    """ 
    def __myContig(self,line):
        response = re.match("CONTIG=<" + self.__id_pattern + ",\s*" + "length=[0-9]*" 
                             + "(,\s* (.*))*",line, re.IGNORECASE)
        
        # se non è well-formed stampo errore
        if (not response):
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", CONTIG malformed!")
    
    
    
    """
       controllo la riga sample
       ##SAMPLE=<ID=S_ID,Genomes=G1_ID;G2_ID; ...;GK_ID,Mixture=N1;N2; ...;NK,Description=S1;S2; ...;SK>
       Return 0 = va bene
       -1 = errore nel formato
    """ 
    def __mySample(self,line):
        response = re.match("SAMPLE=<" + self.__id_pattern + ",\s*" 
                             + "Genomes=" + self.__general_match + "(;" + self.__general_match + ")*" + ",\s*" 
                             + "Mixture=(.[0-9]|[0-9].)" + "(;(.[0-9]|[0-9].))*" + ",\s*"
                             + "Description=\"[^\"]+" + "(;[^\"]+)*" + "\"" 
                             , line, re.IGNORECASE)
        
        # se non è well-formed stampo errore
        if (not response):
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", SAMPLE malformed!")
    
    
    """
       controllo la riga pedigreeDB
       ##pedigreeDB=<url> oppure
    """ 
    def __myPedigreedb(self,line):
        if(self.__linePedigree is not None):
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", pedigree already define!")
        # controllo se il formato è valido
        response_url = re.match("PEDIGREEDB=<(?P<url>(.*))>", line, re.IGNORECASE)
        
        #se il formato non corrisponde mi fermo
        if(not response_url):
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", PEDEGREEDB malformed!")
        
        # se il formato corrisponde controllo che la url sia valida
        getUrlResponse = response_url.group("url")
        
        
        # se non è well-formed stampo errore
        if (not getUrlResponse):
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", PEDEGREEDB malformed!")
        
        self.__linePedigree = self.__actualLine
        
        
        
        
    """
       controllo la riga pedigree
       ##PEDIGREE=<Name_0=G0-ID,Name_1=G1-ID,...,Name_N=GN-ID>
    """ 
    def __myPedigree(self,line):
        if(self.__linePedigree is not None):
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", pedigree already define!")
            
        # controllo se il formato è valido
        response = re.match("PEDIGREE=<" + self.__general_match + "=" + self.__general_match 
                            + "(,\s*" + self.__general_match + "=" + self.__general_match 
                            +")*>", line, re.IGNORECASE)
        
            
        #se il formato non corrisponde mi fermo
        if(not response):
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", PEDIGREE malformed!")
        else: # se corrisponde controllo che i name non siano ripetuti
            # prendo la stringa compresa tra <>
            indexMin = line.find("<") 
            indexMax = line.find(">")
            text = line[indexMin+1:indexMax]
            
            #splitto tutte le stringhe divise da , per ottenere i sottogruppi
            couple = text.strip().split(',')
            
            name = []
            # di ogni sottogruppo elimino la parte dall'= alla fine
            # e salvo la parte iniziale nella lista name
            [name.append(single[:single.find("=")]) for single in couple]
    
            #creo una lista passando da un set, cosi elimino possibili doppioni
            listNoRepeat = list(set(name)) 
            
            # se le due liste non sono lunghe ugali c'erano dei doppioni nei nomi --> errore
            if(len(name) != len(listNoRepeat)):
                sys.exit("Error at line: " + str(self.__actualLine + 1) + ", PEDIGREE malformed!")
        
            self.__linePedigree = self.__actualLine
    
    
    
    def __is_date(self, string, fuzzy=False):
        """
        Return whether the string can be interpreted as a date.
    
        :param string: str, string to check for date
        :param fuzzy: bool, ignore unknown tokens in string if True
        """
        try: 
            parse(string, fuzzy=fuzzy)
            return True
    
        except ValueError:
            return False
        
        
                