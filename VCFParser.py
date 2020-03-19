"""
    TODO: file compresso
    data alt
    
"""

import re # import per l'utilizzo delle regex
import sys # import per chiudere il programma
from urllib.parse import urlparse  # import per controllare l'url
from pathlib import Path # import per prendere il path del file nell'url
from dateutil.parser import parse # per controllare se una stringa è una data
import logging, sys # per visualizzare le stampe di debug
import gzip # per leggere i file gzip


class VCFParser:
    # parametri da inizializzare 
    filename = ""
    zipped = True
    
    # Definisco i pattern per le regex comuni come privati 
    __general_match = "[a-zA-Z0-9_]+"
    __id_pattern = "ID=(?P<id>" + __general_match + ")"
    __desc_pattern = "Description=\"[^\"]+\""
    __num_patter = "Number=(?P<num>([0-9]+|A|R|G|.))"
    
    # mi salvo le righe in cui trovo fileformat (deve essere la prima),
    # fileDate, reference, pedigree/pedigreeDB
    __lineFileFormat = None
    __lineFileDate = None
    __lineReference = None
    __linePedigree = None
    
    __actualLine = 0 # riga in cui sono arrivato a leggere (per print error)
    
    # conto le righe vuote per sistemare il conteggio per le print error
    __emptyLine = 0 
    
    # ultimo valore del campo pos trovato cosi controllo che siano in ordine crescente
    __lastPosValue = -1 
    

    # valori ID dei metadati FILTER trovati
    __lstIdFilterValue = []
    # valori ID dei metadati CONTIG trovati
    __lstIdContigValue = []
    # valori ID dei metadati ALT trovati
    __lstIdAltValue = []
    # valori del campo chrom trovati, mi serve perchè le pos devono essere
    # in ordine crescente SOLO per lo stesso chrom e non in assoluto
    # inoltre possono essere usati nel campo alt
    __lstChromValue = []
    # valori pos del campo dati
    __lstPosValue = []
    
    
    # valori accettati per lo starti id del campo alt (metaline)
    __lstAltValueAccepted = ["DEL", "INS", "DUP", "INV", "CNV", "DUP:TANDEM", "DEL:ME", "INS:ME"]
    
    __lenHeaderDataLine = -1
    
    __FormatDataAllowed = False
    
    # valori del campo info (metaline)
    __dictInfo = {
        'AA': 'String', 'AC': 'Integer', 'AF': 'Float', 'AN': 'Integer',
        'BQ': 'Float', 'CIGAR': 'String', 'DB': 'Flag', 'DP': 'Integer',
        'END': 'Integer', 'H2': 'Flag', 'H3': 'Flag', 'MQ': 'Float',
        'MQ0': 'Integer', 'NS': 'Integer', 'SB': 'String', 'SOMATIC': 'Flag',
        'VALIDATED': 'Flag', '1000G': 'Flag',
    
        # Keys used for structural variants
        'IMPRECISE': 'Flag', 'NOVEL': 'Flag', 'SVTYPE': 'String',
        'SVLEN': 'Integer', 'CIPOS': 'Integer', 'CIEND': 'Integer',
        'HOMLEN': 'Integer', 'HOMSEQ': 'String', 'BKPTID': 'String',
        'MEINFO': 'String', 'METRANS': 'String', 'DGVID': 'String',
        'DBVARID': 'String', 'DBRIPID': 'String', 'MATEID': 'String',
        'PARID': 'String', 'EVENT': 'String', 'CILEN': 'Integer',
        'DPADJ': 'Integer', 'CN': 'Integer', 'CNADJ': 'Integer',
        'CICN': 'Integer', 'CICNADJ': 'Integer'
    }
   
    # se nella riga dati, campo alt, ci sono dei breakends formati da
    # chrom:pos li aggiungo e li controllo a fine file perchè devo prima scorrere tutti
    # i chrom e le pos
    # nella lista salvo [(stringa, line), (stringa2, line2)...]
    __lstBreakendsValue = []
    
    
    
    
    """ 
        inizializzo la classe
    """
    def __init__(self, filename):
        self.filename = filename
        self.zipped = filename.endswith("gz")
        # inizializzo le stampe di debug
        # level = INFO -> no print, level= DEBUG -> print
        logging.basicConfig(stream=sys.stdout, level=logging.INFO)
        

    """
        funzione che legge tutto il contenuto di un file e richiama il parser   
        filename = nome del file
        zipped = specifica se il file è in formato compresso oppure no
    """
    def parseFile(self):
        # TODO: gestione della compressione
        
        # apro e leggo tutto il contenuto del file
        if(self.zipped == False):
            with open(self.filename, 'r') as vcf_file:
                lines = vcf_file.readlines()
        else:
            lines = []
            with gzip.open(self.filename, "rt") as vcf_file:
                for line in vcf_file:
                    lines.append(line)
                
        vcf_file.close() # chiudo il file
        self.__parser(lines) # richiamo il parser
    
    
    """
        parso tutto il file riga per riga! se incontro un errore il programma
        stampa il messaggio e si blocca
    """
    def __parser(self, allLines):
        metadataAllows = True # variabile booleana che mi dice se sono ancora nella sezione ## oppure sono nella sezione dati
        
            
        for line in allLines:
            if(line.strip() == ""):
                self.__emptyLine = self.__emptyLine +1
            # se sono nella sezione metadati controllo che la stringa inizi con ##
            elif(metadataAllows):
                # guardo che tipo di riga sto leggendo metadati o no
                isMetadata = line.startswith('##')
                 
                # Se ho trovato una stringa di metadati
                # tolgo i ## iniziali e passo la stringa al parser ad hoc per i metadati
                if (isMetadata == True):
                    self.__parseMetadata(line[2:])
                    self.__actualLine = self.__actualLine + 1 + self.__emptyLine
                    self.__emptyLine = 0 # risetto le righe vuote perchè gia sommate
                # altrimenti controllo che sia l'header della sezione dati
                else: 
                    metadataAllows = False
                    line = self.__white_cleaner(line)
                    # richiamo la funzione
                    # se è tutto a posto proseguo altrimenti la funzione
                    # bloccherà il programma
                    self.__myHeaderMetadata(line) 
                    self.__actualLine = self.__actualLine + 1 + self.__emptyLine
                    self.__emptyLine = 0 # risetto le righe vuote perchè gia sommate
            else:
                # ripulisco la linea da tab, spazi etc..
                clearLine = self.__white_cleaner(line)
                lst = clearLine.split(" ")
                self.__parseData(lst)
                self.__actualLine = self.__actualLine + 1 + self.__emptyLine
                self.__emptyLine = 0 # risetto le righe vuote perchè gia sommate
         
        # una volta lette tutte le righe controllo gli eventuali valori breakends di alt
        self.__matchBeetwenParenthesis()       
        print("GREAT! VCF IS WELL FORMED!!!")

    
    """
        parser per le righe di metadati
    """
    def __parseMetadata(self, line):
        """
            1) cerco il primo = dato che ogni stringa è una coppia key-value
            2) faccio uno "switch" per capire in quale metadato mi trovo
        """
        indexEqual = line.find("=") 
        
        # se non trovo = la stringa è mal formata
        if(indexEqual == -1):
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", invalid metadata line. Missing = character!")
        
        info = line[:indexEqual]
        # metto la riga in maiuscolo, cosi matcho anche con le righe in minuscolo
        info = info.upper()
        if (info == "FILEFORMAT"):
            self.__myFileformat(line)
        elif (info == "FILEDATE"):
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
        elif (info == "REFERENCE"):
            self.__myReference(line)
        elif (info == "SOURCE"):
            self.__mySource(line)
        elif (info == "PHASING"):
            self.__myPhasing(line)
        else:
            print("Warning: unknown metadata line starting with: " + info + ", at line: " 
                  + str(self.__actualLine + 1))
            
        
    """
       funzione che controlla che la prima riga sia del tipo 
       ##fileformat=VCFvNumber.Number
    """
    def __myFileformat(self, line):
        # debug printing
        logging.debug("   __myFileformat() parsing line: " + line)
        
        
        # controllo che sia la prima volta che incontro la riga fileFormat
        if(self.__lineFileFormat is not None):
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", fileformat already define!")
        elif(self.__actualLine != 0):
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", fileformat must be defined at line 0")
        else:
            # controllo che sia well-formed tramite regex
            response =  re.match("fileformat=VCFv[0-9].[0-9]$", line, re.IGNORECASE) 
            
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
        # debug printing
        logging.debug("   __myFiledate() parsing line: " + line)
        
        
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
        # debug printing
        logging.debug("   __myFilter() parsing line: " + line)
                
        response = re.match("FILTER=<" + self.__id_pattern + ",\s*" + self.__desc_pattern + ">$" ,line, re.IGNORECASE)
        
        # se non è well-formed stampo errore
        if (not response):
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", FILTER malformed!")
        # se è well formed mi salvo l'id in quanto devo confrontarlo con quello dei campi data
        else:
            self.__lstIdFilterValue.append(response.group("id").upper())
    
    """
       controllo la riga info
       ##info=<ID=ID,Number=number,Type=type,Description="description",Source="source",Version="version">
    """ 
    def __myInfo(self,line):
        # debug printing
        logging.debug("   __myInfo() parsing line: " + line)
        
        response = re.match("INFO=<" + self.__id_pattern + ",\s*" + self.__num_patter + ",\s*" +
                         "Type=(?P<type>(Integer|Float|Flag|Character|String)),\s*" + 
                         self.__desc_pattern + "(,\s*Source=(\")*[a-zA-Z0-9_\s]+(\")*)*" +
                         "(,\s*Version=(\")*[0-9]+(\")*)*>$" ,line, re.IGNORECASE)
        
        # se well-formed stampo eventuale warning flag-num
        if (response):
            self.__updateInfoDictionary(response.group("id"), response.group("type"), response.group("num"))
            # controllo che se type = flag allora number = 0.. E' UNO SHOULD BE, QUINDI NON OBBLIGATORIO
            if (response.group("type") == "Flag" and response.group("num") != "0"):
                print("Warning: number should be zero with Flag type at line: " + str(self.__actualLine + 1))
        else:
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", INFO malformed!")
    
              
    """ 
        controlla se il campo id è nel dizionario, se c'è già controlla che 
        anche il type sia definito allo stesso modo!
        Accetto info line multiple definite dall'utente
        NEL DIZIONARIO SALVO ID = [TIPO, NUMERO]
        in quanto poi nelle righe data mi serve sapere quanti valori possono essere
        aggiunti per ciascun id!
    """
    def __updateInfoDictionary(self,id,Type, num):
        # controllo se è già definito e se il tipo è lo stesso di quello già presente
        if (id in self.__dictInfo):
            if (Type != self.__dictInfo[id]): 
                sys.exit("Error at line: " + str(self.__actualLine + 1) + ", " + id + " id already define. Check type!")
            self.__dictInfo[id] = [Type, num]
            
        else: # se non è definito lo aggiungo come tipo
            self.__dictInfo[id] = [Type, num]
        
        
    """
       controllo la riga format
       ##info=<ID=ID,Number=number,Type=type,Description="description">
    """ 
    def __myFormat(self,line):
        # debug printing
        logging.debug("   __myFormat() parsing line: " + line)       

        response = re.match("FORMAT=<" + self.__id_pattern + ",\s*" + self.__num_patter + ",\s*" +
                         "Type=(?P<type>(Integer|Float|Character|String)),\s*" + 
                         self.__desc_pattern + ">$", line, re.IGNORECASE)
        
        # se non è well-formed stampo errore
        if (not response):
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", FORMAT malformed!")
            
        self.__FormatDataAllowed = True
    
    """
       controllo la riga alt
       ##info=<ID=ID;ID,Description="description">
    """ 
    def __myAlt(self,line):
        # debug printing
        logging.debug("   __myAlt() parsing line: " + line)   
        
        # aggiungo una seconda parte al patter ID per accettare una colon-separated list
        response = re.match("ALT=<ID=(?P<id>" + self.__general_match  +  "(:" + self.__general_match + ")*)"+ ",\s*" + self.__desc_pattern + ">$", line, re.IGNORECASE)
        
        # se non è well-formed stampo errore
        if (not response):
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", ALT malformed!")
    
        """
            la regex ha dato positive match, controllo che il campo id sia valido
            deve iniziare con qualcosa presente nella lista __lstAltValueAccepted
        """
        startOk = False
        index = 0
        idString = response.group("id")
        while startOk == False and index < len(self.__lstAltValueAccepted):
            stringStart = self.__lstAltValueAccepted[index]
            # controllo se l'inizio matcha con una delle possibilità
            if(idString.startswith(stringStart)):
                startOk = True
            index = index + 1 # passo alla stringa successiva
            
            
        # se non inizia con uno dei caratteri consentiti mi fermo
        if (startOk == False):
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", ALT ID malformed!")
        
        # se l'id è valido ed è composto da CNV stampo un warning
        if (idString.startswith("CNV")):
            print("Warning! At line: " +  str(self.__actualLine + 1) + 
                  " CNV category should not be used when a more specific category can be applied")
        self.__lstIdAltValue.append("<" + response.group("id") +">")    
    
    """
       controllo la riga assembly
       ##assembly=url to fasta file
    """ 
    def __myAssembly(self,line):
        # debug printing
        logging.debug("   __myAssembly() parsing line: " + line) 
        
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
        response = re.match("(.*).(fasta|fa)$", urlFile, re.IGNORECASE)
        
        # se non è well-formed stampo errore
        if (not response):
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", ASSEMBLY malformed!")
    
    
    """
       controllo la riga contig
       ##contig=<ID=ctg1,lenght=number,....>
    """ 
    def __myContig(self,line):
        # debug printing
        logging.debug("   __myContig() parsing line: " + line)        
        
        response = re.match("CONTIG=<" + self.__id_pattern + ",\s*" + "length=[0-9]+" 
                             + "(,(\s)*(.*))*>$", line, re.IGNORECASE)
        
        # se non è well-formed stampo errore
        if (not response):
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", CONTIG malformed!")
        
        # aggiungo l'id alla lista, mi servirà per i campi data
        self.__lstIdContigValue.append(response.group("id").upper())
    
    
    """
       controllo la riga sample
       ##SAMPLE=<ID=S_ID,Genomes=G1_ID;G2_ID; ...;GK_ID,Mixture=N1;N2; ...;NK,Description=S1;S2; ...;SK>
    """ 
    def __mySample(self,line):
        # debug printing
        logging.debug("   __mySample() parsing line: " + line)      
        
        response = re.match("SAMPLE=<" + self.__id_pattern + ",\s*" 
                             + "Genomes=" + self.__general_match + "(;" + self.__general_match + ")*" + ",\s*" 
                             + "Mixture=(.[0-9]|[0-9].)" + "(;(.[0-9]|[0-9].))*" + ",\s*"
                             + self.__desc_pattern +">$"
                             , line, re.IGNORECASE)
        
        # se non è well-formed stampo errore
        if (not response):
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", SAMPLE malformed!")
    
    
    """
       controllo la riga pedigreeDB
       ##pedigreeDB=<url>
    """ 
    def __myPedigreedb(self,line):
        # debug printing
        logging.debug("  __myPedigreedb() parsing line: " + line)           
        
        if(self.__linePedigree is not None):
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", pedigree already define!")
        # controllo se il formato è valido
        response_url = re.match("PEDIGREEDB=<(?P<url>(.*))>$", line, re.IGNORECASE)
        
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
        # debug printing
        logging.debug("   __myPedigree() parsing line: " + line)        
        
        if(self.__linePedigree is not None):
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", pedigree already define!")
            
        # controllo se il formato è valido
        response = re.match("PEDIGREE=<" + self.__general_match + "=" + self.__general_match 
                            + "(,\s*" + self.__general_match + "=" + self.__general_match 
                            +")*>$", line, re.IGNORECASE)
        
            
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
    
    
    """
        Return whether the string can be interpreted as a date.
    
        :param string: str, string to check for date
        :param fuzzy: bool, ignore unknown tokens in string if True
    """
    def __is_date(self, string, fuzzy=False):
        try: 
            parse(string, fuzzy=fuzzy)
            return True
    
        except ValueError:
            return False
    
    
    """
       controllo la riga reference
       ##REFERENCE=...
    """ 
    def __myReference(self,line):
        # debug printing
        logging.debug("   __myReference() parsing line: " + line)      
        
        if(self.__lineReference is not None):
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", reference already define!")
            
            
        #response = re.match("reference=file:///[\w|/|-]+.(fa|fasta)$", line, re.IGNORECASE)
        response = re.match("reference=[\w:/._-]+$", line, re.IGNORECASE);
        # se non è well-formed stampo errore
        if (not response):
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", REFERENCE malformed!")


    """
       controllo la riga source
       ##SOURCE=...
    """ 
    def __mySource(self,line):
         # debug printing
        logging.debug("   __mySource() parsing line: " + line)      
          
        #response = re.match("reference=file:///[\w|/|-]+.(fa|fasta)$", line, re.IGNORECASE)
        response = re.match("source=[\w:/._-]+$", line, re.IGNORECASE);
        # se non è well-formed stampo errore
        if (not response):
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", SOURCE malformed!")
        
        
    """
       controllo la riga phasing
       ##PHASING=...
    """ 
    def __myPhasing(self,line):
         # debug printing
        logging.debug("   __myPhasing() phasing line: " + line)      
          
        #response = re.match("reference=file:///[\w|/|-]+.(fa|fasta)$", line, re.IGNORECASE)
        response = re.match("phasing=[\w:/._-]+$", line, re.IGNORECASE);
        # se non è well-formed stampo errore
        if (not response):
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", PHASING malformed!")
            
            
            
    """
        controllo se la riga è l'intestazione delle righe dati
        #CHROM POS ID REF ALT QUAL FILTER INFO ........
    """    
    def __myHeaderMetadata(self,line):
         # debug printing
        logging.debug("   __myHeaderMetadata() parsing line: " + line)      
        
        checkIfHeader = line.startswith('#CHROM')
        
        if(checkIfHeader == False):
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ",missing HEADER DATALINE!")
            
        if(self.__FormatDataAllowed == True):
            regexFormatIfAllowed = "(((FORMAT))((\s)+([a-zA-Z0-9_])+)*){0,1}"
        else:
            regexFormatIfAllowed = ""
        response = re.match("#CHROM(\s)+POS(\s)+ID(\s)+REF(\s)+ALT(\s)+"
                            +"QUAL(\s)+FILTER(\s)+INFO(\s)+"
                            + regexFormatIfAllowed + "$"
                            , line, re.IGNORECASE)
        
        # se non è well-formed stampo errore, faccio piu controlli per stampe piu chiare
        if (not response and regexFormatIfAllowed == "" and "FORMAT" in line):
            sys.exit("Error at line: " + str(self.__actualLine + 1) 
                     + ", HEADER DATALINE malformed! Too add FORMAT field you must add format metaline")
        elif(not response):
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", HEADER DATALINE malformed!")
        
        self.__lenHeaderDataLine = len(line.split(" "))
        
        
    
    """ 
        funzione che serve per sostituire tutti i doppi spazi/tab con un singolo
        spazio
    """
    def __white_cleaner(self,line):
        line = line.strip()
        line = line.replace('\t', ' ') # sostituisco i tab con gli spazi
        
        # sostituisco i doppi spazi, ciclo per sistemare anche se ci sono 3 o piu spazi
        while '  ' in line:
            line = line.replace('  ', ' ')
            
            
        return line
    
    
    """
        prende in input una lista contenente i vari campi data
        gia splittati!
        Ricorda: i primi campi sono fissi, posso richiamare direttamente
        la funzione corretta:
        CHROM POS ID REF ALT QUAL FILTER INFO
    """
    def __parseData(self,lst):
        if (self.__lenHeaderDataLine != len(lst)):
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", missing same data!")
        
        self.__myChrom(lst[0])
        self.__myPos(lst[1])
        self.__myId(lst[2])
        self.__myRef(lst[3])
        numDataAllowed = self.__myDataAlt(lst[4])
        self.__myQual(lst[5])
        self.__myDataFilter(lst[6])
        self.__myDataInfo(lst[7],numDataAllowed)
        
    """
        controlla che il campo chrome della riga dati sia well-formed e valido
        
    """   
    def __myChrom(self,chrom):
        # debug printing
        logging.debug("   __myChrom() parsing chrom: " + chrom)
        
        response = re.match('(([a-z0-9_]+$)|(.$)|(<([a-z0-9]+|.)>$)|.$)', chrom, re.IGNORECASE)
        
        # se non è well-formed stampo errore
        if (not response):
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", CHROM malformed!")
        
        # se è il primo chrom che incontro lo aggiungo alla lista e parto con l'ordine crescente per le pos
        if(len(self.__lstChromValue) == 0):
            self.__lstChromValue.append(chrom)
            self.__lastPosValue = -1
        # se trovo un nuovo chrom parto di nuovo con l'ordine crescente
        elif(self.__lstChromValue[-1] != chrom):
            self.__lstChromValue.append(chrom)
            self.__lastPosValue = -1  
        
    """
        controlla che il campo pos della riga dati sia well-formed e valido
        I campi pos devono essere in ordine crescente, quindi confronto il valore
        che trovo con il precedente
    """      
    def __myPos(self,pos):
        # debug printing
        logging.debug("   __myPos() parsing pos: " + pos)
        
        response = re.match("([0-9]+|.)$", pos, re.IGNORECASE)
        
        # se non è well-formed stampo errore
        if (not response):
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", POS malformed!")
        # controllo se pos ha un valore o non è noto
        elif(pos != "."):  
            # se il pos è noto allora controllo sia maggiore rispetto al precedente
            # e ne aggiorno il valore, altrimenti esco.
            if (int(pos) < self.__lastPosValue):
                sys.exit("Error at line: " + str(self.__actualLine + 1) + ", POS value must be in ascending order!")
            self.__lastPosValue = int(pos)
            
            # salvo il valore del pos trovato nella lista perchè mi serve per il campo alt
            self.__lstPosValue.append(pos)
            
            
    """
        controlla che il campo id della riga dati sia well-formed e valido
        
    """   
    def __myId(self,id):
        # debug printing
        logging.debug("   __myId() parsing id: " + id)
        
        response = re.match('(([a-z0-9]+$)|(.$))', id, re.IGNORECASE)
        
        # se non è well-formed stampo errore
        if (not response):
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", ID malformed!")


    """
        controlla che il campo ref della riga dati sia well-formed e valido
        
    """   
    def __myRef(self,ref):
        # debug printing
        logging.debug("   __myRef() parsing ref: " + ref)
        
        response = re.match('((A|C|G|T|N)+)|.$', ref, re.IGNORECASE)
        
        # se non è well-formed stampo errore
        if (not response):
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", REF malformed!")           


    """
        controlla che il campo alt della riga dati sia well-formed e valido
        Alt può essere formato da:
        1) una stringa formata da A,C,G,T,N,*
        2) . se valore sconosciuto
        3) <id_alt>
        4) un angle brackets string, vedesi specifiche pdf page 12
        Nei casi 1-4 può essere formata da piu valori suddivisi da ,
    """   
    def __myDataAlt(self,alt):
        # debug printing
        logging.debug("   __myDataAlt() parsing dataAlt: " + alt)
        
        numberOfDataAllowed = 1
        # controllo se sono nel caso due
        if(alt.startswith("<")):
            if(alt in self.__lstIdAltValue == False):
                sys.exit("Error at line: " + str(self.__actualLine + 1) 
                        + ", ALT malformed! No value id defineds")
        elif(alt == "."): # case 2
            return numberOfDataAllowed
        else: # case 1 and 4
            # splitto la stringa ad ogni , e poi confronto 
            lstSplitted = alt.split(",")
            numberOfDataAllowed = len(lstSplitted)
            for val in lstSplitted:
                # controllo se è nel caso 1
                response = self.__matchAltBaseString(val)
                
                # se il caso 1 è falso controllo il 4, se è falso error
                if (not response):
                    """
                        controllo se la stringa parte con [ cerco l'altra [, 
                        matcho tutto dopo [ con __matchAltBaseString
                        e quello tra [] controllo sia chrom:pos
                    """
                    beetweenPar = None
                    baseStr = None
                    if(val.startswith("[")):
                        indexPar = val[1:].find("[") +1
                        beetweenPar = val[1:indexPar]
                        baseStr = val[indexPar+1:]
                    elif(val.startswith("]")):
                        indexPar = val[1:].find("]") +1
                        beetweenPar = val[1:indexPar]
                        baseStr = val[indexPar+1:]   
                    elif(val.endswith("]")):
                        indexPar = val.find("]")
                        beetweenPar = val[indexPar+1:len(val)-1]
                        baseStr = val[:indexPar]
                    elif(val.endswith("[")):
                        indexPar = val.find("[")
                        beetweenPar = val[indexPar+1:len(val)-1]
                        baseStr = val[:indexPar]   
                        
                    if (beetweenPar == None or baseStr == None):
                        sys.exit("Error at line: " + str(self.__actualLine + 1) + ", ALT malformed!")
                    else:
                        # matcho tutto quello dopo con baseStr
                        res = self.__matchAltBaseString(baseStr)
                        
                        # creo la coppia e l'aggiungo alla lista
                        couple = (beetweenPar, self.__actualLine)
                        self.__lstBreakendsValue.append(couple)

                        if(res == False):
                            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", ALT malformed!")               
        
        
        
        return numberOfDataAllowed # mi serve sapere da quanti valori è composto il campo alt per gestire il campo data
    
    
    
    """
        prova a matchare la stringa con A,c,g,t,n,.,*
    """     
    def __matchAltBaseString(self, baseStr): 
        response = re.match('[A|C|G|T|N|.|\*]+$', baseStr, re.IGNORECASE)
        
        if(response):
            return True
        
        return False
        
    """
        controllo che la __lstBreakendsValue sia del tipo chromVal:PosVal
        con valori di chrom e pos validi
    """
    def __matchBeetwenParenthesis(self):
        for couple in self.__lstBreakendsValue:
            values = couple[0].split(":") 
       
            if(len(values)!=2): # non ho trovato i : allora non è ben formata
                sys.exit("Error at line: " + str(int(couple[1]) + 1) + ", ALT breakends malformed!")
            # controllo che i due valori siano chrom:pos
            elif(values[0] not in self.__lstChromValue or values[1] not in self.__lstPosValue):
                sys.exit("Error at line: " + str(int(couple[1]) + 1) + ", ALT breakends value not valid!")
    
    """
        controlla che il campo qual della riga dati sia well-formed e valido
        
    """   
    def __myQual(self,qual):
        # debug printing
        logging.debug("   __myQual() parsing qual: " + qual)
        
        response = re.match('([0-9]+|[0-9]+.[0-9]+|.)$', qual, re.IGNORECASE)
        
        # se non è well-formed stampo errore
        if (not response):
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", qual malformed!")               
      
      
    """
        controlla che il campo filter della riga dati sia well-formed e valido
        In questo caso non serve una regex il campo puo essere:
        PASS
        .
        una combinazione dei valori presenti in __lstIdFilterValue
    """   
    def __myDataFilter(self,filter):
        # debug printing
        logging.debug("   __myDataFilter() parsing filter: " + filter)
        
        filter = filter.upper() # metto tutto in maiuscolo per semplicita
        
        # controllo se uguale a PASS
        if(filter == "PASS"):
            return
            
        # controllo se uguale al .
        if(filter == "."):
            return
            
        # controllo se uguale a 0 (meglio non usarlo ma possibile)
        if(filter == "0"):
            print("Warning! Value 0 should not be use as value for filter at line: " + str(self.__actualLine + 1))
            return
        
        # cerco il primo ;
        indexEqual = filter.find(";") 
        
        # se non trovo nessun ; controllo che sia un singolo valore presente nella lista
        if(indexEqual == -1 and filter in self.__lstIdFilterValue):
            return
        else:
            # ci sono dei punti e virgola, controllo che il numero di ; 
            # non superi il numero di filtri presenti nella lista
            # es: lst = q10,s50    filter=q10;s50;q10 NO!
            count = len(re.findall(";", filter)) 
            
            if(count < len(self.__lstIdFilterValue)):
                values = filter.split(";") # splitto tutti i valori
                # controllo che ogni valore sia in list
                for v in values:
                    if (not v in  self.__lstIdFilterValue):
                        sys.exit("Error at line: " + str(self.__actualLine + 1) 
                                 + ", filter malformed. " + v + " not defined as a filter")  
                return        
            else:
                sys.exit("Error at line: " + str(self.__actualLine + 1) 
                                 + ", filter malformed. Too much filter.")        
                    
        # se arrivo qui nessuno dei controlli è risultato positivo, la stringa è mal formata
        sys.exit("Error at line: " + str(self.__actualLine + 1) + ", filter malformed!")       

            
    """
        controlla il campo info della riga dati ogni campo è formato da
        key=value,value2,..;key2=value;key3
        Ogni chiave può avere 0...n valori in base a come definito nelle righe info di metadati
        per recuperare queste informazioni sfrutto __dictInfo
        Se la chiave è presente controllo quanti valori deve avere, se ne ha
        controllo che siano stati tutti inseriti e che rispettino le richieste.
        altNumber mi serve nel caso in cui il numero sia A,R,G cosi posso controllare
        che abbia il numero di valori corretto!
        R = altNumber + 1
        A = altNumber
        G = ????
    """
    def __myDataInfo(self,info, altNumber):
        # debug printing
        logging.debug("   __myDataInfo() parsing info: " + info)
        
        splittedInfo = info.split(";")
       
        for element in splittedInfo:
            lstRes = element.split("=") # splitto l'=
            # controllo se quella chiave è nel dizionario
            if(lstRes[0] in self.__dictInfo):
                # keyValues è una lista [type, number]
                # type = tipo del valore, number = occorrenze
                keyValues = self.__dictInfo[lstRes[0]] # prendo i valori
                
                # se il number è un numero controllo i suoi possibili valori tra 0..n
                if(keyValues[1].isdigit()):
                    # controllo che se non devono esserci valori ci sia solo il campo id
                    if (keyValues[1] == "0" and len(lstRes)==1): 
                        continue
                    elif(keyValues[1] == "0" and len(lstRes)==2 and lstRes[1]=="0"): 
                        continue             
                    # controllo che abbia il numero di valori corretto
                    # se si controllo il tipo sia corretto per ognuno
                    else:
                        lstValues = lstRes[1].split(",")
                        if(int(keyValues[1]) == len(lstValues)):
                            # richiamo il metodo per il controllo
                            self.__ciclyChecker(lstValues, keyValues[0])
                        else:
                            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", " 
                                    + lstRes[0] + " doesn't have the correct number of value! ")
                # se è uguale al . può avere quanti valori vuole ma devono tutti corrispondere
                elif(keyValues[1] == "."):
                    lstValues = lstRes[1].split(",") 
                    self.__ciclyChecker(lstValues, keyValues[0])
                # se numero = A
                elif(keyValues[1] == "A"):
                    lstValues = lstRes[1].split(",") 
                    if(len(lstValues) == altNumber):
                        # richiamo il metodo per il controllo
                        self.__ciclyChecker(lstValues, keyValues[0])
                    else:
                        sys.exit("Error at line: " + str(self.__actualLine + 1) + ", " 
                            + lstRes[0] + " doesn't have the correct number of value! ")
                # se numero = R
                elif(keyValues[1] == "R"):
                    lstValues = lstRes[1].split(",") 
                    if(len(lstValues) == altNumber+1):
                        # richiamo il metodo per il controllo
                        self.__ciclyChecker(lstValues, keyValues[0])
                    else:
                        sys.exit("Error at line: " + str(self.__actualLine + 1) + ", " 
                            + lstRes[0] + " doesn't have the correct number of value! ")       
                # se numero = G
                elif(keyValues[1] == "G"):
                    # TODO!
                    print("TODO")


    """    
        controlla che la variabile passata in input sia del tipo richiesto!
        
        # NOTA: se il tipo è flag accetto tutto in quanto non è descritto
    """
    def __checkValType(self,val, type):
        
        type = type.upper()
        
        # se il tipo è un intero controllo che la stringa sia di digit, se fallisce
        # chiudo il programma
        if(val != "." and type == "INTEGER" and not val.lstrip('-+').isnumeric()):
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", " 
                    + val + " must be an integer value!")
        if(val != "." and type == "FLOAT" and not self.__isfloat(val)):
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", " 
                    + val + " must be a float value!")          
        # se è una stringa la matcho con una regex
        if(val != "." and type == "STRING"):
            response = re.match('[a-z0-9_]+', val, re.IGNORECASE)
            
            if(not response):
                sys.exit("Error at line: " + str(self.__actualLine + 1) + ", " 
                         + val + " must be a string value!")
        if(val != "." and type == "CHARACTER" and (not val.isalnum() or len(val) > 1)):
            sys.exit("Error at line: " + str(self.__actualLine + 1) + ", " 
                    + val + " must be a single character value!")


    """ 
        true se il valore passato è un float
        Considero float anche i valori singoli es: 3 è un float
    """
    def __isfloat(self, value):
        """
            # se non voglio considerare gli interi come float faccio il controllo
            index = info.find(".") 
            indexComma = info.find(",") 
            if (index == -1 and indexComma == -1):
                return false;
        """
        try:
            float(value)
            return True
        except ValueError:
            return False    
        
        
    """
        Data una lista di nomi di ID del campo INFO nelle riga data
        controlla se questi valori corrispondano al type richiesto
    """
    def __ciclyChecker(self, lstValues, type):    
        for val in lstValues:
            self.__checkValType(val, type)