
from register import Register
from number_styles import superscript,subscript

import re

class Isotopformel:
    
    @classmethod
    def lag(cls,formel):

        if isinstance(formel,str):
            isotopformel = {a+b: int(c) for a,b,c in [re.split('(\D+)',atom) for atom in formel.split(" ")] if int(c)>0}

        elif isinstance(formel,dict):
            isotopformel = formel

        elif isinstance(formel,list):
            isotopformel = {a+b: int(c) for a,b,c in [re.split('(\D+)',atom) for atom in formel] if int(c)>0}

        else:
            raise TypeError

        for isotop in isotopformel:
            if isotop not in Register.hent_isotoper():
                raise AttributeError(f"Isotopen '{isotop}' finnes ikke i registeret!")

        return cls(isotopformel)
        
    def __init__(self,isotopformel):
        self._isotopformel = {k:v for k,v in isotopformel.items() if v>0}

    def __str__(self):

        isotoper = ["12C","13C","1H","2H","14N","15N","16O","17O","18O","32S","33S","34S","35Cl","37Cl","79Br","81Br","19F","127I","28Si","29Si","30Si","31P"]
        sortert_formel = {k:v for x in isotoper for k,v in self._isotopformel.items() if k==x}
            
        output = ""
        for isotop,antall in sortert_formel.items():
            nukl_tall,grunnstoff = re.findall('\d+|\D+',isotop)
            output += f"{superscript(nukl_tall)}{grunnstoff}{subscript(antall)} "
        return output

    def beregn_masse(self):
        isotoper = Register.hent_isotoper()

        s = sum(N*isotoper[isotop]["masse"] for isotop,N in self._isotopformel.items())
        return round(s,6)

    def beregn_forekomst(self):

        isotoper = Register.hent_isotoper()    
        
        forekomst = 1
        atomer = {}
        for isotop,antall in self._isotopformel.items():
            grunnstoff = re.findall('\d+|\D+',isotop)[1]

            try:
                atomer[grunnstoff].append(antall)
            except KeyError:
                atomer[grunnstoff] = [antall]

            forekomst *= isotoper[isotop]["forekomst"]**antall

        for grunnstoff,antall in atomer.items():
            n = 1
            for r in antall:
                for i in range(1,r+1):
                    forekomst *= n/i
                    n += 1        
        return forekomst

    @staticmethod
    def har_antall_av_isotop(formel,isotop,antall):

        try:
            forekomst = formel[isotop]
        except KeyError:
            con = (antall=="0") or (antall=="0+")

            if "-" in antall:
                con1 = (antall[-1]=="-")
                con2 = (int(antall.split("-")[0])==0)

                return any((con1,con2))

            return con

        if antall.isdigit():
            return forekomst==int(antall)
        elif antall[-1]=="-":
            return forekomst<=int(antall[:-1])
        elif antall[-1]=="+":
            return forekomst>=int(antall[:-1])
        elif "-" in antall:
            fra, til = antall.split("-")
            return forekomst>=int(fra) and forekomst<=int(til)
        else:
            raise TypeError("antall har feil format")

    def kalkuler_DBE(self): # dobbeltbindingsekvivalenter
        bindinger = Register.hent_bindinger()

        antall_atomer = 0
        antall_bindinger = 0
        for isotop,antall in self._isotopformel.items():

            grunnstoff = re.findall('\d+|\D+',isotop)[1]
            
            antall_bindinger += antall*bindinger[grunnstoff]
            antall_atomer += antall
            
        DBE = antall_bindinger/2 - (antall_atomer-1)

        return DBE

    def score(self):

        isotoper_score = {"1H": 1, "12C": 5, "14N": 3, "16O": 5, "32S": 4, "35Cl": 5, "79Br": 10, "81Br": 10}

        score = 0
        for isotop,antall in self._isotopformel.items():
            try:
                iso_score = isotoper_score[isotop]
                antall = self._isotopformel[isotop]
            except KeyError:
                continue

            score += iso_score*antall
        return score

    def hent_formel(self):
        return self._isotopformel