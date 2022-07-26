
from register import Register
from isotopformel import Isotopformel
from number_styles import subscript

import re

def kombs(n,m,del_komb=[]): # n = sum av isotopene, m = antall unike isotoper

    if m==1:
        yield del_komb+[n-sum(del_komb)]
    if m<=1:
        return

    for i in range(n-sum(del_komb)+1):
        yield from kombs(n,m-1,del_komb + [i])

class Molekylformel:

    @classmethod
    def lag(cls,formel):

        if isinstance(formel,str):
            molekylformel = dict([re.findall('\d+|\D+',atom) for atom in formel.split(" ")])

        elif isinstance(formel,dict):
            molekylformel = formel

        elif isinstance(formel,list):
            molekylformel = dict([re.findall('\d+|\D+',atom) for atom in formel])

        else:
            raise TypeError("Molekylformelen har feil format.")

        molekylformel = {grunnstoff: int(antall) for grunnstoff,antall in molekylformel.items()}
        
        for grunnstoff in molekylformel:
            if grunnstoff not in Register.hent_grunnstoffer():
                raise AttributeError(f"Grunnstoffet '{grunnstoff}' finnes ikke i registeret!")

        return cls(molekylformel)
        
    @classmethod
    def fra_isotopformel(cls,isotopformel):
        formel = isotopformel.hent_formel()

        molekylformel = {}
        for isotop,antall in formel.items():

            grunnstoff = re.findall('\d+|\D+',isotop)[1]

            try:
                molekylformel[grunnstoff] += antall
            except KeyError:
                molekylformel[grunnstoff] = antall

        return cls.lag(molekylformel)


    def __init__(self,molekylformel):
        self._molekylformel = molekylformel
        self._isotopfordeling = None

    def trenger_isotopfordeling(f) :
        def indre(self) :
            if self._isotopfordeling==None:
                self.__finn_isotopfordeling()
            f(self)
        return indre

    def __str__(self):

        grunnstoffer = ["C","H","N","O","S","Cl","Br","F","I","Si","P"]
        sortert_formel = {k:v for x in grunnstoffer for k,v in self._molekylformel.items() if k==x}
        
        output = ""
        for grunnstoff,antall in self._molekylformel.items():
            output += f"{grunnstoff}{subscript(antall)}"
        return output

    def __finn_isotopfordeling(self):
        isotopfordeling = [{}]

        for grunnstoff,n in self._molekylformel.items():
            isotopliste = Register.grunnstoffer[grunnstoff]

            kombinasjoner = [el for el in kombs(n,len(isotopliste))]

            temp = []
            for komb in kombinasjoner:

                isotopformel = {}
                for i,nukl_tall in enumerate(isotopliste):

                    if komb[i]==0:
                        continue

                    isotopformel[nukl_tall+grunnstoff] = komb[i]

                for isotopdel in isotopfordeling:
                    temp.append(isotopformel | isotopdel)

            isotopfordeling = temp

        self._isotopfordeling = [Isotopformel.lag(formel) for formel in isotopfordeling]

    @trenger_isotopfordeling
    def beregn_masse(self):
        m = 0
        for isotopformel in self._isotopfordeling:
            
            masse = isotopformel.beregn_masse()
            forekomst = isotopformel.beregn_forekomst()
            
            m += forekomst*masse
        return m

    @trenger_isotopfordeling
    def print_isotopfordeling(self, n=10):

        fordeling = sorted(self._isotopfordeling, key = lambda x: -x.beregn_forekomst())

        print(f"\nIsotopfordeling til {self}")
        print("------------------------------")
        print(f"De {n}/{len(fordeling)} vanligste isotopsammensetningene:")

        maks_forekomst = max(fordeling, key= lambda x: x.beregn_forekomst()).beregn_forekomst()

        for isotopformel in fordeling[:n]:

            relativ_forekomst = str("{:.6f}".format(100*isotopformel.beregn_forekomst()/maks_forekomst))[:6]
            forekomst = str("{:.6f}".format(100*isotopformel.beregn_forekomst()))[:6]
            masse = "{:.6f}".format(isotopformel.beregn_masse())

            print(f"{forekomst} % | {relativ_forekomst} % | {masse} g/mol | {str(isotopformel)}")

    def kalkuler_DBE(self): # dobbelt-bindingsekvivalenter

        bindinger = Register.hent_bindinger()

        antall_atomer = 0
        antall_bindinger = 0
        for grunnstoff,antall in self._molekylformel.items():

            antall_bindinger += antall*bindinger[grunnstoff]
            antall_atomer += antall
            
        DBE = antall_bindinger/2 - (antall_atomer-1)
        
        return DBE

    @trenger_isotopfordeling
    def hent_isotopfordeling(self):
        return self._isotopfordeling

    