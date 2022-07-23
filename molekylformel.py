
from register import Register
from isotopformel import Isotopformel

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
        
    def __init__(self,molekylformel):
        self._molekylformel = molekylformel
        self._isotopfordeling = self.__finn_isotopfordeling()

    def __str__(self):
        output = "molekylformel: "
        for grunnstoff,antall in self._molekylformel.items():
            output += f"{grunnstoff}{antall} "
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

        return [Isotopformel.lag(formel) for formel in isotopfordeling]

    def beregn_masse(self):
        m = 0
        for isotopformel in self._isotopfordeling:
            
            masse = isotopformel.beregn_masse()
            forekomst = isotopformel.beregn_forekomst()
            
            m += forekomst*masse
        return m

    def hent_isotopformel(self):
        return self._isotopformel

    