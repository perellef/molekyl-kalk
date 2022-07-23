
from register import Register

import re

class Isotopformel:
    
    @classmethod
    def lag(cls,formel):

        if isinstance(formel,str):
            isotopformel = {a+b: int(c) for a,b,c in [re.split('(\D+)',atom) for atom in formel.split(" ")]}

        elif isinstance(formel,dict):
            isotopformel = formel

        elif isinstance(formel,list):
            isotopformel = {a+b: int(c) for a,b,c in [re.split('(\D+)',atom) for atom in formel]}

        else:
            raise TypeError

        for isotop in isotopformel:
            if isotop not in Register.hent_isotoper():
                raise AttributeError(f"Isotopen '{isotop}' finnes ikke i registeret!")

        return cls(isotopformel)
        
    def __init__(self,isotopformel):
        self._isotopformel = isotopformel

    def __str__(self):
        output = "isotopformel: "
        for isotop,antall in self._isotopformel.items():
            output += f"{isotop}{antall} "
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