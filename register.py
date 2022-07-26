import re

class Register:

    with open("isotoper.txt") as txt_file:
        linjer = txt_file.readlines()[1:]

    isotoper = {a: {"masse": float(b), "forekomst": float(c)/100} for a,b,c in [el.strip().split(";") for el in linjer]}

    grunnstoffer = {}
    for isotop in isotoper:
        nukleontall,grunnstoff, = re.findall('\d+|\D+',isotop)
        try:
            grunnstoffer[grunnstoff].append(nukleontall)
        except KeyError:
            grunnstoffer[grunnstoff] = [nukleontall]

    with open("bindinger.txt") as txt_file:
        linjer = txt_file.readlines()[1:]

    bindinger = {a: int(b) for a,b in [el.strip().split(";") for el in linjer]}


    @classmethod
    def hent_isotoper(cls):
        return cls.isotoper

    @classmethod
    def hent_grunnstoffer(cls):
        return cls.grunnstoffer
    
    @classmethod
    def hent_bindinger(cls):
        return cls.bindinger