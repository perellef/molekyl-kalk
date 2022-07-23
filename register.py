import re

class Register:

    with open("isotoper.txt") as txt_file:
        lines = txt_file.readlines()[1:]

    isotoper = {a: {"masse": float(b), "forekomst": float(c)/100} for a,b,c in [el.strip().split(";") for el in lines]}

    grunnstoffer = {}
    for isotop in isotoper:
        nukleontall,grunnstoff, = re.findall('\d+|\D+',isotop)
        try:
            grunnstoffer[grunnstoff].append(nukleontall)
        except KeyError:
            grunnstoffer[grunnstoff] = [nukleontall]

    @classmethod
    def hent_isotoper(cls):
        return cls.isotoper

    @classmethod
    def hent_grunnstoffer(cls):
        return cls.grunnstoffer