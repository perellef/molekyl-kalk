from register import Register
from isotopformel import Isotopformel
from molekylformel import Molekylformel

import inspect

class Formelkalkulator:

    @classmethod
    def __kalkuler_isotopformler(cls, isotoper, fra, til, formel={},vekt=0):
        
        isotop = next(iter(isotoper))

        masse = isotoper[isotop]["masse"]
        begr = isotoper[isotop]["begrensning"]

        del isotoper[isotop]
        
        if len(isotoper)==0:  # fra,til = 5,6 | vekt = 5 | H = 1 
            minst = min(int((fra-vekt+1e-8)//masse),0)
            maks = int((til-vekt+1e-8)//masse)

            for m in range(minst,maks+1):
                if (vekt+m*masse<fra):
                    continue
                if not Isotopformel.har_antall_av_isotop({isotop: m},isotop,begr):
                    continue
                yield {isotop: m}|formel

        else:
            maks = int((til-vekt+1e-8)//masse)

            for m in range(maks+1):

                if not Isotopformel.har_antall_av_isotop({isotop: m},isotop,begr):
                    continue

                yield from cls.__kalkuler_isotopformler(isotoper.copy(), fra, til, {isotop: m}|formel, vekt+m*masse)

    @classmethod
    def fra_masseintervall(cls,fra,til,*,begr={},n=10):

        standard_begrensning = {"1H": "30-","12C": "20-","14N": "5-","16O": "10-","19F": "5-","28Si": "5-","31P": "5-","32S": "5-","34S": "1-","35Cl": "5-","37Cl": "2-","79Br": "5-","81Br": "5-","127I": "5-"}

        utskrift = "Begrensningene er satt til:"

        isotoper = {}
        for isotop in reversed(standard_begrensning):
            if isotop in begr: 
                begrensning = begr[isotop]
            else:
                begrensning = standard_begrensning[isotop]

            utskrift += f"\n{isotop}: {begrensning}"
            
            isotoper[isotop] = Register.isotoper[isotop] | {"begrensning": begrensning}

        utskrift += "\nandre isotoper: 0\n.... \n"

        if inspect.stack()[1][3]!="fra_MS_topper": # om funksjon kalt av brukeren
            print(utskrift)
        
        formler = [el for el in cls.__kalkuler_isotopformler(isotoper,fra,til)] # genererer gyldige isotopformler

        isotopformler = [Isotopformel(formel) for formel in formler]

        isotopformler = [formel for formel in isotopformler if formel.kalkuler_DBE()>=0]
        
        if inspect.stack()[1][3]!="fra_MS_topper": # om funksjon kalt av brukeren
            print(f"Det finnes {len(isotopformler)} slike isotopformler mellom [{fra},{til}]. De {n} mest sannsynlige er:")
            print("----------------------")
            for isotopformel in sorted(isotopformler, key = lambda x: -x.score())[:n]:

                masse = "{:.6f}".format(isotopformel.beregn_masse())
                dbe = isotopformel.kalkuler_DBE()

                print(f"{masse} g/mol | {dbe} DBE | {isotopformel}")

                
        return isotopformler

    @classmethod
    def fra_MS_topper(cls,M_topp,topper,*,begr={},n=10):

        if isinstance(M_topp,int):
            fra = M_topp - 0.5
            til = M_topp + 0.5
        elif isinstance(M_topp,list):
            fra,til = M_topp
        else:
            raise TypeError

        formler = cls.fra_masseintervall(fra,til,begr=begr)

        gyldige_formler = []
        for formel in formler:

            begr = {}
            for isotop,antall in formel.hent_formel().items():
                begr[isotop] = str(antall)+"-"

            ut = False

            for topp in topper:
                if isinstance(topp,int):
                    topp_fra = topp - 0.5
                    topp_til = topp + 0.5
                elif isinstance(topp,list):
                    topp_fra,topp_til = topp

                formler_topp = cls.fra_masseintervall(topp_fra,topp_til,begr=begr)

                if len(formler_topp)==0:
                    ut = True
                    break
            
            if not ut:
                gyldige_formler.append(formel)


        isotopformler = [formel for formel in gyldige_formler]
                
        print(f"Det finnes {len(isotopformler)} isotopformler som passer med MS-toppene. De {n} mest sannsynlige er:")
        print("----------------------")
        for isotopformel in sorted(isotopformler, key = lambda x: -x.score())[:n]:
            masse = "{:.6f}".format(isotopformel.beregn_masse())
            dbe = isotopformel.kalkuler_DBE()

            print(f"{masse} g/mol | {dbe} DBE | {isotopformel}")