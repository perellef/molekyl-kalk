# Molekyl-verktøy

En samlemappe som inneholder en rekke nyttige verktøy relatert molekyl- og isotopformler.

## Isotopformel

Instansiering av isotopformel:

```python
isotop = Molekylformel.lag("12C1 1H3 2H1")
isotop = Molekylformel.lag(["12C1","1H3","2H1"])
isotop = Molekylformel.lag({"12C": 1,"1H": 3,"2H": 1})
```

Nyttige operasjoner:

```python
isotop.beregn_forekomst()
isotop.beregn_masse()
isotop.kalkuler_DBE()
```

## Molekylformel

Instansiering av molekylformel:

```python
ch4 = Molekylformel.lag("C1 H4")
ch4 = Molekylformel.lag(["C1","H4"])
ch4 = Molekylformel.lag({"C": 1,"H": 4})
ch4 = Molekylformel.fra_isotopformel(isotop)
```

Nyttige operasjoner:

```python
ch4.beregn_masse()
ch4.kalkuler_DBE()
ch4.print_isotopfordeling()
ch4.hent_isotopfordeling()
```

## Formelkalkulator

Inneholder verktøy for deduksjon av isotopformel fra masse, enten med eller uten MS-topper.

#### Uten MS-topper

Først må et masseintervall defineres.

```python
fra, til = [120.05, 120.07]  # definerer et masseintervall
```

Dette kan så direkte brukes til å finne isotopformlene som oppfyller betingelsene.

```python
Formelkalkulator.fra_masseintervall(fra,til)
```

Output:
```bashh
Det finnes 28 slike isotopformler mellom [120.05,120.07]. De 10 mest sannsynlige er:
----------------------
120.057515 g/mol | 5.0 DBE | ¹²C₈ ¹H₈ ¹⁶O₁
120.068748 g/mol | 5.0 DBE | ¹²C₇ ¹H₈ ¹⁴N₂
120.066069 g/mol | 0.5 DBE | ¹²C₄ ¹H₁₀ ¹⁴N₁ ¹⁶O₃
120.060889 g/mol | 0.0 DBE | ¹²C₅ ¹H₁₂ ¹⁶O₁ ³²S₁
120.056172 g/mol | 5.5 DBE | ¹²C₆ ¹H₆ ¹⁴N₃
120.053493 g/mol | 1.0 DBE | ¹²C₃ ¹H₈ ¹⁴N₂ ¹⁶O₃
120.058660 g/mol | 1.0 DBE | ¹²C₅ ¹H₉ ¹⁶O₂ ¹⁹F₁
120.058004 g/mol | 0.5 DBE | ¹²C₅ ¹H₁₁ ¹⁴N₁ ³⁵Cl₁
120.060657 g/mol | 0.0 DBE | ¹²C₄ ¹H₁₂ ¹⁶O₂ ²⁸Si₁
120.051971 g/mol | 1.0 DBE | ¹²C₆ ¹H₁₁ ³⁷Cl₁
````

#### Med MS-topper

Definer først M+ signal og de største toppene fra massespekteret. Signalene defineres enten gjennom masseintervall, eller heltall (der intervalet settes til ± 0.5).
```python
M = 76
topper = [61,43,33,15]
```
Fra dette kan isotopformler regnes ut.
```python
Formelkalkulator.fra_MS_topper(M,topper)
```
NB: Merk at m/z-signaler med z ≥ 2 ikke tas hensyn til, og kan medføre feil.
## 
Ved store intervaller eksisterer det ofte mange formler. Det vil da være nyttig å begrense noen av isotopene.


```python
begrensning = {}

begrensning["1H"]  = "8-",   #  ->  høyst 8
begrensning["12C"] = "3+",   #  ->  minst 3
begrensning["14N"] = "1",    #  ->  nøyaktig 1
begrensning["16O"] = "1-3",  #  ->  mellom 1 og 3, inklusivt

Formelkalkulator.fra_masseintervall(masse, begr=begrensning)
```

