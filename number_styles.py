import codecs

def digitToWord(digit):
    digits = {0: 'ZERO', 1: 'ONE', 2: 'TWO', 3: 'THREE', 4: 'FOUR', 5: 'FIVE', 6: 'SIX', 7: 'SEVEN', 8:'EIGHT', 9:'NINE'}
    
    return digits[digit]

def superscript(number):
    output = ""
    for digit in str(number):
        output += r"\N{SUPERSCRIPT " + digitToWord(int(digit)) + "}"
    return codecs.decode(output, 'unicode_escape')

def subscript(number):
    output = ""
    for digit in str(number):
        output += r"\N{SUBSCRIPT " + digitToWord(int(digit)) + "}"
    return codecs.decode(output, 'unicode_escape')