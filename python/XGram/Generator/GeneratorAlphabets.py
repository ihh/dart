## Generator for alphabets
import XGram.Model

class Alphabet:
    def __init__(self):
        pass

class AlphabetDNA:
    
    def __init__(self):
        pass
    
    def buildGrammar(self, alphabet):
        alphabet.setName( "DNA" )
        alphabet.setTokens( ("a", "c", "g", "t") )
        alphabet.setComplements( ("t", "g", "c", "a") )        
        alphabet.setWildcard( "*" )        
        alphabet.addExtension("n", ("a", "c", "g", "t" ))
        alphabet.addExtension("x", ("a", "c", "g", "t" ))        
        alphabet.addExtension("u", ("t",))
        alphabet.addExtension("r", ("a","g"))
        alphabet.addExtension("y", ("c","t"))
        alphabet.addExtension("m", ("a","c"))
        alphabet.addExtension("k", ("g","t"))
        alphabet.addExtension("s", ("c","g"))
        alphabet.addExtension("w", ("a","t"))
        alphabet.addExtension("h", ("a","c","t"))
        alphabet.addExtension("b", ("c","g","t"))
        alphabet.addExtension("v", ("a","c","g"))
        alphabet.addExtension("d", ("a","g","t"))
    
class AlphabetRNA:
    def __init__(self):
        pass
    
    def buildGrammar(self, alphabet):
        alphabet.setName( "RNA" )
        alphabet.setTokens( ("a", "c", "g", "u") )
        alphabet.setComplements( ("u", "g", "c", "a") )        
        alphabet.setWildcard( "*" )        
        alphabet.addExtension("n", ("a", "c", "g", "u" ))
        alphabet.addExtension("x", ("a", "c", "g", "u" ))        
        alphabet.addExtension("t", ("u",))
        alphabet.addExtension("r", ("a","g"))
        alphabet.addExtension("y", ("c","u"))
        alphabet.addExtension("m", ("a","c"))
        alphabet.addExtension("k", ("g","u"))
        alphabet.addExtension("s", ("c","g"))
        alphabet.addExtension("w", ("a","u"))
        alphabet.addExtension("h", ("a","c","u"))
        alphabet.addExtension("b", ("c","g","u"))
        alphabet.addExtension("v", ("a","c","g"))
        alphabet.addExtension("d", ("a","g","u"))

class AlphabetProtein(Alphabet):
    
    def __init__(self):
        pass

    def buildGrammar(self, alphabet):
        alphabet.setName( "Protein" )
        alphabet.setTokens( ("a","r","n","d","c","q","e","g","h","i","l","k","m","f","p","s","t","w","y","v"))
        alphabet.setWildcard( "*" )
        alphabet.addExtension( "x", ("a","r","n","d","c","q","e","g","h","i","l","k","m","f","p","s","t","w","y","v"))
        alphabet.addExtension( "b", ("n", "d") )
        alphabet.addExtension( "z", ("q", "e") )\

        
if __name__ == "__main__":
    
    alphabet = XGram.Model.Alphabet.Alphabet()
    generator = AlphabetProtein()
    generator.buildGrammar( alphabet )

    print alphabet.getGrammar()
    
    alphabet = XGram.Model.Alphabet.Alphabet()
    generator = AlphabetDNA()
    generator.buildGrammar( alphabet )
    
    print alphabet.getGrammar()

    alphabet = XGram.Model.Alphabet.Alphabet()
    generator = AlphabetRNA()
    generator.buildGrammar( alphabet )
    
    print alphabet.getGrammar()
    