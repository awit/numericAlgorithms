#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Mathe 3 - Numerische Algebra - Umweltinformatik WS2013/2014
#Knut Hühne, Simon Brietzke, Laura Moldenhauer, Andreas Werner

#Kompatibilität zwische Ironpython und Python3.3 herstellen
try: input=raw_input
except NameError: pass

Antwort= "\nDas interpolierte Polynom lautet:\np(x)"
Abbild = Antwort

def clear(numlines=320):
        print('\n' * numlines)

def PauseInput():
    eingabe=input("\nFortsetzen mit der Enter-Taste")
    clear()

def output(title,xValues,yValues):
    print(title)
    print("\nStützstellen:")
    ausgabe = []
    for index, xValue in enumerate(xValues):
        ausgabe.append("x%i=%s  " % (index, str(xValue)))
    ausgabe = ''.join(ausgabe)
    print(ausgabe)
    print("Stützwerte:")
    ausgabe = []
    for index, yValue in enumerate(yValues):
        ausgabe.append("y%i=%s  " % (index, str(yValue)))
    ausgabe = ''.join(ausgabe)
    print(ausgabe)
    print


class Polynom:

    def __init__(self, coefficients):
        self.coefficients = coefficients

    def __repr__(self):
        rep = []
        for index, coefficient in enumerate(self.coefficients):
            if index==0 and index!=len(self.coefficients)-1:
                if coefficient<0:
                    rep.append(" - %s" % (str(-1*coefficient)) )
                elif coefficient>0:
                    rep.append(" + %s" % (str(coefficient)) )
            elif index==0 and index==len(self.coefficients)-1:
                if coefficient<0:
                    rep.append(" - %s" % (str(-1*coefficient)) )
                elif coefficient>=0:
                    rep.append(" %s" % (str(coefficient)) )
            elif index==1 and index!=len(self.coefficients)-1:
                if coefficient<0:
                    rep.append(" - %sx" % (str(-1*coefficient)) )
                elif coefficient>0:
                    rep.append(" + %sx" % (str(coefficient)) )
            elif index==1 and index==len(self.coefficients)-1:
                if coefficient<0:
                    rep.append(" - %sx" % (str(-1*coefficient)) )
                elif coefficient>0:
                    rep.append(" %sx" % (str(coefficient)) )
            elif index<len(self.coefficients)-1:
                if coefficient<0:
                    rep.append(" - %sx^%i" % (str(-1*coefficient), index))
                elif coefficient>0:
                    rep.append(" + %sx^%i" % (str(coefficient), index))
            else:
                if coefficient!=0:
                    rep.append(" %sx^%i" % (str(coefficient), index))
        if len(rep)==0:
            rep.append(" 0")
        rep.append(Abbild+'=')
        rep.reverse()
        return ''.join(rep)

    def calculate_at(self, position):
        sol = 0
        for index, coefficient in enumerate(self.coefficients):
            sol += coefficient * position ** index
        return sol

    def __add__(self, summand):
        own_grad = len(self.coefficients)
        summand_grad = len(summand.coefficients)
        if own_grad > summand_grad:
            newCoefficients = list(self.coefficients)
            for i in range(0, summand_grad):
                newCoefficients[i] += summand.coefficients[i]
        else:
            newCoefficients = list(summand.coefficients)
            for i in range(0, own_grad):
                newCoefficients[i] += self.coefficients[i]
        return Polynom(newCoefficients)

    def __mul__(self, number):
        coefficients = [number * float(coefficient) for coefficient in self.coefficients]
        return Polynom(coefficients)


def createPolynomFromNull(nullstellen):
    grad = len(nullstellen) + 1
    coefficients = [0 for pos in range(0, len(nullstellen) + 1)]
    for i in range(0, grad):
        comb = combinate(nullstellen, i)
        for parts in comb:
            coefficients[i] += product(parts)
        #coefficients[i] *= (-1) ** (grad - 1 - i)

    multiplier = 1
    for i in range(0, grad):
        coefficients[i] *= multiplier
        multiplier *= -1

    coefficients.reverse()

    for nullstelle in nullstellen:
        try:
            assert Polynom(coefficients).calculate_at(nullstelle) == 0
        except:
            pass

    return Polynom(coefficients)


def newton(xValues, yValues):
    output("Interpolation nach Newton",xValues, yValues)

    allValues = []
    allValues.append(yValues)

    while len(allValues[-1]) > 1:
        lastValues = allValues[-1]

        newValues = [0 for x in range(len(lastValues) -1)]
        top = len(newValues) - 1
        bottom = len(xValues) -1
        for x in range(len(lastValues) - 1, 0, -1):
            try:
                newValues[x - 1] = float((lastValues[x] - lastValues[x - 1])) / (xValues[bottom] - xValues[top])
            except:
                return ("\nDie Interpolierende kann nicht Bestimmt werden, weil\nStetigkeit und/oder Differenzierbarkeit bzw. paarweise Verschiedenheit\nnicht gewährt ist.")
            top -= 1
            bottom -= 1
        allValues.append(newValues)

    finalPolynomCoefficients = [value[0] for value in allValues]
    print("Newtonkoeffizienten c:")
    ausgabe = []
    for index, Value in enumerate(finalPolynomCoefficients):
        ausgabe.append("c%i=%s  " % (index, str(Value)))
    ausgabe = ''.join(ausgabe)
    print(ausgabe)

    polynoms = []
    polynoms.append(Polynom([1]))   # N0(x) ist immer 1
    for i in range(1, len(xValues)):
        upTo = xValues[:i]
        polynoms.append(createPolynomFromNull(upTo))
    print("Newtonpolynome N(x):")
    ausgabe = []
    global Abbild
    for index, Value in enumerate(polynoms):
        Abbild = ("N%i(x)" %index)
        ausgabe.append("%s" % str(Value))
    Abbild = Antwort
    ausgabe = '\n'.join(ausgabe)
    print(ausgabe)


    for i in range(0, len(polynoms)):
        polynoms[i] = polynoms[i] * finalPolynomCoefficients[i]

    finalPolynom = polynoms[0]
    for i in range(1, len(polynoms)):
        finalPolynom = finalPolynom + polynoms[i]

    for i in range(0, len(xValues)):
        try:
            assert yValues[i] == finalPolynom.calculate_at(xValues[i])
        except:
            pass

    return finalPolynom


def lagrange(xValues, yValues):
    output("Interpolation nach LaGrange",xValues, yValues)
    
    starValues = []
    starPolynoms = []
    for x in xValues:
        copy = list(xValues)
        copy.remove(x)
        pol = createPolynomFromNull(copy)
        starPolynoms.append(pol)
        sol = pol.calculate_at(x)
        starValues.append(sol)

    polynoms = []
    
    global Abbild
    print("Lagrange Polynome:")
    for i in range(0, len(starValues)):
        try:
            polynoms.append(starPolynoms[i] * (1.0 / float(starValues[i]) * yValues[i]))
            Abbild = ("l%i(x)" % i)
            print("d%i=%s , y%i=%s ; %s" % (i,str(starValues[i]),i,str(yValues[i]),str(starPolynoms[i])))
        except:
            return ("\nDie Interpolierende kann nicht Bestimmt werden, weil\nStetigkeit und/oder Differenzierbarkeit bzw. paarweise Verschiedenheit\nnicht gewährt ist.")

    Abbild = Antwort
    finalPolynom = polynoms[0]
    for i in range(1, len(polynoms)):
        finalPolynom += polynoms[i]

    return finalPolynom


def product(collection):
    product = 1
    for element in collection:
        product *= element
    return product


def combinate(inputValues, k):

    del returnValues[:]
    del combinationValues[:]

    if k == 0:
        yield [1]
        raise StopIteration

    for i in range(0, k):
        combinationValues.append(0)

    combinations2(inputValues, len(inputValues), k)

    for x in returnValues:
        yield x


def combinations2(inputValues, n, k):
    global combinationValues

    if n > k:
        combinations2(inputValues, n - 1, k)

    if n >= k:
        k -= 1
        combinationValues[k] = inputValues[n - 1]
        if k > 0:
            combinations2(inputValues, n - 1, k)
        else:
            rep = [value for value in combinationValues]

            returnValues.append(rep)


combinationValues = []
returnValues = []

if __name__ == "__main__":
    clear()
    pol = newton([1, 3, 5], [3, 7, 14])
    print(pol)
    PauseInput()

    pol2 = newton([-2, -1, 0, 1, 2, 3], [-16, 0, 4, 8, 0, 64])
    print(pol2)
    PauseInput()

    lagrange_pol = lagrange([1, 3, 5], [3, 7, 14])
    print(lagrange_pol)
    PauseInput()
    
    weitere = True
    while weitere:
        eigeneEingabe=input("\nMöchten Sie selber Stützstellen und Stützwerte eingeben (J/N)? ")
        if eigeneEingabe.upper()=="J":
            print("Bitte geben Sie die Stützstellen ein (leere Eingabe beendent die Erfassung:")
            eingabe=True
            index=0
            inputx=[]
            while eingabe:
                retry=True
                while retry:
                    s_eingabe=input("x%i=" % index)
                    if s_eingabe!="":
                        try:
                            f_eingabe = float(s_eingabe)
                            retry=False
                            inputx.append(f_eingabe)
                            index+=1
                        except ValueError:
                            print("Bitte überprüfen Sie ihre Eingabe!")
                    if s_eingabe=="":
                        retry=False
                        eingabe=False
            clear()
            print("Bitte geben Sie die Stützwerte ein (leere Eingabe beendent die Erfassung:")
            eingabe=True
            index=0
            inputy=[]
            while eingabe:
                retry=True
                while retry:
                    s_eingabe=input("y%i=" % index)
                    if s_eingabe!="":
                        try:
                            f_eingabe = float(s_eingabe)
                            retry=False
                            inputy.append(f_eingabe)
                            index+=1
                        except ValueError:
                            print("Bitte überprüfen Sie ihre Eingabe!")
                    if s_eingabe=="":
                        retry=False
                        eingabe=False
            clear()
            if len(inputx)==len(inputy) and len(inputx)>0:
                pol3=newton(inputx,inputy)
                print(pol3)
                PauseInput()
                pol3=lagrange(inputx,inputy)
                print(pol3)
            else:
                print("Die Anzahl der angegebenen Stützstellen und Stützwerte \nmuss mindestens 1 und gleich sein!!!")
        elif eigeneEingabe.upper()=="N":
            print("Auf Wiedersehen !!!")
            weitere=False
