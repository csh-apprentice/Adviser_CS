'''
This function produces a random three section melt function. On top of that it uses a random alpha to smooth the melt. The results is dictionary 
with the meltparams. The functions include a total melt so the model renormalizes to the total melt.
'''
import random
def polyEval(poly, x):
    return poly['coeff'][1] * x + poly['coeff'][0]
#
def randomMelt(totalMelt):
    poly1, poly2, poly3 = {'deg': 1}, {'deg': 1}, {'deg': 1}
    alpha = random.uniform(200, 4000)
    filterWithFloatMask = random.uniform(0,1) > 0.5
    meltR = {'alpha': alpha, 'filterWithFloatMask': filterWithFloatMask, 'totalMelt': totalMelt, 
             'numberOfPolynomials': 3, 'poly1': poly1, 'poly2': poly2, 'poly3': poly3}
    # set random contiguous range for polynomials
    poly1['min'], poly1['max'] = 0, random.uniform(0, 0.5)*650 + 300
    poly2['min'], poly2['max'] = poly1['max'], random.uniform(-0.5, 0.5) * 400 + 850
    poly3['min'], poly3['max'] = poly2['max'], 20000
    poly1['coeff'] = [random.uniform(-8, 0), random.uniform(-0.3, -0.1)]
    #
    for polyNew, polyLast in zip([poly2, poly3], [poly1, poly2]):
        slope = random.uniform(-4, -0.5) 
        intercept =  polyEval(polyLast, polyLast['max']) - slope * polyLast['max']
        polyNew['coeff'] = [intercept, slope]
    return meltR

