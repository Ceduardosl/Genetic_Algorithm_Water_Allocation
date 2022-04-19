#%%
import numpy as np

def fitness_measure(inputs, pop):
    #Cálculo da aptidão de cada indíviduo contido na população
    fitness = np.sum(pop*inputs, axis = 1)

    return (fitness)



#%%
if __name__ == "__main__":
    # Y = w1x1 + w2x2 + w3x3 + w4x4 + w5x5 + w6x6
    inputs = [4, -2, 3.5, 5, -11, -4.7]
    
    genes = 6 #Coeficientes para otimizar

    num_chrom = 8 #num chromossomes
 
    pop = np.random.uniform(low = -4, high = 4, size = (8,6))

    print(fitness_measure(inputs, pop))
#%%