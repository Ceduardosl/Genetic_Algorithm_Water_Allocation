# Adapted of:
# Ahmed Gad Tutorial
# https://towardsdatascience.com/genetic-algorithm-implementation-in-python-5ab67bb124a6

__author__ = ["Carlos Eduardo Sousa Lima"]
__license__ = "GPL"
__version__ = "1.0"
__email__ = "eduardolima@alu.ufc.br"
#%%
import numpy as np
import pandas as pd
from sympy import cosine_transform

def fitness_measure(inputs, pop):
    #Cálculo da aptidão de cada indíviduo contido na população
    fitness = np.sum(pop*inputs, axis = 1)

    return (fitness)

def select_individuals(pop, fitness, num_parents):
    
    #Criar matriz vazia para melhores indivíduos
    best_individuals = np.empty((num_parents, pop.shape[1]))

    for i in range(num_parents):
        max_id = np.where(fitness == np.max(fitness))[0][0]
        best_individuals[i, :] = pop[max_id, :]
        # Sem a linha de comando abaixo, será selecionado sempre o melhor indivíduo
        # Compromentedno o processo de crossover
        fitness[max_id] = np.NINF
    
    return (best_individuals)

def crossover(best_individuals, num_offspring):

    offsprings = np.empty(num_offspring)

    # num_offspring[0] = num de descendentes
    # num_offspring[1] = num de genes dos descendentes

    for i in range(num_offspring[0]):
        
        #randomico a partir de 1, sempre haverá cruzamento
        split_point = np.random.randint(1, num_offspring[1])

        id_parent1 = i % best_individuals.shape[0]
        id_parent2 = (i + 1) % best_individuals.shape[0] 

        offsprings[i, 0:split_point] = best_individuals[id_parent1, 0:split_point]
        offsprings[i, split_point:] = best_individuals[id_parent2, split_point:]
    
    return (offsprings)

def mutation(offsprings, mutation_rate):
    count = 0
    for i in range(offsprings.shape[0]):
        mutated_gene_id = np.random.randint(1, offsprings.shape[1])
        if np.random.random() < mutation_rate:
            offsprings[i, mutated_gene_id] = offsprings[i, mutated_gene_id] + np.random.uniform(low = -1, high = 1, size = None)
            print("Mutation occurred in the {}th offspring in the {}th gene".format(i+1, mutated_gene_id))
            count += 1                
    if count == 0:
        print("No Mutation ocurred")
    return (offsprings)

#%%
if __name__ == "__main__":
    # Y = w1x1 + w2x2 + w3x3 + w4x4 + w5x5 + w6x6
    inputs = [4, -2, 3.5, 5, -11, -4.7]
    
    # Cada Cromossomos possui 6 genes
    num_genes = 6 #Coeficientes para otimizar

    # Número de Cromossomos, cada cromossomo associado a um indivíduo
    num_chrom = 8

    # Número de indivíduos selecionados para o crossover
    num_parents = 4

    # Taxa de mutação
    mutation_rate = 0.4

    # População Inicial - Matriz (8x6)
    pop = np.random.uniform(low = -4, high = 4, size = (num_chrom, num_genes))

    num_generation = 15

    for generation in range(num_generation):

        print("#####-----{}th Geration-----#####\n".format(generation + 1))
        print("Population")
        print(pop)
        print("\n")

        fitness = fitness_measure(inputs, pop)
        print("Fitness fo individuals: ")
        print("Best Fitness = {}".format(np.max(fitness)))
        print(fitness)
        print("\n")


        best_individuals = select_individuals(pop, fitness, num_parents)
        print("Best individuals:")
        print(best_individuals)
        print("\n")

        n_off = (num_chrom - num_parents, num_genes)
        # o num de descendentes é igual ao num de indivíduos não selecionados
        offsprings = crossover(
            best_individuals,
            num_offspring = (num_chrom - num_parents, num_genes)
            )
        print("Offsprings")
        print(offsprings)
        print("\n")

        mutated_offsprings = mutation(offsprings, mutation_rate)
        print("\nMutated Offsprings")
        print(mutated_offsprings)
        print("\n")

        pop[0:num_parents, :] = best_individuals
        pop[num_parents:, :] = mutated_offsprings
#%%

    
