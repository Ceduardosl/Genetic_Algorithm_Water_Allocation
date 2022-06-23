"""
Optimal Water Allocation based on Genetic Algorithm
"""

__author__ = ["Carlos Eduardo Sousa Lima"]
__license__ = "GPL"
__version__ = "1.0"
__email__ = "eduardolima@alu.ufc.br"


#%%
import numpy as np
import pandas as pd
import os

def fitness_measure(weights, inputs, pop, v_aloc):
    #Cálculo da aptidão de cada indíviduo contido na população

    # fitness = 1/np.sum(weights * np.power(inputs-pop,2), axis = 1)
    """
    Removi o quadrado da diferença, pois a função não percebia quando
    alocava mais água do que o solicitado, fazendo com que atribuisse
    valores positivos e altos de eficiência, perpetuando para os herdeiros
    essa característica
    """
    fitness = 1/np.sum(weights * np.power(inputs-pop,2), axis = 1)
    over_id = np.where(pop[:,0] + pop[:,1] > v_aloc)

    fitness[over_id] = np.NINF
    # print(over_id)
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
        mutated_gene_id = np.random.randint(0, offsprings.shape[1])
        if np.random.random() < mutation_rate:
            """
            No primeiro caso, um gene específico sofre uma mutação aleatória
            No segundo caso, um cromossomo específico sofre alteração em toda sua composição
            mutation_param[0] = número de genes
            mutation_param[1] = valor limite da soma
            """
            offsprings[i, mutated_gene_id] = offsprings[i, mutated_gene_id] + np.random.uniform(-0.1, 0.1)*offsprings[i, mutated_gene_id]
            
            #offsprings[i, mutated_gene_id] = offsprings[i, mutated_gene_id] + np.random.uniform(-offsprings[i, mutated_gene_id]/2, offsprings[i, mutated_gene_id]/2)
            # print("Mutation occurred in the {}th offspring in the {}th gene".format(i+1, mutated_gene_id+1))
            # offsprings[i, ] = np.random.dirichlet(np.ones(mutation_param[0])) * mutation_param[1]
            # print("Mutation occurred in the {}th offspring".format(i+1))
            count += 1                
    #if count == 0:
        #print("No Mutation ocurred")
    return (offsprings)

def generate_pop(v_aloc, shape):

    for i in range(shape[1]):
        if i == 0:
            x = np.random.uniform(0.25, 0.75, (shape[0],1))
        else:
            x = np.concatenate((x, np.random.uniform(0.25, 0.75, (shape[0],1))), axis = 1)
    pop = x * v_aloc

    return (pop) 

# def generate_pop2(target, shape):
#     """
#     Função para gerar dois valores aleatórios dentrod e um range de possíveis valores
#     de forma que a soma dos valores seja igual a um valor alvo
#     Importante: Organizar as demanads de acordo com a mais solicitada
#     """
#     x = np.random.dirichlet(np.ones(shape[1]), shape[0])
#     pop = np.fliplr(np.sort(x * target))
#     return pop

#%%
if __name__ == "__main__":

    data = pd.read_csv("{}/Outputs/Annual_reg_65_garan_0.80.csv".format(os.getcwd()))
    reg_ts = data[["Data", "vol_reg"]]
    reg = 780
    
    demandas = {
        "urbana":[1, 0.1340*reg],
        "rural":[1, 0.0685*reg],
        "animal":[1, 0.1240*reg],
        "agricultura":[3, 0.6690*reg],
        "industria":[1, 0.0045*reg]
    }
    # Organizar os inputs em ordem crescente de demanda
    inputs = [demandas["agricultura"][1], demandas["industria"][1]]
    weights = [demandas["agricultura"][0], demandas["industria"][0]]
    rate = [0.6690, 0.0045]

    num_genes, num_chrom, num_parents,  mutation_rate = [2, 8, 4, 0.60]
    n_off = (num_chrom - num_parents, num_genes)
    # População Inicial - Matriz (8x2)
    num_generation = 10000

#%%
for i in range(0, len(reg_ts)):
    aloc_dict = {
        "urbana": 0,
        "rural": 0,
        "animal": 0,
        "agricultura": 0,
        "industria":0 
        }
    if reg_ts.iloc[i,1] == 780:
        print("{} - Todos os usos atentidos".format(reg_ts.iloc[i,0]))
        aloc_dict["urbana"] = 0.1340*reg
        aloc_dict["rural"] = 0.0685*reg
        aloc_dict["animal"] = 0.1240*reg
        aloc_dict["agricultura"] = 0.6690*reg
        aloc_dict["industria"] = 0.0045*reg
        print(aloc_dict)
        print("Vazão Alocada = {:.2f} hm³\n".format(sum(aloc_dict.values())))

    elif reg_ts.iloc[i,1] < (demandas["urbana"][1] + demandas["rural"][1] + demandas["animal"][1]):
        print("{} - Usos prioritários parcialmente ou não atendidos".format(reg_ts.iloc[i,0]))
        print("{:.2f} hm³ para alocar entre os usos prioritários, que exigem {:.2f} hm³\n".format(reg_ts.iloc[i,1], demandas["urbana"][1] + demandas["rural"][1] + demandas["animal"][1]))
    else:
        print("{} - Otimizar usos não prioritários".format(reg_ts.iloc[i,0]))
        aloc_dict["urbana"] = 0.1340*reg
        aloc_dict["rural"] = 0.0685*reg
        aloc_dict["animal"] = 0.1240*reg
        v_aloc = reg_ts.iloc[i,1] - (aloc_dict["urbana"] + aloc_dict["rural"] + aloc_dict["animal"])
        print("Volume de {:.2f} hm³ para alocar entre usos não prioritários".format(v_aloc))
        pop = generate_pop(v_aloc, [num_chrom, num_genes])
        

        for generation in range(num_generation):
            #print("#####-----{}th Geration-----#####\n".format(generation + 1))
            fitness = fitness_measure(weights, inputs, pop, v_aloc)
            best_individuals = select_individuals(pop, fitness, num_parents)
            offsprings = crossover(
                best_individuals,
                num_offspring = (num_chrom - num_parents, num_genes)
                ) 
            mutated_offsprings = mutation(offsprings, mutation_rate)
            pop[0:num_parents, :] = best_individuals
            pop[num_parents:, :] = mutated_offsprings
            #pop = np.where(pop > v_aloc, 0.90 * v_aloc, pop)
            for i in range(pop.shape[1]):
                pop[:,i] = np.where(pop[:,i] > inputs[i], 0.9*inputs[i], pop[:,i])
        print(pop)
        print("fitness")
        print(fitness_measure(weights, inputs, pop, v_aloc))
        print("\n")

#%%