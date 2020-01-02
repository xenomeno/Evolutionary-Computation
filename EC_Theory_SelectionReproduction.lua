dofile("Statistics.lua")
dofile("Bitmap.lua")
dofile("Graphics.lua")
dofile("EC_Common.lua")
dofile("GA_Common.lua")

local RUNS                          = 400
local POPULATION_SIZE               = 50
local MAX_GENERATIONS               = 80
local MAX_BIRTHS                    = 10000

local GENOME_REALS                  = 2
local GENOME_REAL_MIN               = -3.0
local GENOME_REAL_MAX               = 4.0
local GENOME_BITS                   = 17
local GENOME_MAX_VALUE              = math.pow(2, GENOME_BITS) - 1
local GENOTYPE_LENGTH               = 2 * GENOME_BITS

local IMAGE_WIDTH                   = 1000
local IMAGE_HEIGHT                  = 1000
local IMAGE_FILENAME_FIGURE_6_30    = "EC/Figure6_30.bmp"
local IMAGE_FILENAME_FIGURE_6_31    = "EC/Figure6_31.bmp"
local IMAGE_FILENAME_FIGURE_6_32    = "EC/Figure6_32.bmp"
local IMAGE_FILENAME_FIGURE_6_33    = "EC/Figure6_33.bmp"
local IMAGE_FILENAME_FIGURE_6_34    = "EC/Figure6_34.bmp"
local IMAGE_FILENAME_FIGURE_6_35    = "EC/Figure6_35.bmp"
local IMAGE_FILENAME_FIGURE_6_36    = "EC/Figure6_36.bmp"
local IMAGE_FILENAME_FIGURE_6_37    = "EC/Figure6_37.bmp"
local IMAGE_FILENAME_FIGURE_6_38    = "EC/Figure6_38.bmp"
local IMAGE_FILENAME_FIGURE_6_39    = "EC/Figure6_39.bmp"
local IMAGE_FILENAME_FIGURE_6_40    = "EC/Figure6_40.bmp"
local IMAGE_FILENAME_FIGURE_6_41    = "EC/Figure6_41.bmp"
local IMAGE_FILENAME_FIGURE_6_42    = "EC/Figure6_42.bmp"
local IMAGE_FILENAME_FIGURE_6_43    = "EC/Figure6_43.bmp"
local IMAGE_FILENAME_FIGURE_6_44    = "EC/Figure6_44.bmp"
local IMAGE_FILENAME_FIGURE_6_45    = "EC/Figure6_45.bmp"
local IMAGE_FILENAME_FIGURE_6_46    = "EC/Figure6_46.bmp"
local IMAGE_FILENAME_FIGURE_6_47    = "EC/Figure6_47.bmp"
local IMAGE_FILENAME_FIGURE_6_48    = "EC/Figure6_48.bmp"
local IMAGE_FILENAME_FIGURE_6_49    = "EC/Figure6_49.bmp"
local IMAGE_FILENAME_FIGURE_6_49    = "EC/Figure6_49.bmp"
local IMAGE_FILENAME_FIGURE_6_50    = "EC/Figure6_50.bmp"
local IMAGE_FILENAME_FIGURE_6_51    = "EC/Figure6_51.bmp"
local IMAGE_FILENAME_FIGURE_6_52    = "EC/Figure6_52.bmp"
local IMAGE_FILENAME_FIGURE_6_53    = "EC/Figure6_53.bmp"

local function RealToBinary(x)
  return math.floor(GENOME_MAX_VALUE * (x - GENOME_REAL_MIN) / (GENOME_REAL_MAX - GENOME_REAL_MIN))
end

local function BinaryToReal(x)
  return GENOME_REAL_MIN + (GENOME_REAL_MAX - GENOME_REAL_MIN) * x / GENOME_MAX_VALUE
end

local function Fitness(genotype)
  local result = 0
  for bit_pos = 1, GENOTYPE_LENGTH, GENOME_BITS do
    local bits = ExtractBitstring(genotype, bit_pos, GENOME_BITS)
    local x = BinaryToReal(bits[1])
    result = result + x * x
  end
  
  return result
end

local function GenRandomPop(size)
  local pop = {}
  while #pop < size do
    local bitstring = GenRandomBitstring(GENOTYPE_LENGTH)
    local genotype = PackBitstring(bitstring)
    table.insert(pop, {genotype = genotype, fitness = Fitness(genotype)})
  end
  
  return pop
end

local function IsPopConverged(pop)
  local genotype = pop[1].genotype
  for i = 2, #pop do
    if GetCommonBits(genotype, pop[i].genotype) < GENOTYPE_LENGTH then
      return false
    end
  end
  
  return true
end

local function EvalPop(pop)
  local min_fitness, max_fitness
  local repro_successes = 0
  local total_fitness, total_fitness_change = 0, 0
  for _, ind in ipairs(pop) do
    local fitness = ind.fitness
    min_fitness = (not min_fitness or fitness < min_fitness) and fitness or min_fitness
    max_fitness = (not max_fitness or fitness > max_fitness) and fitness or max_fitness
    total_fitness = total_fitness + fitness
    local fitness_change = ind.parent_fitness and (ind.fitness - ind.parent_fitness) or 0
    total_fitness_change = total_fitness_change + fitness_change
    repro_successes = repro_successes + ((fitness_change > 0) and 1 or 0)
  end
  pop.min_fitness, pop.max_fitness, pop.avg_fitness = min_fitness, max_fitness, total_fitness / #pop
  pop.avg_fitness_change, pop.repro_successes = total_fitness_change / #pop, repro_successes
end

local function RMS_Hamming(pop)
  local rms = 0
  for i, ind1 in ipairs(pop) do
    local genotype1 = ind1.genotype
    for j = i + 1, #pop do
      rms = rms + GENOTYPE_LENGTH - GetCommonBits(genotype1, pop[j].genotype)
    end
  end
  
  return (2 * rms) / (#pop * (#pop - 1))
end

local function UniformStochasticSelection(parents, size)
  local offsprings = {}
  for i = 1, size do
    local parent = parents[math.random(1, #parents)]
    table.insert(offsprings, parent)
  end
  
  return offsprings
end

local function GenPermutation(size)
  local permutation = {}
  for i = 1, size do
    permutation[i] = i
  end
  for i = 1, #permutation do
    local j = math.random(1, #permutation)
    permutation[i], permutation[j] = permutation[j], permutation[i]
  end
  
  return permutation
end

local function UniformDeterministicSelection(parents, size)
  local perm_idx, permutation = 1, {}
  local offsprings = {}
  while #offsprings < size do
    if perm_idx > #permutation then
      permutation = GenPermutation(#parents)
      perm_idx = 1
    end
    local parent_idx = permutation[perm_idx]
    perm_idx = perm_idx + 1
    table.insert(offsprings, parents[parent_idx])
  end
  
  return offsprings
end

local function BinaryTournamentSelection(parents, size)
  local m = #parents
  local offsprings = {}
  while #offsprings < size do
    local ind1 = parents[math.random(1, m)]
    local ind2 = parents[math.random(1, m)]
    table.insert(offsprings, (ind1.fitness > ind2.fitness) and ind1 or ind2)
  end
  
  return offsprings
end

local function TruncateSelection(parents, size)
  table.sort(parents, function(ind1, ind2) return ind1.fitness > ind2.fitness end)
  local offsprings = {}
  while #offsprings < size do
    local parent = parents[math.random(1, size)]
    table.insert(offsprings, parent)
  end
  
  return offsprings
end

local function CrossoverUniformCreator(rate)
  return function(parent1, parent2)
    return BinaryCrossoverUniform(parent1, parent2, rate)
  end
end

local function MutationCreator(rate)
  return function(genotype)
    return BinaryMutation(genotype, rate)
  end
end

local function Crossover_Blend(parent1, parent2, r)
  local offspring = {}
  for i, allele in ipairs(parent1) do
    offspring[i] = r * allele + (1.0 - r) * parent2[i]
  end
  
  return offspring
end

local function CrossoverBlendCreator(r)
  return function(parent1, parent2) return Crossover_Blend(parent1, parent2, r) end
end

local function EV(m, n, init_pop, on_new_gen, crossover, mutation, selection_parent, selection_survival, overlapping)
  local parents = init_pop(m)
  EvalPop(parents)
  if on_new_gen then
    local result = on_new_gen(0, parents)
    if result == "converged" then
      return 0, parents
    end
  end
  local gen = 1
  while gen <= MAX_GENERATIONS do
    local pool = selection_parent(parents, n)
    local offsprings = {}
    if crossover then
      for i = 1, #pool, 2 do
        local parent1, parent2 = pool[i], (i < #pool) and pool[i + 1] or pool[i]
        local offspring1, offspring2 = crossover(parent1.genotype, parent2.genotype)
        table.insert(offsprings, {genotype = offspring1, fitness = Fitness(offspring1), parent_fitness = parent1.fitness})
        if offspring2 then
          table.insert(offsprings, {genotype = offspring2, fitness = Fitness(offspring2), parent_fitness = parent1.fitness})
        end
      end
    else
      offsprings = pool
    end
    if mutation then
      for i, ind in ipairs(offsprings) do
        local mutant = table.copy(ind, "deep")
        mutant.parent_fitness = ind.fitness
        mutation(mutant.genotype)
        mutant.fitness = Fitness(mutant.genotype)
        offsprings[i] = mutant
      end
    end
    if selection_survival then
      if overlapping then
        for _, ind in ipairs(parents) do
          table.insert(offsprings, ind)
        end
      end
      parents = selection_survival(offsprings, m)
    else
      parents = offsprings
    end
    EvalPop(parents)
    if on_new_gen then
      local result = on_new_gen(gen, parents)
      if result == "converged" then
        break
      end
    end
    gen = gen + 1
  end
  
  return gen, parents
end

local function EV_Real(m, n, on_new_gen, blend, selection_parent, selection_survival, overlapping)
  local function fitness(ind)
    return ind[1] * ind[1] + ind[2] * ind[2]
  end

  local parents = {}
  while #parents < m do
    local ind = {}
    for i = 1, 2 do
      ind[i] = GENOME_REAL_MIN + math.random() * (GENOME_REAL_MAX - GENOME_REAL_MIN)
    end
    ind.fitness = fitness(ind)
    table.insert(parents, ind)
  end
  EvalPop(parents)
  if on_new_gen then
    on_new_gen(0, parents)
  end
  
  for gen = 1, MAX_GENERATIONS do
    local pool = selection_parent(parents, 2 * n)
    local offsprings = {}
    for i = 1, #pool, 2 do
      local parent1, parent2 = pool[i], pool[i + 1]
      local offspring = blend(parent1, parent2)
      offspring.fitness = fitness(offspring)
      table.insert(offsprings, offspring)
    end
    if selection_survival then
      if overlapping then
        for _, ind in ipairs(parents) do
          table.insert(offsprings, ind)
        end
      end
      parents = selection_survival(offsprings, m)
    else
      parents = offsprings
    end
    EvalPop(parents)
    if on_new_gen then
      on_new_gen(gen, parents)
    end
  end
  
  return gen, parents
end

local function GA(m, n, crossover, mutation, on_new_gen, survival_selection)
  return EV(m, n, GenRandomPop, on_new_gen, crossover, mutation, BinaryTournamentSelection, survival_selection)
end

local function EP(m, crossover, mutation, on_new_gen)
  return EV(m, m, GenRandomPop, on_new_gen, crossover, mutation, UniformStochasticSelection, TruncateSelection, "overlapping")
end

local function EA1(m, crossover, mutation, on_new_gen, n, survival_selection)
  return GA(m, n or m, crossover, mutation, on_new_gen, survival_selection)
end

local function EA1a(m, crossover, mutation, on_new_gen, n)
  return EV(m, n, GenRandomPop, on_new_gen, crossover, mutation, UniformDeterministicSelection, BinaryTournamentSelection)
end

local function EA2(m, crossover, mutation, on_new_gen)
  return EP(m, crossover, mutation, on_new_gen)
end

local function EA2a(m, crossover, mutation, on_new_gen, n)
  return EV(m, n, GenRandomPop, on_new_gen, crossover, mutation, UniformStochasticSelection, TruncateSelection, "overlapping")
end

local function EA2b(m, crossover, mutation, on_new_gen, n)
  return EV(m, n, GenRandomPop, on_new_gen, crossover, mutation, BinaryTournamentSelection, TruncateSelection, "overlapping")
end

local function EA1_Real(m, blend, on_new_gen)
  return EV_Real(m, m, on_new_gen, blend, BinaryTournamentSelection)
end

local function EA2_Real(m, blend, on_new_gen)
  return EV_Real(m, m, on_new_gen, blend, UniformStochasticSelection, TruncateSelection, "overlapping")
end

local s_Seeds = {}
for run = 1, Max(RUNS, 50) do
  s_Seeds[run] = math.random(1, 100000)
end

local function ExecRuns(EA, runs, points, crossover, mutation, points_fitness)
  for run = 1, runs do
    math.randomseed(s_Seeds[run])
    EA(POPULATION_SIZE, crossover, mutation, function(gen, pop)
      local rms = RMS_Hamming(pop)
      local entry = points[gen] or {x = gen, y = 0}
      points[gen] = entry
      entry.y = entry.y + rms
      if points_fitness then
        local entry = points_fitness[gen] or {x = gen, y = 0}
        points_fitness[gen] = entry
        entry.y = entry.y + pop.avg_fitness
      end
    end)
  end
end

local function DiversityCrossoverEA(EA, filename)
  local graphs = {name_x = "Generations", name_y = "Average Diversity of the Population", funcs = {}}
  local name_05 = string.format("%s: UniCros 0.5", EA == EA1 and "EA1" or "EA2")
  local name_03 = string.format("%s: UniCros 0.3", EA == EA1 and "EA1" or "EA2")
  local name_01 = string.format("%s: UniCros 0.1", EA == EA1 and "EA1" or "EA2")
  local name_no = string.format("%s: No Crossover", EA == EA1 and "EA1" or "EA2")
  local points_05 = {color = RGB_GREEN}
  local points_03 = {color = RGB_CYAN}
  local points_01 = {color = RGB_RED}
  local points_no = {color = RGB_ORANGE}
  graphs.funcs[name_05] = points_05
  graphs.funcs[name_03] = points_03
  graphs.funcs[name_01] = points_01
  graphs.funcs[name_no] = points_no
  
  ExecRuns(EA, RUNS, points_05, CrossoverUniformCreator(0.5))
  ExecRuns(EA, RUNS, points_03, CrossoverUniformCreator(0.3))
  ExecRuns(EA, RUNS, points_01, CrossoverUniformCreator(0.1))
  ExecRuns(EA, RUNS, points_no)
  NormalizeGraphs(graphs, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_x = 0, min_y = 0, div_x = 5, div_y = 4})
  bmp:WriteBMP(filename)
end

local function GetGens(EA, runs, crossover, mutation)
  local total_gens, total_fitness, total_success_crossovers = 0, 0, 0
  for run = 1, runs do
    math.randomseed(s_Seeds[run])
    local gens, pop = EA(POPULATION_SIZE, crossover, mutation, function(gen, pop)
      return IsPopConverged(pop) and "converged"
    end)
    total_gens = total_gens + gens
    total_fitness = total_fitness + pop[1].fitness
  end
  
  return total_gens / runs, total_fitness / runs
end

local function GenerationsToFixedPointAndAvgFitness(filename, filename2)
  local graphs = {name_x = "Uniform Crossover Bias Setting", name_y = "Average Generations to Fixed-Point Convergence", funcs = {}}
  local graphs2 = {name_x = "Uniform Crossover Bias Setting", name_y = "Average Population Fitness at Fixed-Point Convergence", funcs = {}}
  local name_ea1 = string.format("EA-1 Non-Overlap(m=n=%d) Binary Tournament", POPULATION_SIZE)
  local name_ea2 = string.format("EA-2 Overlap(m=n=%d) Truncation", POPULATION_SIZE)
  local points_ea1_gens, points_ea1_fitness = {color = RGB_GREEN}, {color = RGB_GREEN}
  local points_ea2_gens, points_ea2_fitness = {color = RGB_CYAN}, {color = RGB_CYAN}
  graphs.funcs[name_ea1], graphs2.funcs[name_ea1] = points_ea1_gens, points_ea1_fitness
  graphs.funcs[name_ea2], graphs2.funcs[name_ea2] = points_ea2_gens, points_ea2_fitness
  
  for rate = 0.0, 0.5, 0.05 do
    local crossover = (rate > 0.0) and CrossoverUniformCreator(rate)
    local gens1, fitness1 = GetGens(EA1, RUNS, crossover)
    local gens2, fitness2 = GetGens(EA2, RUNS, crossover)
    print(string.format("Crossover rate: %.2f, EA1 gens/fitness: %.2f/%.2f, EA2 gens/fitness: %.2f/%.2f", rate, gens1, fitness1, gens2, fitness2))
    table.insert(points_ea1_gens, {x = rate, y = gens1})
    table.insert(points_ea2_gens, {x = rate, y = gens2})
    table.insert(points_ea1_fitness, {x = rate, y = fitness1})
    table.insert(points_ea2_fitness, {x = rate, y = fitness2})
  end
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_y = true, skip_KP = true, min_x = 0, min_y = 0, div_x = 10, div_y = 8, max_y = 80})
  bmp:WriteBMP(filename)
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs2, {int_y = true, skip_KP = true, min_x = 0, min_y = 0, div_x = 10, div_y = 4, max_y = 34})
  bmp:WriteBMP(filename2)
end

local function ExecRunsSuccessRates(name, EA, runs, points, points_fchange, points_favg, crossover, mutation)
  for run = 1, runs do
    math.randomseed(s_Seeds[run])
    EA(POPULATION_SIZE, crossover, mutation, function(gen, pop)
      local entry = points[gen] or {x = gen, y = 0}
      points[gen] = entry
      entry.y = entry.y + pop.repro_successes / #pop
      local entry_fchange = points_fchange[gen] or {x = gen, y = 0}
      points_fchange[gen] = entry_fchange
      entry_fchange.y = entry_fchange.y - pop.avg_fitness_change
      local entry_favg = points_favg[gen] or {x = gen, y = 0}
      points_favg[gen] = entry_favg
      entry_favg.y = entry_favg.y + pop.avg_fitness
      if run % 100 == 0 and gen == MAX_GENERATIONS then
        print(string.format("%s, Run: %d/%d, Gen: %d, Success rate: %.2f, Avg Fitness Change: %.2f, Avg Fitness: %.2f", name, run, runs, gen, entry.y / run, entry_fchange.y / run, entry_favg.y / run))
      end
    end)
  end
end

local function SuccessRateEA(EA, filename, filename2, filename3)
  local graphs = {name_x = "Generations", name_y = "Average Number of Reproductive Successes", funcs = {}}
  local graphs2 = {name_x = "Generations", name_y = "Average Parent-Offspring Fitness Change", funcs = {}}
  local graphs3 = {name_x = "Generations", name_y = "Average Fitness of the Population", funcs = {}}
  local name_05 = string.format("%s: UniCros 0.5", EA == EA1 and "EA1" or "EA2")
  local name_04 = string.format("%s: UniCros 0.4", EA == EA1 and "EA1" or "EA2")
  local name_03 = string.format("%s: UniCros 0.3", EA == EA1 and "EA1" or "EA2")
  local name_02 = string.format("%s: UniCros 0.2", EA == EA1 and "EA1" or "EA2")
  local name_01 = string.format("%s: UniCros 0.1", EA == EA1 and "EA1" or "EA2")
  local points_05, points_05_fchange, points_05_favg = {color = RGB_GREEN}, {color = RGB_GREEN}, {color = RGB_GREEN}
  local points_04, points_04_fchange, points_04_favg = {color = RGB_CYAN}, {color = RGB_CYAN}, {color = RGB_CYAN}
  local points_03, points_03_fchange, points_03_favg = {color = RGB_RED}, {color = RGB_RED}, {color = RGB_RED}
  local points_02, points_02_fchange, points_02_favg = {color = RGB_ORANGE}, {color = RGB_ORANGE}, {color = RGB_ORANGE}
  local points_01, points_01_fchange, points_01_favg = {color = RGB_WHITE}, {color = RGB_WHITE}, {color = RGB_WHITE}
  graphs.funcs[name_05], graphs2.funcs[name_05], graphs3.funcs[name_05] = points_05, points_05_fchange, points_05_favg
  graphs.funcs[name_04], graphs2.funcs[name_04], graphs3.funcs[name_04] = points_04, points_04_fchange, points_04_favg
  graphs.funcs[name_03], graphs2.funcs[name_03], graphs3.funcs[name_03] = points_03, points_03_fchange, points_03_favg
  graphs.funcs[name_02], graphs2.funcs[name_02], graphs3.funcs[name_02] = points_02, points_02_fchange, points_02_favg
  graphs.funcs[name_01], graphs2.funcs[name_01], graphs3.funcs[name_01] = points_01, points_01_fchange, points_01_favg
  
  ExecRunsSuccessRates("Rate 0.5", EA, RUNS, points_05, points_05_fchange, points_05_favg, CrossoverUniformCreator(0.5))
  ExecRunsSuccessRates("Rate 0.4", EA, RUNS, points_04, points_04_fchange, points_04_favg, CrossoverUniformCreator(0.4))
  ExecRunsSuccessRates("Rate 0.3", EA, RUNS, points_03, points_03_fchange, points_03_favg, CrossoverUniformCreator(0.3))
  ExecRunsSuccessRates("Rate 0.2", EA, RUNS, points_02, points_02_fchange, points_02_favg, CrossoverUniformCreator(0.2))
  ExecRunsSuccessRates("Rate 0.1", EA, RUNS, points_01, points_01_fchange, points_01_favg, CrossoverUniformCreator(0.1))
  NormalizeGraphs(graphs, RUNS)
  NormalizeGraphs(graphs2, RUNS)
  NormalizeGraphs(graphs3, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_x = 0, min_y = 0, div_x = 8, div_y = 6, max_y = 0.6})
  bmp:WriteBMP(filename)
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs2, {int_x = true, skip_KP = true, min_x = 0, min_y = 0, div_x = 8, div_y = 6, max_y = 0.6})
  bmp:WriteBMP(filename2)
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs3, {int_x = true, skip_KP = true, min_x = 0, min_y = 0, div_x = 8, div_y = 6, max_y = 0.6})
  bmp:WriteBMP(filename3)
end

local function ExecRuns_Real(EA, runs, points, blend)
  for run = 1, runs do
    math.randomseed(s_Seeds[run])
    EA(POPULATION_SIZE, blend, function(gen, pop)
      local entry = points[gen] or {x = gen, y = 0}
      points[gen] = entry
      entry.y = entry.y + pop.avg_fitness
    end)
  end
end

local function BlendingFactorEA(EA, filename)
  local graphs = {name_x = "Generations", name_y = string.format("%s Real Genome Average Population Fitness", EA == EA1_Real and "EA1" or "EA2"), funcs = {}}
  local name_05 = string.format("Blending 0.5")
  local name_03 = string.format("Blending 0.3")
  local name_01 = string.format("Blending 0.1")
  local points_05 = {color = RGB_GREEN}
  local points_03 = {color = RGB_CYAN}
  local points_01 = {color = RGB_RED}
  graphs.funcs[name_05] = points_05
  graphs.funcs[name_03] = points_03
  graphs.funcs[name_01] = points_01
  
  ExecRuns_Real(EA, RUNS, points_05, CrossoverBlendCreator(0.5))
  ExecRuns_Real(EA, RUNS, points_03, CrossoverBlendCreator(0.3))
  ExecRuns_Real(EA, RUNS, points_01, CrossoverBlendCreator(0.1))
  NormalizeGraphs(graphs, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_x = 0, min_y = 0, div_x = 5, div_y = 4})
  bmp:WriteBMP(filename)
end

local function DiversityMutationEA(EA, diversity_rates, fitness_rates, filename, filename2)
  local graphs = {name_x = "Generations", name_y = string.format("%s: Average Diversity of the Population", EA == EA1 and "EA1" or "EA2"), funcs = {}}
  local graphs2 = {name_x = "Generations", name_y = string.format("%s: Average Fitness of the Population", EA == EA1 and "EA1" or "EA2"), funcs = {}}
  local name_div1 = string.format("Mutation %.2f", diversity_rates[1])
  local name_div2 = string.format("Mutation %.2f", diversity_rates[2])
  local name_div3 = string.format("Mutation %.2f", diversity_rates[3])
  local name_div4 = string.format("Mutation %.2f", diversity_rates[4])
  local name_div_no = string.format("No Mutation")
  local name_fit1 = string.format("Mutation %.2f", fitness_rates[1])
  local name_fit2 = string.format("Mutation %.2f", fitness_rates[2])
  local name_fit3 = string.format("Mutation %.2f", fitness_rates[3])
  local name_fit4 = fitness_rates[4] and string.format("Mutation %.2f", fitness_rates[4])
  local name_fit_no = string.format("No Mutation")
  local points_div1, points_fit_1 = {color = RGB_GREEN}, {color = RGB_GREEN}
  local points_div2, points_fit_2 = {color = RGB_CYAN}, {color = RGB_CYAN}
  local points_div3, points_fit_3 = {color = RGB_RED}, {color = RGB_RED}
  local points_div4, points_fit_4 = {color = RGB_ORANGE}, {color = RGB_ORANGE}
  local points_div_no, points_fit_no = {color = RGB_WHITE}, {color = RGB_WHITE}
  graphs.funcs[name_div1], graphs2.funcs[name_fit1] = points_div1, points_fit_1
  graphs.funcs[name_div2], graphs2.funcs[name_fit2] = points_div2, points_fit_2
  graphs.funcs[name_div3], graphs2.funcs[name_fit3] = points_div3, points_fit_3
  graphs.funcs[name_div4] = points_div4
  if name_fit4 then
    graphs2.funcs[name_fit4] = points_fit_4
  end
  graphs.funcs[name_div_no], graphs2.funcs[name_fit_no] = points_div_no, points_fit_no
  
  ExecRuns(EA, RUNS, points_div1, nil, MutationCreator(diversity_rates[1]), points_fit_1)
  ExecRuns(EA, RUNS, points_div2, nil, MutationCreator(diversity_rates[2]), points_fit_2)
  ExecRuns(EA, RUNS, points_div3, nil, MutationCreator(diversity_rates[3]), points_fit_3)
  ExecRuns(EA, RUNS, points_div4, nil, MutationCreator(diversity_rates[4]), points_fit_4)
  ExecRuns(EA, RUNS, points_div_no, nil, nil, points_fit_no)
  NormalizeGraphs(graphs, RUNS)
  NormalizeGraphs(graphs2, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_x = 0, min_y = 0, div_x = 8, div_y = 4})
  bmp:WriteBMP(filename)
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs2, {int_x = true, skip_KP = true, min_x = 0, min_y = 0, div_x = 8, div_y = 8})
  bmp:WriteBMP(filename2)
end

local function ExecSuccessRates(EA, runs, points, mutation)
  for run = 1, runs do
    math.randomseed(s_Seeds[run])
    EA(POPULATION_SIZE, nil, mutation, function(gen, pop)
      if gen == 0 then return end
      local entry = points[gen] or {x = gen, y = 0}
      points[gen] = entry
      entry.y = entry.y + pop.repro_successes / #pop
    end)
  end
end

local function SuccessRatesMutationEA(EA, filename)
  local graphs = {name_x = "Generations", name_y = string.format("%s: Average Number of Reproduction Successes", EA == EA1 and "EA1" or "EA2"), funcs = {}}
  local name_05 = string.format("Mutation 0.5")
  local name_03 = string.format("Mutation 0.3")
  local name_02 = string.format("Mutation 0.2")
  local name_01 = string.format("Mutation 0.1")
  local name_005 = string.format("Mutation 0.05")
  local points_05 = {color = RGB_GREEN}
  local points_03 = {color = RGB_CYAN}
  local points_02 = {color = RGB_RED}
  local points_01 = {color = RGB_ORANGE}
  local points_005 = {color = RGB_WHITE}
  graphs.funcs[name_05] = points_05
  graphs.funcs[name_03] = points_03
  graphs.funcs[name_02] = points_02
  graphs.funcs[name_01] = points_01
  graphs.funcs[name_005] = points_005
  
  ExecSuccessRates(EA, RUNS, points_05, MutationCreator(0.5))
  ExecSuccessRates(EA, RUNS, points_03, MutationCreator(0.3))
  ExecSuccessRates(EA, RUNS, points_02, MutationCreator(0.2))
  ExecSuccessRates(EA, RUNS, points_01, MutationCreator(0.1))
  ExecSuccessRates(EA, RUNS, points_005, MutationCreator(0.05))
  NormalizeGraphs(graphs, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_x = 0, min_y = 0, div_x = 8, div_y = 4})
  bmp:WriteBMP(filename)
end

local function ExecAvgFitness(EA, runs, points, crossover, mutation, m, n)
  if m and n then
    MAX_GENERATIONS = MAX_BIRTHS / (m + n)
  end
  for run = 1, runs do
    math.randomseed(s_Seeds[run])
    EA(m or POPULATION_SIZE, crossover, mutation, function(gen, pop)
      local entry = points[gen] or {x = (m and n) and (gen + 1) * (m + n) or gen, y = 0}
      points[gen] = entry
      entry.y = entry.y + pop.avg_fitness
    end, m ~= n and n, m ~= n and UniformStochasticSelection)
  end
end

local function CrossoverMutationEA(EA, filename)
  local graphs = {name_x = "Generations", name_y = string.format("%s: Average Fitness of the Population", EA == EA1 and "EA1" or "EA2"), funcs = {}}
  local name_crossover = string.format("No Mutation, Crossover 0.2")
  local name_mutation = string.format("Mutation 0.01, No Crossover")
  local name_crossover_mutation = string.format("Mutation 0.001, Crossover 0.1")
  local points_crossover = {color = RGB_GREEN}
  local points_mutation = {color = RGB_CYAN}
  local points_crossover_mutation = {color = RGB_WHITE}
  graphs.funcs[name_crossover] = points_crossover
  graphs.funcs[name_mutation] = points_mutation
  graphs.funcs[name_crossover_mutation] = points_crossover_mutation
  
  ExecAvgFitness(EA, RUNS, points_crossover, CrossoverUniformCreator(EA == EA1 and 0.2 and 0.5))
  ExecAvgFitness(EA, RUNS, points_mutation, nil, MutationCreator(EA == EA1 and 0.01 or 0.05))
  ExecAvgFitness(EA, RUNS, points_crossover_mutation, CrossoverUniformCreator(EA == EA1 and 0.1 or 0.4), MutationCreator(EA == EA1 and 0.001 or 0.01))
  NormalizeGraphs(graphs, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_x = 0, min_y = 0, div_x = 8, div_y = 4})
  bmp:WriteBMP(filename)
end

local function CompareEA(EA, crossover_rate, mutation_rate, pop_sizes, filename)
  local graphs = {name_x = string.format("Births: Uniform Crossover = %.2f, Mutation = %.3f", crossover_rate, mutation_rate), name_y = string.format("%s: Average Fitness of the Population", EA == EA1 and "EA1" or "EA2"), funcs = {}}
  local colors = {RGB_GREEN, RGB_CYAN, RGB_WHITE, RGB_ORANGE, RGB_RED, RGB_BLUE}
  local names, points = {}, {}
  for i, sizes in ipairs(pop_sizes) do
    names[i] = string.format("m=%d,n=%d", sizes[1], sizes[2])
    points[i] = {color = colors[i]}
    graphs.funcs[names[i]] = points[i]
  end
  
  local crossover = CrossoverUniformCreator(crossover_rate)
  local mutation = MutationCreator(mutation_rate)
  for i, sizes in ipairs(pop_sizes) do
    ExecAvgFitness(EA, RUNS, points[i], crossover, mutation, sizes[1], sizes[2])
  end
  MAX_GENERATIONS = 80
  NormalizeGraphs(graphs, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_x = 0, min_y = 0, div_x = 8, div_y = 4})
  bmp:WriteBMP(filename)
end

local function CompareEAs(crossover_rate, mutation_rate, pop_sizes, filename)
  local graphs = {name_x = string.format("Births: Uniform Crossover = %.2f, Mutation = %.3f", crossover_rate, mutation_rate), name_y = string.format("Average Fitness of the Population"), funcs = {}}
  local colors = {RGB_GREEN, RGB_CYAN, RGB_WHITE, RGB_ORANGE, RGB_RED, RGB_BLUE}
  local names, points = {}, {}
  for i, sizes in ipairs(pop_sizes) do
    names[i] = string.format("%s: m=%d,n=%d", sizes[4], sizes[1], sizes[2])
    points[i] = {color = colors[i]}
    graphs.funcs[names[i]] = points[i]
  end
  
  local crossover = CrossoverUniformCreator(crossover_rate)
  local mutation = MutationCreator(mutation_rate)
  for i, sizes in ipairs(pop_sizes) do
    ExecAvgFitness(sizes[3], RUNS, points[i], crossover, mutation, sizes[1], sizes[2])
  end
  MAX_GENERATIONS = 80
  NormalizeGraphs(graphs, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_x = 0, min_y = 0, div_x = 8, div_y = 4})
  bmp:WriteBMP(filename)
end

DiversityCrossoverEA(EA1, IMAGE_FILENAME_FIGURE_6_30)
DiversityCrossoverEA(EA2, IMAGE_FILENAME_FIGURE_6_31)
GenerationsToFixedPointAndAvgFitness(IMAGE_FILENAME_FIGURE_6_32, IMAGE_FILENAME_FIGURE_6_33)
SuccessRateEA(EA1, IMAGE_FILENAME_FIGURE_6_34, IMAGE_FILENAME_FIGURE_6_35, IMAGE_FILENAME_FIGURE_6_36)
SuccessRateEA(EA2, IMAGE_FILENAME_FIGURE_6_37, IMAGE_FILENAME_FIGURE_6_38, IMAGE_FILENAME_FIGURE_6_39)
BlendingFactorEA(EA1_Real, IMAGE_FILENAME_FIGURE_6_40)
BlendingFactorEA(EA2_Real, IMAGE_FILENAME_FIGURE_6_41)
DiversityMutationEA(EA1, {0.3, 0.2, 0.1, 0.05}, {0.3, 0.1, 0.05}, IMAGE_FILENAME_FIGURE_6_42, IMAGE_FILENAME_FIGURE_6_44)
DiversityMutationEA(EA2, {0.5, 0.3, 0.1, 0.05}, {0.3, 0.2, 0.1, 0.05}, IMAGE_FILENAME_FIGURE_6_43, IMAGE_FILENAME_FIGURE_6_45)
SuccessRatesMutationEA(EA1, IMAGE_FILENAME_FIGURE_6_46)
SuccessRatesMutationEA(EA2, IMAGE_FILENAME_FIGURE_6_47)

MAX_GENERATIONS = 120
GENOTYPE_LENGTH = 4 * GENOME_BITS

CrossoverMutationEA(EA1, IMAGE_FILENAME_FIGURE_6_48)
CrossoverMutationEA(EA2, IMAGE_FILENAME_FIGURE_6_49)
CompareEA(EA1, 0.1, 0.001, {{10, 10}, {10, 20}, {20, 20}, {20, 40}, {40, 40}}, IMAGE_FILENAME_FIGURE_6_50)
CompareEA(EA2, 0.4, 0.01, {{10, 10}, {10, 20}, {10, 30}, {20, 20}, {20, 40}, {40, 40}}, IMAGE_FILENAME_FIGURE_6_52)
CompareEAs(0.1, 0.001, {{20, 20, EA1, "EA1"}, {20, 40, EA1, "EA1"}, {20, 40, EA1a, "EA1a"}}, IMAGE_FILENAME_FIGURE_6_51)
CompareEAs(0.4, 0.01, {{40, 40, EA2, "EA2"}, {40, 20, EA2a, "EA2a"}, {40, 1, EA2a, "EA2a"}, {40, 1, EA2b, "EA2b"}}, IMAGE_FILENAME_FIGURE_6_53)
