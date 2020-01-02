dofile("Statistics.lua")
dofile("Bitmap.lua")
dofile("Graphics.lua")
dofile("EC_Common.lua")
dofile("GA_Common.lua")

local RUNS                          = 400
local MAX_GENERATIONS               = 1000
local MAX_BIRTHS                    = 1000

local MIN_X                         = 0.0
local MAX_X                         = 4.0
local LANDSCAPE_A_PEAK_X            = 3.5
local LANDSCAPE_B_PEAK_X            = 0.5
local MAX_FITNESS                   = 10.0

local GENOME_BITS                   = 16
local GENOME_MAX_VALUE              = math.pow(2, GENOME_BITS) - 1
local GENOME_PRECISION              = (MAX_X - MIN_X) / GENOME_MAX_VALUE

local IMAGE_WIDTH                   = 1000
local IMAGE_HEIGHT                  = 1000
local IMAGE_FILENAME_FIGURE_7_1     = "EC/Figure7_1.bmp"
local IMAGE_FILENAME_FIGURE_7_2     = "EC/Figure7_2.bmp"
local IMAGE_FILENAME_FIGURE_7_3     = "EC/Figure7_3.bmp"

local function RealToBinary(x)
  return math.floor(MAX_X * (x - MIN_X) / (MAX_X - MIN_X))
end

local function BinaryToReal(x)
  return MIN_X + (MAX_X - MIN_X) * x / GENOME_MAX_VALUE
end

local function LandscapeA(x)
  if x < LANDSCAPE_A_PEAK_X then
    return x * MAX_FITNESS / LANDSCAPE_A_PEAK_X
  else
    return (MAX_X - x) * MAX_FITNESS / (MAX_X - LANDSCAPE_A_PEAK_X)
  end
end

local function LandscapeB(x)
  if x < LANDSCAPE_B_PEAK_X then
    return x * MAX_FITNESS / LANDSCAPE_B_PEAK_X
  else
    return (MAX_X - x) * MAX_FITNESS / (MAX_X - LANDSCAPE_B_PEAK_X)
  end
end

local s_Landscape = LandscapeA

local function Fitness(genotype)
  return s_Landscape(BinaryToReal(genotype[1]))
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
  pop.min_fitness, pop.max_fitness, pop.avg_fitness = min_fitness or 0, max_fitness or 0, total_fitness / #pop
  pop.avg_fitness_change, pop.repro_successes = total_fitness_change / #pop, repro_successes
end

local function GenRandomPop(size, on_new_birth)
  local pop = {}
  while #pop < size do
    local bitstring = GenRandomBitstring(GENOME_BITS)
    local genotype = PackBitstring(bitstring)
    local ind = {genotype = genotype, fitness = Fitness(genotype)}
    table.insert(pop, ind)
    if on_new_birth then
      on_new_birth(0, pop, #pop, ind, {})
    end
  end
  
  return pop
end

local function UniformStochasticSelection(parents, size)
  local offsprings = {}
  for i = 1, size do
    local parent = parents[math.random(1, #parents)]
    table.insert(offsprings, parent)
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

local function GetFitnessProportionalProb(parents)
  local prob, total_prob = {}, 0.0
  for i, ind in ipairs(parents) do
    total_prob = total_prob + ind.fitness
    prob[i] = total_prob
  end
  prob.total = total_prob
  
  return prob
end

local function FitnessProportionalSelection(parents, size)
  local prob = GetFitnessProportionalProb(parents)
  local offsprings = {}
  while #offsprings < size do
    local rand = math.random() * prob.total
    local i = table.binary_search(prob, rand)
    table.insert(offsprings, parents[i])
  end
  
  return offsprings
end

local function Crossover1Point(parent1, parent2)
  local offspring1, offspring2 = table.copy(parent1), table.copy(parent2)
  Binary1PointCrossover(offspring1, offspring2, math.random(1, parent1.bits))
  
  return offspring1, offspring2
end

local function MutationCreator(rate)
  return function(genotype)
    return BinaryMutation(genotype, rate)
  end
end

local function EV(m, n, on_new_birth, on_new_gen, crossover, mutation, selection_parent, selection_survival, overlapping)
  local parents = GenRandomPop(m, on_new_birth)
  EvalPop(parents)
  if on_new_gen then
    local result = on_new_gen(0, parents)
    if result == "converged" then
      return 0, parents
    end
  end
  local gen, births = 1, #parents
  while gen <= MAX_GENERATIONS and births < MAX_BIRTHS do
    local pool = selection_parent(parents, n)
    local offsprings = {}
    if crossover then
      for i = 1, #pool, 2 do
        local parent1, parent2 = pool[i], (i < #pool) and pool[i + 1] or pool[i]
        local offspring1, offspring2 = crossover(parent1.genotype, parent2.genotype)
        local new_ind = {genotype = offspring1, parent_fitness = parent1.fitness}
        if mutation then
          mutation(new_ind.genotype)
        end
        new_ind.fitness = Fitness(new_ind.genotype)
        table.insert(offsprings, new_ind)
        births = births + 1
        if on_new_birth then
          on_new_birth(gen, parents, births, new_ind, offsprings)
        end
        if offspring2 then
          local new_ind = {genotype = offspring2, parent_fitness = parent1.fitness}
          if mutation then
            mutation(new_ind.genotype)
          end
          new_ind.fitness = Fitness(new_ind.genotype)
          table.insert(offsprings, new_ind)
          births = births + 1
          if on_new_birth then
            on_new_birth(gen, parents, births, new_ind, offsprings)
          end
        end
      end
    else
      offsprings = {}
      for _, ind in ipairs(pool) do
        table.insert(offsprings, ind)
        births = births + 1
        if on_new_birth then
          on_new_birth(gen, parents, births, ind, offsprings)
        end
      end
    end
    if not crossover and mutation then
      for i, ind in ipairs(offsprings) do
        local mutant = table.copy(ind, "deep")
        mutant.parent_fitness = ind.fitness
        mutation(mutant.genotype)
        mutant.fitness = Fitness(mutant.genotype)
        offsprings[i] = mutant
        births = births + 1
        if on_new_birth then
          on_new_birth(gen, parents, births, mutant, offsprings)
        end
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

local function GA(m, n, crossover, mutation, on_new_birth, on_new_gen)
  return EV(m, n, on_new_birth, on_new_gen, crossover, mutation, FitnessProportionalSelection)
end

local function ES(m, n, crossover, mutation, on_new_birth, on_new_gen)
  return EV(m, n, on_new_birth, on_new_gen, crossover, mutation, UniformStochasticSelection, TruncateSelection, "overlapping")
end

local s_Seeds = {}
for run = 1, RUNS do
  s_Seeds[run] = math.random(1, 100000)
end

local function ExecRuns(EA, points, m, n, crossover, mutation, landscape_change_check)
  for run = 1, RUNS do
    math.randomseed(s_Seeds[run])
    s_Landscape = LandscapeA
    EA(m, n, crossover, mutation, function(gen, pop, births, new_ind, offsprings)
      if landscape_change_check(births) then
        s_Landscape = (s_Landscape == LandscapeB) and LandscapeA or LandscapeB
        for _, ind in ipairs(pop) do
          ind.fitness = Fitness(ind.genotype)
        end
        for _, ind in ipairs(offsprings) do
          ind.fitness = Fitness(ind.genotype)
        end
      end
      EvalPop(pop)
      EvalPop(offsprings)
      local best = (pop.max_fitness > offsprings.max_fitness) and pop.max_fitness or offsprings.max_fitness
      local entry = points[births] or {x = births, y = 0}
      points[births] = entry
      entry.y = entry.y + best
    end)
  end
end

local function BestInPop(filename, landscape_change_check)
  local graphs = {name_x = "Trials", name_y = "Fitness: Best in Population", funcs = {}}
  local name_ga = string.format("GA")
  local name_es = string.format("ES")
  local points_ga = {color = RGB_GREEN}
  local points_es = {color = RGB_CYAN}
  graphs.funcs[name_ga] = points_ga
  graphs.funcs[name_es] = points_es
  
  ExecRuns(ES, points_es, 1, 10, Crossover1Point, MutationCreator(1 / GENOME_BITS), landscape_change_check)
  ExecRuns(GA, points_ga, 10, 10, Crossover1Point, MutationCreator(1 / GENOME_BITS), landscape_change_check)
  NormalizeGraphs(graphs, RUNS)
  
  local bmp = Bitmap.new(IMAGE_WIDTH, IMAGE_HEIGHT, RGB_BLACK)
  DrawGraphs(bmp, graphs, {int_x = true, skip_KP = true, min_x = 0, min_y = 0, div_x = 10, div_y = 5, max_y = MAX_FITNESS})
  bmp:WriteBMP(filename)
end

BestInPop(IMAGE_FILENAME_FIGURE_7_1, function(births) return births == 500 end)
BestInPop(IMAGE_FILENAME_FIGURE_7_2, function(births) return births > 1 and births % 300 == 0 end)
BestInPop(IMAGE_FILENAME_FIGURE_7_3, function(births) return births > 1 and births % 100 == 0 end)
